#!/usr/bin/env python
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import concurrent.futures
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import matplotlib.pyplot as plt
import pyBigWig
import seaborn as sns
import numpy as np
import pandas as pd
import os
import glob
import pybedtools
import gffutils
from Bio import SeqIO
from tqdm import trange
import sys
from collections import defaultdict, OrderedDict
from pandas.errors import EmptyDataError
pd.options.mode.chained_assignment = None  # default='warn'
pd.set_option('display.max_rows', 500)

class Density:
    def values(self, chrom, start, end, strand):
        return 0
    
class ReadDensity(Density):
    """
    ReadDensity class
    Attributes:
        self.pos(positive *.bw file)
        self.neg(negative *.bw file)
    """

    def __init__(self, pos, neg, name=None):
        try:
            self.pos = pyBigWig.open(pos)
            self.neg = pyBigWig.open(neg)
            self.name = name if name is not None else pos.replace(
                'fwd', '*'
            ).replace(
                'rev', '*'
            )

        except Exception as e:
            print("couldn't open the bigwig files!")
            print(e)

    def values(self, chrom, start, end, strand):
        """
        Parameters
        ----------
        chrom : basestring
            (eg. chr1)
        start : int
            0-based start (first position in chromosome is 0)
        end : int
            1-based end (last position is not included)
        strand : str
            either '+' or '-'
        Returns
        -------
        densites : list
            values corresponding to density over specified positions.
        """

        try:
            if strand == "+":
                return list(pd.Series(self.pos.values(chrom, start, end)).fillna(0))
            elif strand == "-":
                return list(pd.Series(self.neg.values(chrom, start, end)).fillna(0))
            else:
                print("Strand neither + or -")
                return 1
        except RuntimeError:
            # usually occurs when no chromosome exists in the bigwig file
            return list(pd.Series([np.NaN] * abs(start - end)).fillna(0))

def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

def get_c_positions_and_coverage_in_window(bedtool, rdd, genome_fa):
    """
    Given a bedtool (BedTool) of regions, return the positions (0-based) of C's (or G's on neg strand).
    """
    d = {}
    pos = bedtool.to_dataframe().iloc[0]
    chrom = pos.chrom
    start = pos.start
    end = pos.end
    strand = pos.strand
    sequences = bedtool.sequence(fi=genome_fa, s=False, name=True) 
    with open(sequences.seqfn) as f:
        for record in SeqIO.parse(f, "fasta"):
            if strand == '+':
                relpos = find(record.seq.upper(), 'C')
            elif strand == '-':
                relpos = find(record.seq.upper(), 'G')
            # print(record.seq)
            # print(relpos)
    abspos = ["{}:{}".format(chrom, start + p) for p in relpos]
    coverage = rdd.values(chrom=chrom, start=start, end=end, strand=strand)
    coverage = [np.abs(c) for c in coverage]  # doesn't matter for how we're making bigwigs, but just to be sure. 
    c_coverage = [coverage[p] for p in relpos]
    for p, c in zip(abspos, c_coverage):
        d[p] = c
    return d


def get_total_c_coverage(row, rdd, genome_fa):
    bedtool = pybedtools.BedTool(
        [pybedtools.create_interval_from_list(
            [row['chrom'], row['start'], row['end'], "window", "0", row['strand']]
        )]
    )
    return sum_all_c_coverage_in_window(bedtool, rdd, genome_fa)


def sum_all_c_coverage_in_window(bedtool, rdd, genome_fa):
    all_coverage = 0
    c_positions = get_c_positions_and_coverage_in_window(bedtool, rdd, genome_fa)
    for pos, cov in c_positions.items():
        all_coverage += cov
    return int(all_coverage)


## helper func
def get_gene_positions_dict(db):
    """
    Returns the genic positions for each gene in a gencode annotation database.
    ie. {"gene_id":{'chrom':chr1,'start':0,'end':100,'strand':"+"}}
    """
    gene_positions = defaultdict()
    genes = db.features_of_type('gene')
    for gene in genes:
        gene_id = gene.attributes['gene_id'][0] if type(gene.attributes['gene_id']) == list else gene.attributes['gene_id']
        gene_positions[gene_id] = {'chrom':gene.seqid, 'start':gene.start-1, 'end':gene.end, 'strand':gene.strand}
    return gene_positions


def read_and_filter_editing_sites(sites_file, conf, v=False):
    annotated_headers = [
        'chrom','start','end','conf','edit_frac','strand',
        'gene_id','gene_name','region','annotation_string'
    ]
    sites = pd.read_csv(sites_file, names=annotated_headers, sep='\t')
    if v:
        print("before confidence filtering: ",sites.shape[0])

    sites = sites[sites['conf'] >=conf]
    # sites = sites[['chrom','start','end','edit_frac','gene_id','region','strand']]
    if v:
        print("after confidence filtering: ",sites.shape[0])
    return sites

## pilot window func

def create_chrom_sizes_dict(chrom_sizes):
    """
    Creates a chrom sizes dictionary. Useful for identifying 
    chromosomal boundaries in which we cannot exceed.
    """
    chrom_sizes_dict = {}
    with open(chrom_sizes, 'r') as f:
        for line in f:
            chrom, size = line.strip('\n').split('\t')
            chrom_sizes_dict[chrom] = int(size)
    return chrom_sizes_dict

def create_window_interval(center_site, flank, chrom_sizes_dict):
    """
    Given a position (center_site), and a flanking distance,
    return a window interval.
    window size will be flank + 1 + flank (ie. 2 + 1 + 2 = 5 if flank = 2)
    """
    window_start = center_site.start - flank
    window_start = 0 if window_start < 0 else window_start
    
    window_end = center_site.start + flank + 1
    window_end = chrom_sizes_dict[center_site.chrom] if window_end > chrom_sizes_dict[center_site.chrom] else window_end
    
    return pybedtools.create_interval_from_list(
        [
            center_site.chrom, 
            str(window_start), 
            str(window_end), 
            center_site.name, 
            center_site.score, 
            center_site.strand
        ]
    )

def create_window_intervals(sites, flank, chrom_sizes_dict):
    windows = []
    for site in sites:
        windows.append(create_window_interval(site, flank, chrom_sizes_dict))
    return pybedtools.BedTool(windows)
    

def get_regions_to_count(cds_file, three_prime_utr_file, five_prime_utr_file):
    # (2) Reads in a pre-parsed set of regions files (expect merged/nonoverlapping coordinates of exons corresponding to CDS/UTR)
    # This reads in a dataframe but indexes by 'name' because this is 100x faster to filter than bedtools.filter()
    regions_to_count = pd.DataFrame()
    
    if cds_file is not None:
        cds = pd.read_csv(
            cds_file, sep='\t', names=[
                'chrom','start','end','name','score','strand'
            ], 
            index_col='name'
        )
        regions_to_count = pd.concat([regions_to_count, cds])
    else: 
        cds = None
    if three_prime_utr_file is not None:
        three_prime_utr = pd.read_csv(
            three_prime_utr_file, sep='\t', names=[
                'chrom','start','end','name','score','strand'
            ], 
            index_col='name'
        )
        regions_to_count = pd.concat([regions_to_count, three_prime_utr])
    else:
        three_prime_utr = None
    if five_prime_utr_file is not None:
        five_prime_utr = pd.read_csv(
            five_prime_utr_file, sep='\t', names=[
                'chrom','start','end','name','score','strand'
            ], 
            index_col='name'
        )
        regions_to_count = pd.concat([regions_to_count, five_prime_utr])
    else:
        five_prime_utr = None
    return regions_to_count

## pilot scoring func


def score_edits_exons(
    annotated_edits_file, 
    bg_edits_file,
    output_file, 
    output_file_summed,
    conf, 
    genome_fa, 
    chrom_sizes_file, 
    rdd,
    cds_file,
    three_prime_utr_file,
    five_prime_utr_file
):
    """
    1. Reads and filters our (annotated) editing site (fg). The "name" (edit_frac) MUST contain the edited,totalcov for each site.
    2. Reads in a pre-parsed set of regions files (expect merged/nonoverlapping coordinates of exons corresponding to CDS/UTR)
        A. Concatenates this list so we don't necessarily care which regions are CDS vs UTR etc. just that these regions 
           will define where we calculate coverages on.
    ##### TODO ##### 3. Filter out bg_edits_file from fg_edits_file. 
    4. For every gene, generate a (non-overlapping) BedTool containing all the regions within that gene we want to count.
    5. Subsets our (annotated) editing list (fg) to get only the edits across one gene, for every gene. If a gene has no edit sites, pass. 
        A. For this step, we're relying on what's been annotated by annotator. So we are only counting edits that are unambiguously assigned
           (edits at a given position that overlaps multiple genes in the same region are not counted). 
    6. Intersect with all edits from (3) to collect all edits that exist within the window.
    7. Add up all the edited-reads and total-reads across edited sites and calculate the "edit/editedc" fraction.
    8. Calculate the coverage across all C's in each window
    9. Since we actually want per-gene coverage info instead of per-exon coverage, I'll sum these up here.
    10. Save the summed coverage info per-gene
    11. Save the coverage info per-exon
    """
    chrom_sizes_dict = create_chrom_sizes_dict(chrom_sizes_file)
    # (1) Reads and filters our (annotated) editing site (fg). The "name" (edit_frac) MUST contain the edited,totalcov for each site.
    fg = read_and_filter_editing_sites(annotated_edits_file, conf)
    
    all_scores_df = pd.DataFrame(
        columns = [
            'chrom','start','end','name','score',
            'strand','edit_coverage','editable_coverage',
            'edited_over_edited_c','all_c_coverage','edited_over_all_c'
        ]
    )
    all_scores = []
    regions_to_count = get_regions_to_count(
        cds_file=cds_file,
        three_prime_utr_file=three_prime_utr_file,
        five_prime_utr_file=five_prime_utr_file
    )
    fg.set_index('gene_id', inplace=True)
          
    progress = trange(len(set(fg.index)))
    for gene_id in set(fg.index): # set(fg['gene_id']):
        print("Trying: {}".format(gene_id))
        try:
            # (4) For every gene, generate a (non-overlapping) BedTool containing all the regions within that gene we want to count.
            bedtool = pybedtools.BedTool.from_dataframe(
                regions_to_count.loc[[gene_id]].reset_index()[
                    ['chrom','start','end','name','score','strand']
                ]
            ).sort().merge(c=(4,5,6), o=('distinct','sum','distinct'))
            # (5) Subsets our (annotated) editing list (fg) to get only the edits across one gene, for every gene. If a gene has no edit sites, pass.
            
            fg_sites_in_region = fg.loc[[gene_id]].reset_index() # fg[fg['gene_id']==gene_id]

            if fg_sites_in_region.shape[0] >= 1:
                # thickStart = edited #
                fg_sites_in_region.loc[:,'thickStart'] = fg_sites_in_region['edit_frac'].apply(
                    lambda x: int(x.split(',')[0])
                )
                # thickEnd = total coverage # 
                fg_sites_in_region.loc[:,'thickEnd'] = fg_sites_in_region['edit_frac'].apply(
                    lambda x: int(x.split(',')[1])
                )
                fg_sites_in_region.loc[:,'name'] = fg_sites_in_region.loc[:,'gene_id'] + \
                    "|" + fg_sites_in_region.loc[:,'region']
                fg_sites_in_region = fg_sites_in_region[
                    ['chrom','start','end','name','conf','strand','thickStart','thickEnd']
                ]
                # (4) Filter out bg_edits_file from fg_edits_file. 
                fg_prefiltered_sites_bedtool = pybedtools.BedTool.from_dataframe(fg_sites_in_region) 
                if bg_edits_file is not None:
                    bg_sites_bedtool = pybedtools.BedTool(bg_edits_file)
                    fg_sites_bedtool = fg_prefiltered_sites_bedtool.sort().intersect(bg_sites_bedtool.sort(), s=True, v=True)
                else:
                    fg_sites_bedtool = fg_prefiltered_sites_bedtool.sort()
                
                if len(fg_sites_bedtool) > 0:
                    # (6) Intersect with all edits from (3) to collect all edits that exist within the window.
                    intersected_edits = bedtool.intersect(
                        fg_sites_bedtool, s=True, wa=True, loj=True
                    ).to_dataframe(
                        names=[
                            'chrom','start','end','name','score','strand',
                            'edit_chrom','edit_start','edit_end','edit_name',
                            'edit_score','edit_strand','edit_coverage','editable_coverage'
                        ]
                    )
                    # We can to calculate across all exons, regardless of whether or not we found an edit there. 
                    try:
                        intersected_edits['edit_coverage'].replace({'.':0}, inplace=True)
                        intersected_edits['editable_coverage'].replace({'.':0}, inplace=True)
                    except TypeError:  # it will complain that dtype of edit_coverage and editable_coverage isn't a string
                        pass

                    # pybedtools does some weird things with mixed dtypes, so let's make sure we convert to INT here (so 1+1 = 2, not "11")
                    intersected_edits['edit_coverage'] = intersected_edits['edit_coverage'].astype(int)
                    intersected_edits['editable_coverage'] = intersected_edits['editable_coverage'].astype(int)

                    if intersected_edits.shape[0] > 0:
                        # (7) Add up all the edited-reads and total-reads across edited sites and calculate the "edit/editedc" fraction.

                        # blockCount is the "number of reads supporting an edit site"
                        summed_edits = pd.DataFrame(
                            intersected_edits.groupby(
                                ['chrom','start','end','name','score','strand']
                            )['edit_coverage'].sum()
                        ).reset_index()
                        # editable_coverage (blockSizes) is the "total number of reads at the edited site"
                        summed_total_coverage = pd.DataFrame(
                            intersected_edits.groupby(
                                ['chrom','start','end','name','score','strand']
                            )['editable_coverage'].sum()).reset_index()

                        # df contains all coverage info per-exon. 
                        df = pd.merge(
                            summed_edits, 
                            summed_total_coverage, 
                            how='outer', 
                            left_on=['chrom','start','end','name','score','strand'],
                            right_on=['chrom','start','end','name','score','strand']
                        )
                        df['edited_over_edited_c'] = df['edit_coverage']/df['editable_coverage']
                        # (8) Calculate the coverage across all C's in each window
                        df['all_c_coverage'] = df.apply(get_total_c_coverage, args=(rdd, genome_fa, ), axis=1)
                        df['edited_over_all_c'] = df['edit_coverage']/df['all_c_coverage']

                        df.fillna(0, inplace=True)
                        # reorder columns to match
                        df = df[[
                            'chrom','start','end','name','score',
                            'strand','edit_coverage','editable_coverage',
                            'edited_over_edited_c','all_c_coverage','edited_over_all_c'
                        ]]
                        # print("DF: ", df)
                        all_scores.append(df)
                
        except KeyError as e: # this takes care of ambiguous genes (ie. edits that were annotated with more than one gene)
            pass
        progress.update(1)
        pybedtools.cleanup()
        
    progress = trange(len(all_scores), desc='Concatenating all scores...')
    for s in all_scores:
        all_scores_df = pd.concat([all_scores_df, s])
        progress.update(1)
            
    # (9) Since we actually want per-gene coverage info instead of per-exon coverage, I'll sum these up here
    region_edited = pd.DataFrame(all_scores_df.groupby(['name'])['edit_coverage'].sum())
    region_edited_c = pd.DataFrame(all_scores_df.groupby(['name'])['editable_coverage'].sum())
    region_all_c = pd.DataFrame(all_scores_df.groupby(['name'])['all_c_coverage'].sum())
    
    all_regions = pd.merge(region_edited, region_edited_c, how='outer', left_index=True, right_index=True)
    all_regions = pd.merge(all_regions, region_all_c, how='outer', left_index=True, right_index=True)
    all_regions['edited_over_edited_c'] = all_regions['edit_coverage']/all_regions['editable_coverage']
    all_regions['edited_over_all_c'] = all_regions['edit_coverage']/all_regions['all_c_coverage']
    # (10) Save the summed coverage info per-gene
    all_regions.sort_values(by=['edited_over_all_c'], ascending=False).to_csv(
        output_file_summed, 
        sep='\t', 
        index=True, 
        header=True
    )
    # (11) Save the coverage info per-exon
    all_scores_df.sort_values(by=['chrom','start','end','strand']).to_csv(
        output_file, 
        sep='\t', 
        index=False, 
        header=True
    )
    
def main():
    parser = ArgumentParser(
        description="Creates windows surrounding edit sites and calculates the read coverage across edits and editable sites."
    )
    parser.add_argument(
        "--annotated_edits_file", 
        help="input file from joined edits (conf score should be 4th column, \
        it should be annotated using my annotator script.)", 
        required=True, 
        default=None
    )
    parser.add_argument(
        "--bg_edits_file", 
        help="background edits we want to subtract.", 
        required=False, 
        default=None
    )
    parser.add_argument("--conf", help="conf score", required=False, type=float, default=0.)
    parser.add_argument("--gtfdb", help="GFFUtils DB file", required=True, default=None)
    parser.add_argument("--output_file", help="output bedfile", required=False, default=None)
    parser.add_argument("--output_file_summed", help="output bedfile merged by gene", required=False, default=None)
    parser.add_argument("--cds_file", help="CDS coordinates (merged) for each gene in BED6 format.", required=False, default=None)
    parser.add_argument("--three_prime_utr_file", help="three prime UTR coordinates (merged) for each gene in BED6 format.", required=False, default=None)
    parser.add_argument("--five_prime_utr_file", help="five prime UTR coordinates (merged) for each gene in BED6 format.", required=False, default=None)
    parser.add_argument("--genome_fa", help="genome fasta file", required=True, default=None)
    parser.add_argument("--chrom_sizes_file", help="genome chrom.sizes file", required=True, default=None)
    parser.add_argument("--pos_bw", help="bigwig file (positive)", required=True, default=None)
    parser.add_argument("--neg_bw", help="bigwig file (negative)", required=True, default=None)
    args = parser.parse_args()
    annotated_edits_file = args.annotated_edits_file
    bg_edits_file = args.bg_edits_file
    conf = args.conf
    gtfdb = args.gtfdb
    pos_bw = args.pos_bw
    neg_bw = args.neg_bw
    five_prime_utr_file = args.five_prime_utr_file
    three_prime_utr_file = args.three_prime_utr_file
    cds_file = args.cds_file
    
    rdd = ReadDensity(
        pos=pos_bw,
        neg=neg_bw,  
    )
    
    output_file = "{}.{}.exons.bed".format(
        annotated_edits_file,
        conf,
    ) if args.output_file is None else args.output_file
    output_file_summed = "{}.{}.exons_merged.bed".format(
        annotated_edits_file,
        conf,
    ) if args.output_file_summed is None else args.output_file_summed
    
    genome_fa = args.genome_fa
    chrom_sizes_file = args.chrom_sizes_file

    score_edits_exons(
        annotated_edits_file=annotated_edits_file, 
        bg_edits_file=bg_edits_file,
        output_file=output_file, 
        output_file_summed=output_file_summed,
        conf=conf, 
        genome_fa=genome_fa,
        chrom_sizes_file=chrom_sizes_file,
        rdd=rdd,
        cds_file=cds_file,
        five_prime_utr_file=five_prime_utr_file,
        three_prime_utr_file=three_prime_utr_file
    )

if __name__ == '__main__':
    main()