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
from tqdm import tqdm
import multiprocessing
import shutil
import time
from collections import defaultdict, OrderedDict
pd.options.mode.chained_assignment = None  # default='warn'
import sys

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
    NOTE: Only works with one continuous window. Will break if there is more than one interval. 
    """
    try:
        assert bedtool.to_dataframe().shape[0] == 1
    except AssertionError:
        print("This function does not work with multiple regions: {}".format(bedtool.to_dataframe()))
        sys.exit(1)
        
    d = {}
    pos = bedtool.to_dataframe().iloc[0]
    chrom = pos.chrom
    start = pos.start
    end = pos.end
    strand = pos.strand
    # print("{}:{}-{}:{}".format(chrom, start, end, strand))
    sequences = bedtool.sequence(fi=genome_fa, s=False, name=True) 
    # print(sequences)
    with open(sequences.seqfn) as f:
        for record in SeqIO.parse(f, "fasta"):
            # print(record.seq.upper())
            if strand == '+':
                relpos = find(record.seq.upper(), 'C')
            elif strand == '-':
                relpos = find(record.seq.upper(), 'G')
            else:
                print("Strand error: {}".format(strand))
                sys.exit(1)
    abspos = ["{}:{}".format(chrom, start + p) for p in relpos]
    coverage = rdd.values(chrom=chrom, start=start, end=end, strand=strand)
    coverage = [np.abs(c) for c in coverage]  # doesn't matter for how we're making bigwigs, but just to be sure. 
    c_coverage = [coverage[p] for p in relpos]
    for p, c in zip(abspos, c_coverage):
        d[p] = c
    return d


def get_total_c_coverage(row, rdd, genome_fa):
    """
    Just wraps sum_all_c_coverage_in_window() to work with pandas row
    """
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


def read_and_filter_editing_sites(sites_file, conf, src):
    
    if src == 'confeditswap':
        # If the edit coverage is in the 'score' column formatted as: 1,5 where 1 is edited coverage, 5 is the total 
        header = [
            'chrom','start','end','name','score','strand'
        ]
        sites = pd.read_csv(sites_file, names=header, sep='\t')
        sites = sites[sites['name'] >= conf]
        sites['thickStart'] = sites['score'].apply(lambda x: int(x.split(',')[0]))
        sites['thickEnd'] = sites['score'].apply(lambda x: int(x.split(',')[1]))
        sites = sites[['chrom','start','end','name','score','strand','thickStart','thickEnd']]
        return pybedtools.BedTool.from_dataframe(sites).sort()
    elif src == 'editconfswap':
        # If the edit coverage is in the 'score' column formatted as: 1,5 where 1 is edited coverage, 5 is the total 
        header = [
            'chrom','start','end','score','name','strand'
        ]
        sites = pd.read_csv(sites_file, names=header, sep='\t')
        sites = sites[sites['name'] >= conf]
        sites['thickStart'] = sites['score'].apply(lambda x: int(x.split(',')[0]))
        sites['thickEnd'] = sites['score'].apply(lambda x: int(x.split(',')[1]))
        sites = sites[['chrom','start','end','name','score','strand','thickStart','thickEnd']]
        return pybedtools.BedTool.from_dataframe(sites).sort()
    elif src == 'sailor':
        # If the edit coverage is embedded in the 'name' as: 13|G>A|0.076923077 where 13 is the total coverage, 0.07.. is edit %
        header = [
            'chrom','start','end','name','score','strand',
        ]
        
    elif src == 'brianv1':
        # If the edit coverage is embedded in the 'name' as: 13|G>A|0.076923077 where 13 is the total coverage, 0.07.. is edit %
        header = [
            'chrom','start','end','name','score','strand',
            'gene_id','gene_name','region','annotation_string'
        ]
    else:
        print("Unsupported annotation, for now.")
        sys.exit(1)
        
    sites = pd.read_csv(sites_file, names=header, sep='\t')
    sites['thickStart'] = sites['name'].apply(
        lambda x: int(int(x.split('|')[0]) * float(x.split('|')[2]))
    )
    sites['thickEnd'] = sites['name'].apply(
        lambda x: int(x.split('|')[0])
    )
    sites = sites[sites['score'] >= conf]
    sites = sites[['chrom','start','end','name','score','strand','thickStart','thickEnd']]
    return pybedtools.BedTool.from_dataframe(sites).sort()

def editc(bed_file, sites_file, conf, rdd, genome_fa, src):
    tqdm.pandas()
    
    sites = read_and_filter_editing_sites(sites_file=sites_file, conf=conf, src=src)
    bedtool = pybedtools.BedTool(bed_file).sort()
    
    intersected_edits = bedtool.intersect(
        sites, s=True, wa=True, loj=True
    ).to_dataframe(
        names=[
            'chrom','start','end','name','score','strand',
            'edit_chrom','edit_start','edit_end','edit_name',
            'edit_score','edit_strand','edit_coverage','editable_coverage'
        ]
    )
    
    intersected_edits['edit_coverage'].fillna(-1, inplace=True)
    intersected_edits['editable_coverage'].fillna(-1, inplace=True)
    # We can to calculate across all exons, regardless of whether or not we found an edit there. 
    try:
        intersected_edits['edit_coverage'].replace({'.':0}, inplace=True)
        intersected_edits['editable_coverage'].replace({'.':0}, inplace=True)
    except TypeError:  # it will complain that dtype of edit_coverage and editable_coverage isn't a string
        pass
    # pybedtools does some weird things with mixed dtypes, so let's make sure we convert to INT here (so 1+1 = 2, not "11")
    intersected_edits['edit_coverage'] = intersected_edits['edit_coverage'].astype(int)
    intersected_edits['editable_coverage'] = intersected_edits['editable_coverage'].astype(int)
    
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
        )['editable_coverage'].sum()
    ).reset_index()
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
    df['all_c_coverage'] = df.progress_apply(get_total_c_coverage, args=(rdd, genome_fa, ), axis=1)
    df['edited_over_all_c'] = df['edit_coverage']/df['all_c_coverage']
    df.fillna(0, inplace=True)
    # reorder columns to match
    df = df[[
        'chrom','start','end','name','score',
        'strand','edit_coverage','editable_coverage',
        'edited_over_edited_c','all_c_coverage','edited_over_all_c'
    ]]
    return df

def editc_parallel(bed_df, sites_file, conf, pos_bw, neg_bw, genome_fa, src, chrom, i):
    tqdm.pandas(desc=chrom, position=i)
    rdd = ReadDensity(
        pos=pos_bw,
        neg=neg_bw,  
    )
    sites = read_and_filter_editing_sites(sites_file=sites_file, conf=conf, src=src)
    bedtool = pybedtools.BedTool.from_dataframe(bed_df).sort()
    
    intersected_edits = bedtool.intersect(
        sites, s=True, wa=True, loj=True
    ).to_dataframe(
        names=[
            'chrom','start','end','name','score','strand',
            'edit_chrom','edit_start','edit_end','edit_name',
            'edit_score','edit_strand','edit_coverage','editable_coverage'
        ]
    )
    
    intersected_edits['edit_coverage'].fillna(-1, inplace=True)
    intersected_edits['editable_coverage'].fillna(-1, inplace=True)
    # We can to calculate across all exons, regardless of whether or not we found an edit there. 
    try:
        intersected_edits['edit_coverage'].replace({'.':0}, inplace=True)
        intersected_edits['editable_coverage'].replace({'.':0}, inplace=True)
    except TypeError:  # it will complain that dtype of edit_coverage and editable_coverage isn't a string
        pass
    # pybedtools does some weird things with mixed dtypes, so let's make sure we convert to INT here (so 1+1 = 2, not "11")
    intersected_edits['edit_coverage'] = intersected_edits['edit_coverage'].astype(int)
    intersected_edits['editable_coverage'] = intersected_edits['editable_coverage'].astype(int)
    
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
        )['editable_coverage'].sum()
    ).reset_index()
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
    df['all_c_coverage'] = df.progress_apply(get_total_c_coverage, args=(rdd, genome_fa, ), axis=1)
    df['edited_over_all_c'] = df['edit_coverage']/df['all_c_coverage']
    df.fillna(0, inplace=True)
    # reorder columns to match
    df = df[[
        'chrom','start','end','name','score',
        'strand','edit_coverage','editable_coverage',
        'edited_over_edited_c','all_c_coverage','edited_over_all_c'
    ]]
    pybedtools.cleanup()
    return df

def main():
    parser = ArgumentParser(
        description="Creates windows surrounding edit sites and calculates the read coverage across edits and editable sites."
    )
    parser.add_argument(
        "--bed_file", 
        help="", 
        required=True, 
        default=None
    )
    parser.add_argument(
        "--edits_file", 
        help="Combine the fwd and rev outputs from SAILOR.", 
        required=True, 
        default=None
    )
    parser.add_argument("--conf", help="pre-process, remove any edit prior to edit/C calculation.", required=False, default=0, type=float)
    parser.add_argument("--output_file", help="output bedfile", required=False, default=None)
    parser.add_argument("--genome_fa", help="genome fasta file", required=True, default=None)
    parser.add_argument("--pos_bw", help="bigwig file (positive)", required=True, default=None)
    parser.add_argument("--neg_bw", help="bigwig file (negative)", required=True, default=None)
    parser.add_argument("--editsrc", help="src", required=False, default='confeditswap')
    parser.add_argument("--threads", help="number of threads (Default: 8)", required=False, default=8, type=int)
    args = parser.parse_args()
    bed_file = args.bed_file
    sites_file = args.edits_file
    pos_bw = args.pos_bw
    neg_bw = args.neg_bw
    conf = args.conf
    genome_fa = args.genome_fa
    editsrc = args.editsrc
    threads = args.threads
    
    rdd = ReadDensity(
        pos=pos_bw,
        neg=neg_bw,  
    )
    
    output_file = "{}_editc.txt".format(
        os.path.splitext(bed_file)[0],
    ) if args.output_file is None else args.output_file
    
    genome_fa = args.genome_fa

    bed_df = pd.read_csv(bed_file, sep='\t', names=['chrom','start','end','name','score','strand'])
    bed_df.sort_values(by=['chrom','start','end','strand'], inplace=True)
    print("Total line num: {}".format(bed_df.shape[0]))
    if bed_df.shape[0] <= threads:
        editc_df = editc(
            bed_file=bed_file,
            sites_file=sites_file,
            conf=conf,
            rdd=rdd,
            genome_fa=genome_fa,
            src=editsrc
        )
    else:
        df_split = []
        for chromosome in set(bed_df['chrom']):
            df_split.append(bed_df[bed_df['chrom']==chromosome])
        # bed_df, sites_file, conf, rdd, genome_fa, src
        args_iter = []
        print("Splitting into {} parts: ".format(len(df_split)))
        i = 0 # position for tqdm
        for ds in df_split:
            chrom = list(set(ds['chrom']))[0]
            print("Sizeof part {}: {}".format(chrom, ds.shape[0]))
            args_iter.append((ds, sites_file, conf, pos_bw, neg_bw, genome_fa, editsrc, chrom, i))
            i += 1
        with multiprocessing.Pool(processes=threads) as pool:
            editc_df = pd.concat(pool.starmap(editc_parallel, args_iter))
            pool.close()
            pool.join()
            
        editc_df.drop_duplicates(inplace=True)

    editc_df.to_csv(output_file, sep='\t', index=False, header=False)
    
if __name__ == '__main__':
    main()