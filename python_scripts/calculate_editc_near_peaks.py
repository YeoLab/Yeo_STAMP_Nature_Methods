#!/usr/bin/env python
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import pandas as pd
import numpy as np
import os
import pybedtools
import warnings
import subprocess
import sys
import pandas.api.types as ptypes
from collections import defaultdict
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.set_option('display.max_colwidth', 500)
pd.set_option('display.max_columns', 500)
from tqdm import trange

import score_edits_total_coverage as se


## helper func
def get_ratio_real_rand(index, df):
    """
    Calculates the ratio of % of real peaks that meet edit/c 
    cutoffs (blue bars) to % of random peaks that meet the 
    same cutoff (orange bars) in our main figure. 
    """
    fraction = df.loc[index]['fraction']/df.loc[index]['mean_rand_fraction']
    return fraction


def counts_to_rpkm(feature_counts_table):
    """
    RPKM transformation func
    """
    counts = feature_counts_table.iloc[:,5:]
    lengths = feature_counts_table['Length']
    mapped_reads = counts.sum()
    return (counts * pow(10,9)).div(mapped_reads, axis=1).div(lengths, axis=0)


def counts_to_tpm(counts_table, skip_col=5):
    """
    simple function that converts a featureCounts pandas Dataframe
    into a TPM dataframe.
    
    :param counts_table: pandas.DataFrame() 
        either a featureCounts table (first five cols contain non-count info,
        the rest contain raw counts) or a generic counts table (use skip_col=0
        in this case)
    :return tpm: pandas.DataFrame
    """
    rpkm = counts_to_rpkm(counts_table)
    tpm = rpkm.div(rpkm.sum())*pow(10,6)
    return tpm


def read_and_create_SAF(bed_file, output_dir):
    """
    Create a SAF file of just 3'UTRs and CDS. The 3'UTR/CDS 
    length is going to be what RPKM/TPM is calculated upon.
    Expects a BED file of regions with geneids in the 'name'
    column.
    """
    try:
        regions_file_for_random_subset_SAF = os.path.join(
            output_dir, 
            os.path.splitext(os.path.basename(bed_file))[0] + ".SAF"
        )

        utr_bed = pd.read_csv(
            bed_file, 
            sep='\t', 
            names=['chrom','start','end','name','score','strand']
        )
        utr_bed['start1base'] = utr_bed['start'] + 1
        utr_saf = utr_bed[['name','chrom','start1base','end','strand']]
        utr_saf.to_csv(
            regions_file_for_random_subset_SAF, 
            sep='\t', 
            index=False, 
            header=False
        )
        return regions_file_for_random_subset_SAF
    except Exception:
        print("Had trouble creating SAF file.")
        sys.exit(1)

        
def read_and_filter_annotated_edits_file(edit_file, bg_edit_file, conf=0):
    """
    This function expects an annotated (brianv1) edits file that has been 
    modified such that the 'score' column contains a comma-delimited list 
    of edited,total coverage for each edit site. "edited" is the number of 
    "high" quality (defined by bcftools call) C>T conversions, while 
    "total coverage" is total quality coverage. We can modify this function 
    later, but for now we're just looking at edits that were found to be 
    within the 3'UTR or CDS region.
    
    This function will also filter for low-confidence edit sites if applicable, 
    and filter out sites that were also found in a background (ApoControl) set. 
    
    Returns a bedtool of filtered edit sites whose thickStart contains edited
    coverage and thickEnd contains total. 
    """
    fg_sites_in_region = pd.read_csv(
        edit_file, 
        sep='\t', 
        names=['chrom','start','end','conf','score','strand','geneid','genename','region','annotation']
    )
    
    print("Number of total edits: {}".format(fg_sites_in_region.shape[0]))
    fg_sites_in_region = fg_sites_in_region[
        (fg_sites_in_region['region']=='CDS') | (fg_sites_in_region['region']=='3utr')
    ]
    fg_sites_in_region = fg_sites_in_region[
        fg_sites_in_region['conf']>=conf
    ]
    print("Number of edits after region filtering: {}".format(fg_sites_in_region.shape[0]))

    # thickStart = edited #
    fg_sites_in_region.loc[:,'thickStart'] = fg_sites_in_region['score'].apply(
        lambda x: int(x.split(',')[0])
    )
    # thickEnd = total coverage # 
    fg_sites_in_region.loc[:,'thickEnd'] = fg_sites_in_region['score'].apply(
        lambda x: int(x.split(',')[1])
    )

    fg_sites_in_region = fg_sites_in_region[
        ['chrom','start','end','geneid','conf','strand','thickStart','thickEnd']
    ]

    edits = pybedtools.BedTool.from_dataframe(fg_sites_in_region).sort()
    unfiltered_edits = pybedtools.BedTool.from_dataframe(fg_sites_in_region).sort()
    
    print("Before removing bg (apo) edits: {}".format(unfiltered_edits.to_dataframe().shape[0]))
    if bg_edit_file is not None:
        edits = unfiltered_edits.intersect(bg_edit_file, s=True, v=True)
    else:
        edits = unfiltered_edits
    print("After removing bg (apo) edits: {}".format(edits.to_dataframe().shape[0]))

    return edits


def run_featurecounts(regions_file_for_random_subset_SAF, bam_file, output_dir):
    """
    Using a SAF file generated from some annotation, runs featureCounts. 
    We don't actually use this other than to gauge expression, so I'm going to be
    pretty permissive and count everything. 
    """
    try:
        output_counts = os.path.join(output_dir, '{}.counts.txt'.format(os.path.basename(bam_file).split('.')[0]))
        output_stdout = os.path.join(output_dir, '{}.counts.out'.format(os.path.basename(bam_file).split('.')[0]))
        # Define and run featureCounts to generate counts table of BAM file
        cmd = "module load subreadfeaturecounts;featureCounts "
        cmd += "-a {} ".format(regions_file_for_random_subset_SAF) 
        cmd += "-F SAF -M -O -s 2 -o {} ".format(output_counts)
        cmd += "{}".format(bam_file)
        
        with open(output_stdout, 'w') as o:
            subprocess.check_call(cmd, shell=True, stdout=o)
        # Transform counts into TPM and find the lowest TPM where we can detect an edit.
        return output_counts
    except Exception as e:
        print("Had trouble running featureCounts: {}".format(e))
        sys.exit(1)

        
def filter_tpm(counts_file, rand_bed, bam_file, edited_genes, output_dir):
    """
    bam_file: used as a column name only. This function takes a counts file 
    from featureCounts and performs TPM normalization. Uses the set of genes 
    that were edited and finds the lower bound TPM required to call an edit. 
    Then, filter genes that do not meet this requirement from our pool of 
    random 3'UTR+CDS (rand_bed) and write to a file. 
    
    params: 
    
    counts_file: string
    rand_bed: string
    bam_file: string
    edited_genes: set
    output_dir: string
    """
    try:
        regions_file_for_random_subset_expressed = os.path.join(
            output_dir, 
            os.path.splitext(os.path.basename(rand_bed))[0] + ".expressed.bed"
        )
        
        counts_table = pd.read_csv(counts_file, index_col=0, skiprows=1, sep='\t')
        tpm = counts_to_tpm(counts_table)        
        edited_genes_tpm = tpm.loc[edited_genes]
        edited_genes_tpm.dropna(inplace=True) # drops ambiguously annotated edits, usually comma-delimited (ie. ENSG0000001,ENSG0000002)
        
        # edge case - most (all?) genes with editing should not have a TPM of zero
        print(edited_genes_tpm[edited_genes_tpm[bam_file]==0].shape[0])
        edited_genes_tpm = edited_genes_tpm[edited_genes_tpm[bam_file]>0]
        tpm_cutoff = edited_genes_tpm[bam_file].min()
        print(tpm_cutoff)
        # Select only the genes in our rand_bed file that have at least this TPM
        bed = pd.read_csv(
            rand_bed, 
            sep='\t', 
            names=['chrom','start','end','name','score','strand']
        )
        # expressed = list(expressed_genes.index)
        expressed = list(tpm[tpm[bam_file]>=tpm_cutoff].index)
        print("Number of expressed genes with TPM cutoff of {}: {}".format(tpm_cutoff, len(expressed)))
        bed_expressed = bed[bed['name'].isin(expressed)]
        bed_expressed.to_csv(regions_file_for_random_subset_expressed, sep='\t', index=False, header=False)
        return regions_file_for_random_subset_expressed
    
    except Exception as e:
        print("Had trouble filtering by TPM: {}".format(e))
        sys.exit(1)
    
    
def read_and_filter_rand_regions(rand_bed, edited_genes, bam_file, peak_bedtool, motif_bedtool, d, chrom_sizes, output_dir, full_prefix):
    """
    MAIN function for generating our background set of regions which are 
    used as a pool for randomly placing peaks in. This requires the following steps:
    1. Start with a set of regions (all 3'UTR+CDS regions as defined by Gencode)
    2. Calculate any expressed 3'UTR+CDS using the BAM file from standard RNASeq pipeline.
    3. Intersect with annotated edits and find the lowest expressed gene that has an edit.
        A. Use this as a lower-bound TPM and remove any non-expressed gene from our pool.
    4. Remove any peak + slopped region so we don't accidentally randomly place a region inside a real peak.
    5. Remove any motif + slopped region so we don't accidentally randomly place a region next to a known motif.
    """
    try:
        regions_file_for_random_subset_SAF = read_and_create_SAF(
            bed_file=rand_bed, 
            output_dir=output_dir
        )
        output_counts = run_featurecounts(
            regions_file_for_random_subset_SAF=regions_file_for_random_subset_SAF, 
            bam_file=bam_file,
            output_dir=output_dir
        )
        regions_file_for_random_subset_expressed = filter_tpm(
            counts_file=output_counts, 
            rand_bed=rand_bed, 
            bam_file=bam_file, 
            edited_genes=edited_genes, 
            output_dir=output_dir
        )
        filtered_regions_for_random_selection_motif_and_peak_filtered = regions_file_for_random_subset_expressed + ".filtered.bed"
        
        # Must not be near a real peak.
        avoid_these_real_peak_regions = peak_bedtool.slop(b=d, g=chrom_sizes).sort()
        avoid_these_real_peak_regions = avoid_these_real_peak_regions.saveas(
            os.path.join(output_dir, '{}.peak_slopped.bed'.format(full_prefix))
        )
        # Must not contain any motif
        avoid_these_motifs_regions = motif_bedtool.slop(b=d, g=chrom_sizes).sort()
        avoid_these_motifs_regions.saveas(
            os.path.join(output_dir, '{}.tmp_motif_slopped.bed'.format(full_prefix))
        )

        # Combine the two peak/motif bedfiles. Write this to a tmp file that we'll use to remove from our random pool.
        avoid_these_regions = pybedtools.BedTool.from_dataframe(pd.concat(
            [
                avoid_these_real_peak_regions.to_dataframe(), 
                 avoid_these_motifs_regions.to_dataframe()
            ]
        )).sort()
        avoid_these_regions.saveas(os.path.join(output_dir, '{}.regions_to_avoid.bed'.format(full_prefix)))
        
        # this is the bedtool containing 3'UTRs and CDS for expressed genes:
        unfiltered_regions_for_random_selection = pybedtools.BedTool(
            regions_file_for_random_subset_expressed
        ).sort()
        
        filtered_regions_for_random_selection = unfiltered_regions_for_random_selection.subtract(avoid_these_regions, s=True)
        print(
            "After filtering out real regions, the regions_file_for_random_subset_expressed looks like this: {}".format(
                filtered_regions_for_random_selection.to_dataframe().shape[0]
            )
        )
        
        filtered_regions_for_random_selection.saveas(filtered_regions_for_random_selection_motif_and_peak_filtered)
        return filtered_regions_for_random_selection_motif_and_peak_filtered
    except Exception as e:
        print("Had trouble filtering random bed file: {}".format(e))
        sys.exit(1)

        
def read_and_return_peak_bedtool(fn, l10p, l2fc, annotated='motif'):
    """
    Because we're dealing with peaks from various sources/annotations, 
    I'm using this function to universally read in and return peaks.
    
    annotated: string
        Use this to define how peaks were annotated. 
        Can be: "ericv1", "ericv2", "motif", or "no"
    """
    if annotated == 'ericv1': # this is an 'eric-style' annotation
        peaks = pd.read_csv(
            fn, names=[
                'chrom','start','end','l10p','l2fc','strand','annotation','geneid'
            ], 
            sep='\t'
        )
        peaks['region'] = peaks['annotation'].apply(lambda x: x.split('|')[0])
        peaks = peaks[(peaks['region']=='CDS') | (peaks['region']=='3utr')]
        peaks = peaks[(peaks['l10p']>=l10p) & (peaks['l2fc']>=l2fc)]
        peaks = peaks[['chrom','start','end','geneid','l2fc','strand']]
    elif annotated == 'ericv2':
        peaks = pd.read_csv(
            fn, 
            names=[
                'chrom','start','end','l10p','l2fc','strand',
                'annotation','annotation2','geneid','genename','region_w_overlap'
            ], 
            sep='\t'
        )
        peaks['region'] = peaks['region_w_overlap'].apply(lambda x: x.split('|')[0])
        peaks = peaks[(peaks['region']=='CDS') | (peaks['region']=='3utr')]
        peaks = peaks[(peaks['l10p']>=l10p) & (peaks['l2fc']>=l2fc)]
        peaks = peaks[['chrom','start','end','geneid','l2fc','strand']]
    elif annotated == 'motif': # I pre-parse the eCLIP peaks in another notebook and in doing so, I lose the l10p and l2fc information. That's okay, I filter in that notebook prior.
        peaks = pd.read_csv(fn, names=['chrom','start','end','l10p','l2fc','strand'], sep='\t')
    elif annotated == 'no': # basic input-normalized BED6
        peaks = pd.read_csv(fn, names=['chrom','start','end','l10p','l2fc','strand'], sep='\t')
        peaks = peaks[(peaks['l10p']>=l10p) & (peaks['l2fc']>=l2fc)]
    # merge neighboring peaks just to make sure we're not overlapping still
    print("Number of peaks after filtering for significance: {}".format(peaks.shape[0]))
    unmerged_peaks = pybedtools.BedTool.from_dataframe(peaks).sort()
    merged_peaks = unmerged_peaks.merge(d=1, c=(4,5,6), o=('collapse','collapse','distinct')).sort()
    print("Number of peaks after merging neighboring peaks: {}".format(merged_peaks.to_dataframe().shape[0]))
    return merged_peaks


def slop_and_merge(peaks, d, g):
    """
    Just a function that wraps a bedtools slop() with merge() 
    to make sure that after extending peak regions that we don't 
    also introduce overlaps. 
    """
    # print("Number of peaks before merge: {}".format(len(peaks)))
    unmerged = peaks.slop(b=d, g=g)
    merged = unmerged.merge(s=True, c=(4,5,6), o=('collapse','collapse','distinct')).sort()
    # print("Number of peaks after merge: {}".format(len(merged)))
    return merged


def plot_histogram_peak_lengths(bedtool, d, output_dir):
    """
    Uses a bedtool of peaks to generate a histogram of lengths. 
    May be useful in determining the size of our edit windows for HOMER calculation.
    Returns the median length of all peak.
    """
    windows = bedtool.to_dataframe()
    windows['length'] = windows['end'] - windows['start']
    std = windows['length'].std()
    
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.axvline(windows['length'].median(), ls='--', label='median: {:.2f}'.format(windows['length'].median()))
    ax.axvline(windows['length'].mean(), label='mean: {:.2f}'.format(windows['length'].mean()))
    ax.hist(windows['length'], alpha=0.3)
    ax.legend()
    fig.suptitle("Peak lengths after slop (d={}) and merging (std={:.2f})".format(d, std))
    fig.savefig(os.path.join(output_dir, 'peak_lengths.svg'))
    plt.close(fig)
    return windows['length'].median()
    
    
def compute_window_and_edit_fractions(peaks, edits, d, pos_bw, neg_bw, chrom_sizes, genome_fa, output_dir):
    
    rdd = se.ReadDensity(
        pos=pos_bw,
        neg=neg_bw,  
    )

    peak_windows = slop_and_merge(peaks, d=d, g=chrom_sizes)
    plot_histogram_peak_lengths(bedtool=peak_windows, d=d, output_dir=output_dir)
    ### 
    # Intersect slopped windows with edits, giving us all edits that overlap the region.
    # Foreach edit position, we preserve edited coverage (number of reads that were 
    # C>T converted at that position) in the thickStart column. We also preserve the 
    # total coverage at each edited position in the thickEnd column. For readability, 
    # I've renamed these columns "edit_coverage" and "editable_coverage". 
    # 
    # @see read_and_filter_annotated_edits_file() to see how we create the edits
    # BedTool from an annotated edits file. 
    ###
    intersected_edits = peak_windows.intersect(
        edits, s=True, wa=True, loj=True
    ).to_dataframe(
        names=[
            'chrom','start','end','name','score','strand',
            'edit_chrom','edit_start','edit_end','edit_name',
            'edit_score','edit_strand','edit_coverage','editable_coverage'
        ]
    )
    ###
    # Peaks that do not contain any intersected edits will return . and . for 
    # edit_coverage and editable_coverage. This also causes the dtypes for these 
    # columns to be non-integer Objects. Convert and coerce to int-types here.
    ###
    if not ptypes.is_numeric_dtype(intersected_edits['edit_coverage']):
        intersected_edits['edit_coverage'].replace({'.':0}, inplace=True)
    if not ptypes.is_numeric_dtype(intersected_edits['editable_coverage']):
        intersected_edits['editable_coverage'].replace({'.':-1}, inplace=True)
    intersected_edits[["edit_coverage", "editable_coverage"]] = intersected_edits[
        ["edit_coverage", "editable_coverage"]
    ].apply(pd.to_numeric)
    assert ptypes.is_numeric_dtype(intersected_edits['edit_coverage'])
    assert ptypes.is_numeric_dtype(intersected_edits['editable_coverage'])
    ###
    # Compress all edit sites within a peak here, summing up the total edit
    # coverage and editable coverage across the entire region. 
    ###

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
    # merge on peak regions, now we have (foreach peak) total edited reads and total coverage across these edited reads.
    df = pd.merge(
        summed_edits, 
        summed_total_coverage, 
        how='outer', 
        left_on=['chrom','start','end','name','score','strand'],
        right_on=['chrom','start','end','name','score','strand']
    )
    df['edited_over_edited_c'] = df['edit_coverage']/df['editable_coverage']

    # (8) Calculate the coverage across ALL C's in each window
    df['all_c_coverage'] = df.apply(se.get_total_c_coverage, args=(rdd, genome_fa, ), axis=1)
    df['edited_over_all_c'] = df['edit_coverage']/df['all_c_coverage']

    # reorder columns to match
    df = df[[
        'chrom','start','end','name','score',
        'strand','edit_coverage','editable_coverage',
        'edited_over_edited_c','all_c_coverage','edited_over_all_c'
    ]]
    ###
    # Don't worry about peaks without coverage, which will end up as NaN types anyway. 
    # Expected behavior divides edited coverage (0) with total coverage (0): 0/0 = NaN
    # Let's also assert regions of 0 total coverage do not have any edit coverage.
    ###
    df.dropna(inplace=True)
    try:
        assert df[(df['edit_coverage']>0) & (df['all_c_coverage']==0)].shape[0] == 0
    except AssertionError:
        print("We found sites that have edit coverage but no total coverage?")
        print(df[(df['edit_coverage']>0) & (df['all_c_coverage']==0)])
        sys.exit(1)
    # Just a QC/sanity check to make sure that we're getting roughly enough random regions with editing. 
    number_windows_with_edits = len(peak_windows) - intersected_edits[intersected_edits['edit_coverage']==0].shape[0]
    
    return df, number_windows_with_edits

def choose_bins(df, n_bins):
    """
    Simply gets the max edit/C fraction and returns equally sized n_bins.
    These will be our cutoffs
    """
    # n = int(df.shape[0]*0.005)  # if we need to remove outliers
    # print("Setting the ranges without the top {} values".format(n))
    min_score = 0
    max_score = df.sort_values('edited_over_all_c', ascending=False)['edited_over_all_c'].max()
    print("max score: ", max_score)
    cutoffs = np.arange(min_score, max_score, max_score/n_bins)
    ticklabels = [round(c, 4) for c in cutoffs]
    return list(cutoffs), ticklabels

def randomize_peaks(bedtool, incl, g='hg19'):
    """
    Uses our peaks and randomly places equally-sized regions along a set of 
    3'UTR+CDS (incl) regions. As of bedtools 2.27, bedtools.shuffle() ignores 
    strand info, so there is a chance these regions will not be placed on 
    a correct 3'UTR/CDS (and removed downstream due to lack of coverage). 
    To try and correct for this, this function doubles the number of regions 
    so all things being equal, should return about the same number of regions 
    with proper coverage for a more accurate background comparison. 
    """
    random1 = bedtool.shuffle(
        genome=g, 
        incl=incl, 
        noOverlapping=True,
        chromFirst=True
    ).sort()
    random2 = bedtool.shuffle(
        genome=g, 
        incl=incl, 
        noOverlapping=True,
        chromFirst=True
    ).sort()

    random = pybedtools.BedTool.from_dataframe(
        pd.concat([random1.to_dataframe(), random2.to_dataframe()])
    ).sort()
    return random


def calculate_zscore(real_bins, rand_bins, output_dir, full_prefix):
    try:
        zscore_table = pd.merge(
            real_bins[['fraction']], 
            rand_bins[['mean_rand_fraction','standard_deviation']], 
            how='inner', 
            left_index=True, 
            right_index=True
        )
        zscore_table['zscore'] = (zscore_table['fraction'] - zscore_table['mean_rand_fraction'])/zscore_table['standard_deviation']
        zscore_table.fillna(0, inplace=True)
        zscore_table.to_csv(os.path.join(output_dir, "{}.zscores.txt".format(
            full_prefix
        )), sep='\t')
        return zscore_table
    except Exception as e:
        print("Had trouble calculating zscores: {}".format(e))
        sys.exit(1)
    
def main():
    parser = ArgumentParser(
        description="This is notebook 15 v4 in script form."
    )
    parser.add_argument(
        "--rand_bed", 
        help="BED file containing regions that may be selected as pseudobinding regions. \
        These regions will be further filtered of motifs (defined by --motif_file) and \
        of true-positives (defined by --peak) prior to random selection.", 
        required=False, 
        default='/projects/ps-yeolab3/bay001/annotations/hg19/gencode_v19/hg19_v19_cds_and_three_prime_utrs.bed'
    )
    parser.add_argument(
        "--motif_file", 
        help="BED file containing motifs to remove from random pool.", 
        required=True, 
    )

    parser.add_argument(
        "--edits_file", 
        help=\
        "Edits file in BED6 format (edit#,cov# must be in the \
        'name' column, comma-delimited)", 
        required=True
    )
    parser.add_argument(
        "--bg_edits_file", 
        help="Background file to subtract edits from. Typically \
        an APOBEC control", 
        required=False,
        default=None
    )
    parser.add_argument(
        "--peak", 
        help="Peak (with or without motif) file in annotated or BED format", 
        required=True
    )
    parser.add_argument(
        "--annotated_type", 
        help="can be: inputnorm, ericv1, ericv2, brianv1, motif", 
        required=True
    )
    parser.add_argument(
        "--d", 
        help="Distance to extend each peak region to (Default: 25nt)", 
        required=False, 
        type=int, 
        default=25
    )
    parser.add_argument(
        "--bam_file", 
        help="BAM file used to calculate the minimal TPM \
        threshold to call random peaks from.", 
        required=False, 
        default=None
    )
    parser.add_argument(
        "--pos_bw", 
        help="positive bigwig file for coverage calculations", 
        required=False, 
        default=None
    )
    parser.add_argument(
        "--neg_bw", 
        help="negative bigwig file for coverage calculations", 
        required=False, 
        default=None
    )
    parser.add_argument(
        "--genome_fa", 
        help="genome fasta file", 
        required=False, 
        default='/projects/ps-yeolab3/bay001/annotations/hg19/hg19.fa'
    )
    parser.add_argument(
        "--chrom_sizes", 
        help="genome fasta file", 
        required=False, 
        default='/projects/ps-yeolab3/bay001/annotations/hg19/hg19.chrom.sizes'
    )
    parser.add_argument(
        "--num_rand_trials", 
        help="Number of random trials to perform. Default: 100", 
        required=False, 
        default=100,
        type=int
    )
    parser.add_argument(
        "--output_dir", 
        help="output directory.", 
        required=True, 
        type=str
    )
    parser.add_argument(
        "--conf", 
        help="dont count edits with low conf (disabled)", 
        required=False, 
        default=0,
        type=float
    )
    parser.add_argument(
        "--l10p", 
        help="l10p threshold for peak", 
        required=False, 
        default=3,
        type=int
    )
    parser.add_argument(
        "--l2fc", 
        help="log2 fold threshold for peak", 
        required=False, 
        default=3,
        type=int
    )
    parser.add_argument(
        "--bins", 
        help="Number of bins to divide bargraphs by", 
        required=False, 
        default=10,
        type=int
    )
    ###
    # Parse arguments
    ###
    args = parser.parse_args()
    rand_bed = args.rand_bed
    motif_file = args.motif_file
    annotated_type = args.annotated_type
    edits_file = args.edits_file
    bg_edits_file = args.bg_edits_file
    peak = args.peak
    d = args.d
    bam_file = args.bam_file
    pos_bw = args.pos_bw
    neg_bw = args.neg_bw
    genome_fa = args.genome_fa
    chrom_sizes = args.chrom_sizes
    num_rand_trials = args.num_rand_trials
    output_dir = args.output_dir
    conf = args.conf
    l10p = args.l10p
    l2fc = args.l2fc
    n_bins = args.bins
    
    edits_file_prefix = ".".join(os.path.splitext(os.path.basename(edits_file))[0].split('.')[:3])
    peak_prefix = os.path.splitext(os.path.basename(peak))[0]
    
    full_prefix = "{}.{}.d{}.c{}.l2fc{}.l10p{}".format(
        edits_file_prefix, peak_prefix, d, conf, l2fc, l10p
    )
    ### 
    # Read in our edits and peaks file
    ###
    edits = read_and_filter_annotated_edits_file(
        edit_file=edits_file, 
        bg_edit_file=bg_edits_file, 
        conf=conf
    )
    
    peaks = read_and_return_peak_bedtool(
        fn=peak,
        l10p=l10p,
        l2fc=l2fc,
        annotated=annotated_type
    )
    
    ###
    # Generate a background dataset
    ###
    motifs = pybedtools.BedTool(motif_file).sort()
    filtered_rand_regions_file = read_and_filter_rand_regions(
        rand_bed=rand_bed,
        edited_genes=set(edits.to_dataframe()['name']),
        bam_file=bam_file,
        peak_bedtool=peaks,
        motif_bedtool=motifs,
        d=d,
        chrom_sizes=chrom_sizes,
        output_dir=output_dir,
        full_prefix=full_prefix
    )
    
    ###
    # Calculate edit/c over real peak regions
    ###
    df_real, target_num = compute_window_and_edit_fractions(
        peaks=peaks, 
        edits=edits, 
        d=d, 
        pos_bw=pos_bw, 
        neg_bw=neg_bw, 
        chrom_sizes=chrom_sizes, 
        genome_fa=genome_fa,
        output_dir=output_dir
    )
    cutoffs, ticks = choose_bins(df_real, n_bins=n_bins)
    
    ###
    # Calculate edit/c over random pseudopeak regions
    ###
    df_rand_window_shapes = []
    df_rand_tot = []
    progress = trange(num_rand_trials)
    for i in range(num_rand_trials):
        random = randomize_peaks(
            bedtool=peaks,
            incl=filtered_rand_regions_file
        )
        df_rand, _ = compute_window_and_edit_fractions(
            peaks=random, 
            edits=edits, 
            d=d, 
            pos_bw=pos_bw, 
            neg_bw=neg_bw, 
            chrom_sizes=chrom_sizes, 
            genome_fa=genome_fa,
            output_dir=output_dir
        )
        df_rand_window_shapes.append(df_rand.shape[0])
        df_rand_tot.append(df_rand)
        progress.update(1)
    
    ###
    # Generate QC plot to make sure we're getting enough random 
    # peaks with enough coverage/edits. @see randomize_peaks() 
    # for a better explanation.
    ###
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.axvline(df_real.shape[0], linestyle='-', alpha=0.7)
    ax = plt.hist(df_rand_window_shapes, alpha=0.3)
    fig.suptitle("Number of pseudo peaks w/ countable reads vs actual")
    fig.savefig(os.path.join(output_dir, '{}.pseudo_peaks_w_countable_reads.svg'.format(
        full_prefix
    )))
    plt.close(fig)
    ###
    # Tidy results into seaborn-compatible dataframe
    ###
    rand_bins = defaultdict(list)

    for t in range(num_rand_trials):
        for cutoff, tick in zip(cutoffs, ticks):
            trial = df_rand_tot[t]
            rand_bins[tick].append(
                (trial[trial['edited_over_all_c']>=cutoff].shape[0])/(trial.shape[0])
            )

    rand_bins = pd.DataFrame(rand_bins).T
    rand_bins.to_csv(os.path.join(output_dir, '{}.random_editc_fractions_per_bin.tsv'.format(
        full_prefix
    )), sep='\t')
    rand_bins_melted = rand_bins.reset_index().melt(
        id_vars='index', value_name='fraction', var_name='trial'
    ).set_index('index')
    rand_bins_melted['class'] = 'rand'
    
    real_bins = defaultdict(dict)
    for cutoff, tick in zip(cutoffs, ticks):
        real_bins[tick]['fraction'] = (df_real[df_real['edited_over_all_c']>=cutoff].shape[0])/(df_real.shape[0])

    real_bins = pd.DataFrame(real_bins).T
    real_bins['class'] = 'real'
    real_bins['trial'] = 0
    real_bins.to_csv(os.path.join(output_dir, '{}.real_editc_fractions_per_bin.tsv'.format(
        full_prefix
    )), sep='\t')
    
    bins = pd.concat([real_bins, rand_bins_melted])
    bins.fillna(0, inplace=True)
    
    ###
    # Seaborn will plot standard deviation for us, but we need 
    # to manually get the numbers for z-score calculation
    ###
    rand_bins['mean_rand_fraction'] = rand_bins.mean(axis=1)
    rand_bins['standard_deviation'] = rand_bins[list(range(0, num_rand_trials))].std(axis=1)
    
    zscore_table = calculate_zscore(
        real_bins=real_bins, 
        rand_bins=rand_bins, 
        output_dir=output_dir, 
        full_prefix=full_prefix
    )
    
    ### 
    # Compute blue vs orange bar @see get_ratio_real_rand()
    ###
    ratio = defaultdict()
    for i in zscore_table.index:
        ratio[i] = get_ratio_real_rand(i, zscore_table)

    ratio = pd.DataFrame(ratio, index=['ratio']).T.reset_index()
    
    bins.reset_index(inplace=True)
    
    ###
    # Plot the main figure. 
    ###
    fig, ax = plt.subplots(figsize=(15, 10))

    ax = sns.barplot(x='index', y='fraction', hue='class', data=bins, alpha=0.3, ci='sd', errwidth=1)
    num_k = int(len(ax.patches)/2)
    for i in range(n_bins):
        # get_x pulls left or right; get_height pushes up or down
        ax.text(
            ax.patches[i].get_x()-.03, 
            ax.patches[i].get_height()+0.01, 
            "Z-score: {:.2f}".format(zscore_table.iloc[i]['zscore']), 
            fontsize=8,
            color='dimgrey'
        )
    ax2 = ax.twinx()
    ratio['ratio'].plot(ax=ax2)
    ax2.set_ylim(0, 20)
    ax.set_xlim(-1, n_bins)
    ax.set_xlabel("edit/c threshold")
    ax.set_ylabel("fraction of peaks containing edit/c scores above threshold")
    ax.set_title("Fraction of real/rand peaks +/- {}bp (n={})\nwith an edit/C score >= threshold.".format(
        d, peaks.to_dataframe().shape[0]
    ))
    fig.savefig(os.path.join(output_dir, "{}.svg".format(
        full_prefix
    )))
    plt.close(fig)
if __name__ == '__main__':
    main()