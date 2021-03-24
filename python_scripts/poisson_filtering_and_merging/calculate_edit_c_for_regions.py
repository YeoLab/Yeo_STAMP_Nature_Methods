#!/usr/bin/env python
from pyfaidx import Fasta

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
from collections import defaultdict, OrderedDict
pd.options.mode.chained_assignment = None  # default='warn'
import sys
import pandas as pd
import re
import json
from multiprocessing import Pool
import time

class ReadDensity:
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

def get_c_positions_and_coverage_in_window(chrom, start, end, strand, rdd):
    """
    Given region specified by chrom start end and strand, return the positions (0-based) of C's (or G's on neg strand).
    """ 
    d = {}
    # print("{}:{}-{}:{}".format(chrom, start, end, strand))
    
    # print(sequences)
    FA = Fasta(fasta, rebuild=False)
    
    try:
        if 'chr' not in chrom:
            sequence = FA['chr{}'.format(chrom)][start:end].seq
        else:
            sequence = FA[chrom][start:end].seq
    except:
        if chrom == 'chrMT':
            sequence = FA['chrM'][start:end].seq
    
    # print(record.seq.upper())
    if strand == '+':
        matches = re.finditer('C', sequence.upper())
    elif strand == '-':
        matches = re.finditer('G', sequence.upper())
    else:
        print("Strand error: {}".format(strand))
        sys.exit(1)
        
    relpos = [match.start() for match in matches]
    abspos = ["{}:{}".format(chrom, start + p) for p in relpos]
    coverage = rdd.values(chrom=chrom, start=start, end=end, strand=strand)
    coverage = [np.abs(c) for c in coverage]  # doesn't matter for how we're making bigwigs, but just to be sure. 
    c_coverage = [coverage[p] for p in relpos]
    for p, c in zip(abspos, c_coverage):
        d[p] = c
    return d

def sum_all_c_coverage_in_window(chrom, start, stop, strand, rdd):
    all_coverage = 0
    c_positions = get_c_positions_and_coverage_in_window(chrom, start, stop, strand, rdd)
    for pos, cov in c_positions.items():
        all_coverage += cov
    return int(all_coverage)
    

def total_conversions_in_window(stamp_sites, chrom, start, end, strand):
    return stamp_sites[(stamp_sites.start_stamp > start) &
                 (stamp_sites.end_stamp < end) & 
                 (stamp_sites.chrom_stamp == chrom) &
                 (stamp_sites.strand_stamp == strand)
                ].num_edited_stamp.sum()


def get_edit_c_info(t):
    region_label, chrom, start, end, strand = t
        
    edit_c_info = {}
        
    d = ReadDensity(
        pos=forward_bw,
        neg=reverse_bw
    )

    total_c_coverage = sum_all_c_coverage_in_window(chrom, start, end, strand, d)
    
    total_conversions = total_conversions_in_window(stamp_sites, chrom, start, end, strand)
        
    edit_c_info[region_label] = {
                            'start': start,
                            'end': end,
                            'chrom': chrom,
                            'strand': strand,
                             'coverage': total_c_coverage,
                             'conversions': total_conversions,
                             'edit_c': total_conversions/total_c_coverage
                            }
    region_level_edit_c = pd.DataFrame(edit_c_info).T
    
    return region_level_edit_c


def load_stamp_sites(annotated_stamp_sites_filepath):
    stamp_sites = pd.read_csv(annotated_stamp_sites_filepath,
                                             sep='\t', names=['chrom_stamp', 'start_stamp', 'end_stamp', 'conf_stamp', 'coverage_stamp', 'strand_stamp', 
                                                              'geneid', 'genename', 'region', 'annot'])
    stamp_sites['num_edited_stamp'] = [int(i.split(',')[0]) for i in stamp_sites.coverage_stamp]
    stamp_sites['total_coverage_stamp'] = [int(i.split(',')[1]) for i in stamp_sites.coverage_stamp]
    return stamp_sites


def load_json_info(input_json_filepath):
    # Load parameters/filepaths for stamp sites, bigwigs, and regions to probe editC in
    with open(input_json) as f:
        data = json.load(f)

    output_folder = data.get('output_folder')
    label = data.get('label')
    annotated_stamp_sites_file = data.get('annotated_stamp_sites_file')
    forward_bw = data.get('forward_bw')
    reverse_bw = data.get('reverse_bw')
    fasta = data.get('fasta')
    regions_for_edit_c = data.get('regions_for_edit_c')

    print("Label: {}".format(label))
    print("Annotated STAMP sites: {}".format(annotated_stamp_sites_file))
    print("Forward bigwig: {}".format(forward_bw))
    print("Reverse bigwig: {}".format(reverse_bw))
    print("Reference being used: {}".format(fasta))
    print("Regions being annotated: {}\m".format(regions_for_edit_c))
    
    # Load STAMP Sites and regions to probe
    print('\n\nLoading annotated STAMP sites...')
    stamp_sites = load_stamp_sites(annotated_stamp_sites_file)
    # Regions should have columns: 
    # region_id, chrom, start, end, strand
    
    print('\n\nLoading regions to annotate with editC information -- should have following columns: region_id, chrom, start, end, strand')
    regions = pd.read_csv(regions_for_edit_c, sep='\t')
    
    return output_folder,label, stamp_sites, forward_bw, reverse_bw, fasta, regions 

def edit_c(output_folder,label, stamp_sites, forward_bw, reverse_bw, fasta, regions_for_edit_c):
    region_ids = list(regions.region_id)
    chroms = list(regions.chrom)
    starts = list(regions.start)
    ends = list(regions.end)
    strands = list(regions.strand)

     # Set up a multi-processing pool to process each region individually and add editC information 
    print("\n\nCalculating editC information for erach region using multi-processing pool... this step could take a while (fastest when ppn = 8)")
    num_processes = 8

    p = Pool(num_processes)
    region_level_edit_c_dfs = p.map(get_edit_c_info, zip(region_ids, chroms, starts, ends, strands))
    p.close()
    p.join()

    out_file = '{}/{}_edit_c_for_all_regions.tsv'.format(output_folder, label)
    print("\n\nOutputting results to {}".format(out_file))
    edit_c_for_all_regions = pd.concat(region_level_edit_c_dfs)
    edit_c_for_all_regions = edit_c_for_all_regions.fillna(int(0))
    edit_c_for_all_regions.to_csv(out_file, sep='\t', index=True, header=True)
    

input_json = sys.argv[1]
print('input is {}'.format(input_json))

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Calculating Edit-C Region Baseline")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
output_folder,label, stamp_sites, forward_bw, reverse_bw, fasta, regions  = load_json_info(input_json)
print('loaded, calculating editC')
edit_c(output_folder,label, stamp_sites, forward_bw, reverse_bw, fasta, regions)