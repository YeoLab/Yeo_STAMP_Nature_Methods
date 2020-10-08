#!/usr/bin/env python
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import glob
import pybedtools
from Bio import SeqIO
from tqdm import trange
from collections import defaultdict, OrderedDict

def filter_edits(scored_edits, conf, score_column, score, output_file):
    df = pd.read_csv(scored_edits, sep='\t')
    print("Number of edits before filtering: {}".format(df.shape[0]))
    df = df[df['score'] >= conf]
    df = df[df[score_column] >= score]
    print("Number of edits after filtering: {}".format(df.shape[0]))
    df['name'] = str(score_column)
    df = df[['chrom', 'start', 'end', 'name', score_column, 'strand']]
    bt = pybedtools.BedTool.from_dataframe(df)
    df = bt.merge(s=True, c=(4,5,6), o=('count','mean','distinct')).to_dataframe()
    print("Number of edits after merging: {}".format(df.shape[0]))
    df['length'] = df['end'] - df['start']
    fig, ax = plt.subplots()
    ax.hist(df['length'], bins=100)
    fig.savefig(output_file + ".peak_lengths.svg")
    df.to_csv(output_file, sep='\t', index=False, header=False)

def main():
    parser = ArgumentParser(
        description="This takes an annotated and scored set of edits " \
        "(output from score_edits_total_coverage.py) and filters based on score column." \
        "Expects each scored edits to have the following columns: " \
        "[chrom, start, end, name, score, strand, edit_coverage, editable_coverage, edited_over_edited_c, " \
        "all_c_coverage, edited_over_all_c]. Output is a BED6 file."
    )
    parser.add_argument("--scored_edits", help="output from score_edits_total_coverage.py", required=True, default=None)
    parser.add_argument("--conf", help="conf filter if we wanted to pre-filter edits by the fourth \"conf\" column (Default: 0)", required=False, type=float, default=0.)
    parser.add_argument("--score_column", help="Default: edited_over_all_c (edit/c)", required=False, type=str, default="edited_over_all_c")
    parser.add_argument("--score", help="Select windows whose score is at least defined by score_column", required=True, type=float, default=0.01)
    parser.add_argument("--output_file", help="output filtered file", required=False, default=None)
    
    args = parser.parse_args()
    
    scored_edits = args.scored_edits
    conf = args.conf
    score_column = args.score_column
    score = args.score
    output_file = "{}.{}.bed".format(
        scored_edits,
        score,
    ) if args.output_file is None else args.output_file
    
    filter_edits(
        scored_edits=scored_edits, 
        conf=conf, 
        score_column=score_column,
        score=score,
        output_file=output_file
    )
    
if __name__ == '__main__':
    main()
