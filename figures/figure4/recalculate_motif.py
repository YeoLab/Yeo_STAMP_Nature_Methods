#!/usr/bin/env python
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import pandas as pd
import os
import subprocess
import glob
import pysam
import sys
from collections import defaultdict

# If findMotifs.pl crashes (due to not finding any known motifs within targets and raising divideByZero errors, let's just read some stand-in file where the p-value is 1)
known_results_null = '/home/bay001/projects/kris_apobec_20200121/scripts/knownResultsNull.txt'

known_results_columns = [
    'Motif Name', 'Consensus', 'P-value', 'Log P-value', 'q-value (Benjamini)', 
    '# of Target Sequences with Motif', '% of Target Sequences with Motif', 
    '# of Background Sequences with Motif', '% of Background Sequences with Motif'
]


def run_homer_findmotifs(real_fa, rand_fa, mknown, output_dir, p):
    cmd = 'findMotifs.pl {} fasta {} -nofacts -p {} -rna -S 20 -len 6 -noconvert -nogo -known -fasta {} -mknown {}'.format(
        real_fa, output_dir, p, rand_fa, mknown
    )
    with open(os.devnull, 'w') as fnull:
        subprocess.check_call(cmd, shell=True)
    try:
        return pd.read_table(os.path.join(output_dir, 'knownResults.txt'), sep='\t')
    except IOError as e:
        print("{}, output_dir: {}".format(e, output_dir))
        return pd.read_table(known_results_null, sep='\t')
    
def get_fasta_from_analyze_motifs(analyze_motifs_output, region):
    """
    Uses the standard output structure from analyze_motifs to return 
    the paths to real and random fasta files used to compute original de-novo motifs.
    """
    basename = os.path.basename(os.path.dirname(analyze_motifs_output))
    real_fa = os.path.join(analyze_motifs_output, 'fasta', '{}.bed.{}.{}.fa'.format(basename, region, 'real'))
    rand_fa = os.path.join(analyze_motifs_output, 'fasta', '{}.bed.{}.{}.fa'.format(basename, region, 'random'))
    try:
        assert os.path.exists(real_fa)
    except AssertionError:
        print(real_fa)
    assert os.path.exists(rand_fa)
    return real_fa, rand_fa
    
    
def get_pvalue_from_knownresults(real_fa, rand_fa, mknown, output_dir, p, label):
    df = run_homer_findmotifs(real_fa, rand_fa, mknown, output_dir, p)
    df.index = [label]
    df.columns = known_results_columns
    return df


def main():
    parser = ArgumentParser(description="Runs HOMER to re-calculate motifs. This script takes outputs from analyze_motifs (from clip analysis) and an existing motif to re-calculate p-values for that motif.")
    parser.add_argument("--input_dirs", help="Outputs from analyze_motifs. Should contain fasta files representing fg and randomized bg peaks", required=True, default=None, nargs='+')
    parser.add_argument("--mknown", help="known motif for which to re-calculate p-values for.", required=True, type=str)
    parser.add_argument("--output_dir", help="output file where the new results will go.", required=True, type=str, default=None)
    parser.add_argument("--p", "--processors", help="Processors for Homer", required=False, type=int, default=4)
    parser.add_argument("--region", help="Region within analyze_motifs/output", required=False, type=str, default='all')
    
    args = parser.parse_args()
    input_dirs = args.input_dirs
    mknown = args.mknown
    output_dir = args.output_dir
    p = args.p
    region = args.region
    merged = pd.DataFrame(columns=known_results_columns)
    for input_dir in input_dirs:
        label = os.path.basename(os.path.dirname(input_dir))
        real_fa, rand_fa = get_fasta_from_analyze_motifs(analyze_motifs_output=input_dir, region=region)
        df = get_pvalue_from_knownresults(
            real_fa=real_fa,
            rand_fa=rand_fa,
            mknown=mknown,
            output_dir=os.path.join(output_dir, os.path.basename(os.path.dirname(input_dir))),
            p=p,
            label=label
        )
        merged = pd.concat([merged, df])
        
    merged.to_csv(os.path.join(output_dir, 'p-values.tsv'), sep='\t')
    
    
if __name__ == '__main__':
    main()
