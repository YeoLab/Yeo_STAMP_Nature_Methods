#!/usr/bin/env python
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import numpy as np
import pandas as pd
import os
import glob
import pysam
import sys
from Bio import SeqIO
from tqdm import trange
from collections import defaultdict

def get_readcount(bam_file):
    """
    Parses a bam file idxstats to get the number of reads.
    The BAM file MUST have an index.
    """
    num_reads = pysam.idxstats(
        bam_file
    ).split('\n')
    nums = {}
    for num in num_reads:
        try:
            chrom, chrlen, mapped, unmapped = num.split('\t')
            nums[chrom] = int(mapped) + int(unmapped)
        except ValueError:
            print(num)
    return pd.DataFrame(nums, index=['num']).T.sum().values[0]

def split_bams_group(bam_file, output_file, barcodes, v=False):
    """
    Basically reads through the whole 10X BAM file and looks for barcodes 
    specified in the 'barcodes' parameter. If it finds a read belonging to that 
    barcode, it will write to a new BAM file containing 
    
    params:
    
    bam_file: string
        path to the possorted bam file (with MD tags! not relevant for this script but absolutely necessary for edit calls)
    output_file: string
        output bam file including (or excluding with -v) the barcodes.
    barcodes: set
        set of barcodes {"TGCGGGTGTGTTCCTC-1", "TGCGGGTGTGTTCCTG-1"...}
    """
    
    
    # For progress bar.
    try:
        progress = trange(get_readcount(bam_file))  # the total number of reads in the bam file

        s = pysam.AlignmentFile(bam_file, "rb")
        split_bam = pysam.AlignmentFile(
            output_file, 
            "wb", 
            template=s
        )
        for read in s:
            if not read.is_unmapped and not read.is_secondary:  # get only primary mapped reads.
                try:
                    barcode = read.get_tag('CB')  # identify read=assigned barcode
                    gene_annotation = read.get_tag('GX')  # 
                    if v:
                        if barcode not in barcodes:
                            split_bam.write(read)
                    else:
                        if barcode in barcodes:
                            split_bam.write(read)
                except KeyError:
                    pass # no barcode or no GX tag, so 10X doesn't count this read and we shouldn't either.
            progress.update(1)
        s.close()
        split_bam.close()
        return 0
    except Exception as e:
        print(e)
        sys.exit(1)

def split_bams(bam_file, output_dir, barcodes):
    """
    Basically reads through the whole 10X bam file and looks for barcodes 
    specified in the 'barcodes' parameter. If it finds a read belonging to that 
    barcode, it will append that read to the dictionary 'reads_dict'.
    
    The result will be a dictionary of something like:
        {
            'TGCGGGTGTGTTCCTC-1': [read1, read2, read3, ... readN],
            'ACGGTTAAGTGGTTAA-1': [read1, read2, ...]
        }
    """

    reads_dict = defaultdict(list)
    # For progress bar.
    progress = trange(get_readcount(bam_file))  # the total number of reads in the bam file
    
    s = pysam.AlignmentFile(bam_file, "rb")
    
    for read in s:
        if (not read.is_unmapped) and (not read.is_secondary) and (not read.is_duplicate):  # get only primary mapped reads.
            try:
                barcode = read.get_tag('CB')  # identify read=assigned barcode
                gene_annotation = read.get_tag('GX')  # 
                if barcode in barcodes:
                    reads_dict[barcode].append(read)
            except KeyError:
                pass # no barcode or no GX tag, so 10X doesn't count this read and we shouldn't either.
        progress.update(1)
    s.close()
    return reads_dict

def main():
    parser = ArgumentParser(description="Compares two bed files")
    parser.add_argument("--possorted_bam_file", help="10X bamfile (WITH MD TAG! use samtools calmd on native 10X output).", required=True, default=None)
    parser.add_argument("--barcodes_file", help="file with a list of barcodes to split", required=True, type=str)
    parser.add_argument("--output_dir", help="output file where the split bam files will all go.", required=False, type=str, default=None)
    parser.add_argument("--group", default=False, action="store_true", help="if specified, reads pertaining to each barcode will be merged into one bam file instead of individually splitting them.")
    parser.add_argument("--v", default=False, action="store_true", help="if specified, the resulting bam file will be one SUBTRACTED of the barcode list. Only applies when --group is specified.")
    
    args = parser.parse_args()
    possorted_bam_file = args.possorted_bam_file
    barcodes_file = args.barcodes_file
    output_dir = args.output_dir if args.output_dir is not None else os.path.dirname(possorted_bam_file)
    group = args.group
    v = args.v
    
    barcodes = pd.read_csv(barcodes_file, names=['barcodes'])
    barcodes = set(barcodes['barcodes']) # searching set is faster.
    
    
    if not group:
        # If we want to split each barcode into its own file.
        print("Splitting the BAM file...")
        reads_dict = split_bams(
            bam_file=possorted_bam_file, 
            output_dir=output_dir, 
            barcodes=barcodes
        )
        print("Saving outputs...")

        samfile = pysam.AlignmentFile(possorted_bam_file, "rb")
        for barcode in reads_dict.keys():
            split_bam = pysam.AlignmentFile(
                os.path.join(
                    output_dir, "{}-{}.bam".format(
                        os.path.splitext(os.path.basename(possorted_bam_file))[0], barcode
                    )
                ), 
                "wb", 
                template=samfile
            )
            for read in reads_dict[barcode]:
                split_bam.write(read)
        split_bam.close()
        samfile.close()
    else:
        # If we want subset the bam (all barcodes in list go into one grouped barcoded bam).
        output_file = os.path.join(
            output_dir, "{}-{}.bam".format(
                os.path.splitext(os.path.basename(possorted_bam_file))[0], os.path.basename(barcodes_file)
            )
        )
        split_bams_group(bam_file=possorted_bam_file, output_file=output_file, barcodes=barcodes, v=v)
        
if __name__ == '__main__':
    main()
