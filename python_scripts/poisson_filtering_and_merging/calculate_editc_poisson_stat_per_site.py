from multiprocessing import Pool
from scipy.stats import poisson
import pandas as pd
import numpy as np
import json
import sys

def apply_edit_c_pvalue_to_site(r, baseline_edit_c):    
    total_coverage = r.total_coverage_stamp
    total_edits = r.num_edited_stamp
    
    baseline_num_edits = baseline_edit_c * total_coverage
    
    return 1 - poisson.cdf(total_edits, baseline_num_edits, loc=0)

def poisson_filter(region_info_tuple):
    region_id, chrom, start, end, strand = region_info_tuple
        
    region_sites = stamp_sites[
        (stamp_sites.chrom_stamp == chrom) &
        (stamp_sites.start_stamp >= start) &
        (stamp_sites.end_stamp <= end) &
        (stamp_sites.strand_stamp == strand)
    ].copy()

    baseline_edit_c = edit_c_output.loc[region_id].edit_c

    if baseline_edit_c > 0:
        # If there were any edits at all in this region, then calculate the p-value for this site
        region_sites['poisson_p'] = region_sites.apply(apply_edit_c_pvalue_to_site, args=(baseline_edit_c,), axis=1)

        adjusted_p = 0.05/len(region_sites)
        region_sites['adjusted_p_cutoff'] = adjusted_p

        region_sites['passed_poisson'] = region_sites.poisson_p < adjusted_p
        region_sites['region_edit_c'] = baseline_edit_c
        region_sites['region_id'] = region_id

        return region_sites[['chrom_stamp', 'start_stamp', 'end_stamp', 'strand_stamp', 'region_id', 'region_edit_c', 'poisson_p', 'adjusted_p_cutoff', 'passed_poisson']]
    else:
        # If there were no edits in this region, then just return an empty dataframe
        return pd.DataFrame(columns=['chrom_stamp', 'start_stamp', 'end_stamp', 'strand_stamp', 'region_id', 'region_edit_c', 'poisson_p', 'adjusted_p_cutoff', 'passed_poisson'])

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
    print('\nLoading annotated STAMP sites...')
    stamp_sites = load_stamp_sites(annotated_stamp_sites_file)
    # Regions should have columns: 
    # region_id, chrom, start, end, strand
    
    print('\nLoading regions to annotate with editC information -- should have following columns: region_id, chrom, start, end, strand')
    regions = pd.read_csv(regions_for_edit_c, sep='\t')
    
    return output_folder,label, stamp_sites, forward_bw, reverse_bw, fasta, regions 

def to_int(i):
    return int(i)

def add_index(r):
    return '{}:{}-{}'.format(r.chrom_stamp, r.start, r.end)

def calculate_poisson_stats(output_folder,label, stamp_sites, forward_bw, reverse_bw, fasta, regions):
    # Spread the information from our regions 
    region_ids = list(regions.region_id)
    chroms = list(regions.chrom)
    starts = list(regions.start)
    ends = list(regions.end)
    strands = list(regions.strand)

    num_processes = 8
    p = Pool(num_processes)
    all_poisson_filtered = p.map(poisson_filter, zip(region_ids, chroms, starts, ends, strands))
    p.close()
    p.join()

    features_with_poisson_info = pd.concat(all_poisson_filtered)

    features_with_poisson_info['start'] = features_with_poisson_info.start_stamp.apply(to_int)
    features_with_poisson_info['end'] = features_with_poisson_info.end_stamp.apply(to_int)

    features_with_poisson_info.index = features_with_poisson_info.apply(add_index, axis=1)
    
    return features_with_poisson_info

input_json = sys.argv[1]
print('input is {}'.format(input_json))

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("~~Calculating Poisson Statistics~~")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
output_folder,label, stamp_sites, forward_bw, reverse_bw, fasta, regions  = load_json_info(input_json)
edit_c_filename = '{}/{}_edit_c_for_all_regions.tsv'.format(output_folder, label)
print('Loading file containing per region baseline edit c information: {}'.format(edit_c_filename))
edit_c_output = pd.read_csv(edit_c_filename, sep='\t', index_col=0)


features_with_poisson_info = calculate_poisson_stats(output_folder,label, stamp_sites, forward_bw, reverse_bw, fasta, regions)
# For sites that overlapped more than one region, be conservative and keep the site if it passes poission test for at least one of those regions
features_with_poisson_info = features_with_poisson_info.sort_values('passed_poisson', ascending=False).drop_duplicates(['chrom_stamp', 'start_stamp', 'end_stamp'])
features_with_poisson_info.to_csv('{}/{}_poisson_filtering_information.tsv'.format(output_folder, label), sep='\t', index=True)



