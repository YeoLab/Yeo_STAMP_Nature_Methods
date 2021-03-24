from multiprocessing import Pool
from scipy.stats import poisson
import pandas as pd
import numpy as np
import json
import pybedtools
import sys

def add_index(r):
    return '{}:{}-{}'.format(r.chrom_stamp, r.start_stamp, r.end_stamp)

def load_stamp_sites(annotated_stamp_sites_filepath):
    stamp_sites = pd.read_csv(annotated_stamp_sites_filepath,
                                             sep='\t', names=['chrom_stamp', 'start_stamp', 'end_stamp', 'conf_stamp', 'coverage_stamp', 'strand_stamp', 
                                                              'geneid', 'genename', 'region', 'annot'])
    stamp_sites['num_edited_stamp'] = [int(i.split(',')[0]) for i in stamp_sites.coverage_stamp]
    stamp_sites['total_coverage_stamp'] = [int(i.split(',')[1]) for i in stamp_sites.coverage_stamp]
    stamp_sites.index = stamp_sites.apply(add_index, axis=1)
    stamp_sites['SiteID'] = stamp_sites.index
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
    print('Loading annotated STAMP sites...')
    stamp_sites = load_stamp_sites(annotated_stamp_sites_file)
    # Regions should have columns: 
    # region_id, chrom, start, end, strand
    
    print('Loading regions to annotate with editC information -- should have following columns: region_id, chrom, start, end, strand')
    regions = pd.read_csv(regions_for_edit_c, sep='\t')
    
    return output_folder,label, stamp_sites, forward_bw, reverse_bw, fasta, regions 


def get_merged_stamp_windows(sites):
    max_distance = 100
    stamp_bedtool = pybedtools.BedTool.from_dataframe(sites[['chrom_stamp', 'start_stamp', 'end_stamp', 'SiteID', 'geneid', 'strand_stamp']])
    stamp_windows_merged = stamp_bedtool.sort().merge(s=True, d=max_distance,
                          c=(4,5, 6), 
                          o=('distinct', 'distinct', 'distinct')).to_dataframe(names=['chrom', 'start', 'end', 'ids', 'ids_2', 'strand'])
    return stamp_windows_merged


def write_print(output_folder, label, string, action='a'):
    p_shrinkers_label = p_shrinkers
    if not filter_with_poisson:
        p_shrinkers_label = 'no_poisson'
    log_file = '{}/{}_filter_merge_log_file_ps_{}_st_{}.txt'.format(output_folder, label, p_shrinkers_label, score_thresholds)
    with open(log_file, action) as f:
        print(string)
        f.write('{}\n'.format(string))
    
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("~~~~~Merging and Filtering...~~~~~~")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
input_json = sys.argv[1]
print('input is {}'.format(input_json))

output_folder,label, stamp_sites, forward_bw, reverse_bw, fasta, regions  = load_json_info(input_json)

poisson_outputs = pd.read_csv('{}/{}_poisson_filtering_information.tsv'.format(output_folder, label), sep='\t', index_col=0)

score_thresholds = [0.9, 0.99, 0.999]

#score_thresholds = [0.7, 0.8, 0.9]

#score_thresholds = [0, 0.75, 0.95]

filter_with_poisson = False
#p_shrinkers = [1, 10, 100]

p_shrinkers = [100]

density_filters = [True, False]

write_print(output_folder, label, '_____________{}_____________'.format(label), action='w')
for p_shrinker in p_shrinkers:
    for score_threshold in score_thresholds:
        for density_filter in density_filters:
            write_print(output_folder, label, '________________________________________')
            write_print(output_folder, label, 'p-boost: {} score threshold : {}, density filter: {}'.format(p_shrinker, score_threshold, density_filter))
            # Before any filters are applied
            num_original_stamp_sites = len(stamp_sites)
            merged_windows = get_merged_stamp_windows(stamp_sites)
            num_original_stamp_merged_windows = len(merged_windows)


            write_print(output_folder, label, 'Before filtering: {} {}'.format(num_original_stamp_sites, num_original_stamp_merged_windows))




            # Filter with Poisson info
            if filter_with_poisson:
                stamp_sites_with_poisson_info = stamp_sites.join(poisson_outputs, how='left', rsuffix='_poisson')
                filtered_stamp_sites_with_poisson_info = stamp_sites_with_poisson_info[
                        stamp_sites_with_poisson_info.poisson_p < stamp_sites_with_poisson_info.adjusted_p_cutoff/p_shrinker]
                num_stamp_sites_after_poisson_filtering = len(filtered_stamp_sites_with_poisson_info)
                merged_windows_after_poisson_filtering = get_merged_stamp_windows(filtered_stamp_sites_with_poisson_info)
                num_merged_windows_after_poisson_filtering = len(merged_windows_after_poisson_filtering)

                write_print(output_folder, label, 'After poisson filter: {} {}'.format(num_stamp_sites_after_poisson_filtering, num_merged_windows_after_poisson_filtering))
            else:
                write_print(output_folder, label, "No poisson filter applied...")
                filtered_stamp_sites_with_poisson_info = stamp_sites
    

            # Filter with score
            stamp_sites_filtered_by_poisson_and_score = filtered_stamp_sites_with_poisson_info[filtered_stamp_sites_with_poisson_info.conf_stamp > score_threshold]
            num_stamp_sites_filtered_by_poisson_and_score = len(stamp_sites_filtered_by_poisson_and_score)
            merged_windows_after_poisson_and_score_filtering = get_merged_stamp_windows(stamp_sites_filtered_by_poisson_and_score)
            num_merged_windows_after_poisson_and_score_filtering = len(merged_windows_after_poisson_and_score_filtering)



            write_print(output_folder, label, 'After poisson and score filters: {} {}'.format(num_stamp_sites_filtered_by_poisson_and_score, num_merged_windows_after_poisson_and_score_filtering))



            # isolated site filter
            if density_filter:
                merged_windows_after_poisson_and_score_filtering['length'] = merged_windows_after_poisson_and_score_filtering['end'].subtract(merged_windows_after_poisson_and_score_filtering['start'])

                num_isolated_sites = len(merged_windows_after_poisson_and_score_filtering[merged_windows_after_poisson_and_score_filtering.length == 1])
                num_sites_after_poisson_and_score_filtering_and_isolation_filter = num_stamp_sites_filtered_by_poisson_and_score - num_isolated_sites
                merged_windows_after_poisson_and_score_filtering_and_isolation_filter = merged_windows_after_poisson_and_score_filtering[merged_windows_after_poisson_and_score_filtering.length > 1]
                num_merged_windows_after_poisson_and_score_filtering_and_isolation_filter = len(merged_windows_after_poisson_and_score_filtering_and_isolation_filter)
                write_print(output_folder, label, 'After poisson, score, and isolation filters: {} {}'.format(num_sites_after_poisson_and_score_filtering_and_isolation_filter, num_merged_windows_after_poisson_and_score_filtering_and_isolation_filter))
            else:
                merged_windows_after_poisson_and_score_filtering_and_isolation_filter = merged_windows_after_poisson_and_score_filtering


            write_print(output_folder, label, 'Total merged windows: {}'.format(len(merged_windows_after_poisson_and_score_filtering_and_isolation_filter)))


            density_string = ''
            if density_filter:
                density_string = 'density'
            else:
                density_string = 'no_density'

            p_shrinkers_label = p_shrinker
            if not filter_with_poisson:
                print('Changing p-shrinker label')
                p_shrinkers_label = 'no_poisson'
        
            merged_windows_after_poisson_and_score_filtering_and_isolation_filter.to_csv('{}/{}_filtered_merged_score_{}_poisson_boost_{}_{}.tsv'.format(output_folder, label,
                                                                                                                                                         score_threshold,
                                                                                                                                                         p_shrinkers_label,
                                                                                                                                                         density_string
                                                                                                                                                        ), sep='\t', header=True, index=False)