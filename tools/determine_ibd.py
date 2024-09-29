import os, sys
import pandas as pd
import numpy as np
import itertools

def remove_haplotype(id_str):
    return id_str.split('.')[0]

def sort_ids_numerically(id1, id2):
    num1 = int(id1[2:])
    num2 = int(id2[2:])
    return (id1, id2) if num1 < num2 else (id2, id1)

def find_shared_segments_new(pair_df):
    id1, id2 = sorted(set(pair_df['id1_clean'].unique()) | set(pair_df['id2_clean'].unique()))
    
    # Separate dataframes for each haplotype combination
    df_00 = pair_df[(pair_df['id1'].str.endswith('.0') & pair_df['id2'].str.endswith('.0'))]
    df_01 = pair_df[(pair_df['id1'].str.endswith('.0') & pair_df['id2'].str.endswith('.1'))]
    df_10 = pair_df[(pair_df['id1'].str.endswith('.1') & pair_df['id2'].str.endswith('.0'))]
    df_11 = pair_df[(pair_df['id1'].str.endswith('.1') & pair_df['id2'].str.endswith('.1'))]
    
    # Combine all segments
    all_segments = pd.concat([df_00, df_01, df_10, df_11])
    all_segments = all_segments.sort_values('start')
    
    # Find unique breakpoints
    breakpoints = sorted(set(all_segments['start']) | set(all_segments['end']))
    
    # Initialize list to store results
    results = []
    
    # Iterate through adjacent breakpoints
    for i in range(len(breakpoints) - 1):
        start, end = breakpoints[i], breakpoints[i+1]
        
        # Check which segments overlap this range
        overlapping_segments = all_segments[(all_segments['start'] <= start) & (all_segments['end'] >= end)]
        
        # Determine IBD status
        haplotype_combinations = set(zip(overlapping_segments['id1'].str[-1], overlapping_segments['id2'].str[-1]))
        
        if len(haplotype_combinations) == 2 and ({('0', '1'), ('1', '0')} == haplotype_combinations or {('1', '0'), ('0', '1')} == haplotype_combinations):
            ibd = 2
        elif len(haplotype_combinations) > 0:
            ibd = 1

        else:
            continue  # Skip segments with no overlap
        
        exact_match = overlapping_segments[(overlapping_segments['start'] == start) & (overlapping_segments['end'] == end)]
        if not exact_match.empty:
            segment_length = exact_match['length'].iloc[0]
        else:
            containing_segment = overlapping_segments.iloc[0]
            total_length = containing_segment['end'] - containing_segment['start']
            segment_length = containing_segment['length'] * (end - start) / total_length
            
        results.append({
            'id1': id1,
            'id2': id2,
            'chrom': pair_df['chrom'].iloc[0],
            'start': start,
            'end': end,
            'length': segment_length,
            'ibd': ibd
        })
    
    result_df = pd.DataFrame(results)
    return result_df

def run_transform(simdir, write_output=True):

    # Initialize an empty DF to concat pairwise chrom results 
    result_columns = {'id1': str,'id2': str,'chrom': int,'start': int,'end': int,'length': float,'ibd': int}
    result_df = pd.DataFrame(columns=result_columns.keys()).astype(result_columns)

    working_dir = f'{simdir}/ss/files' 
    output_file = f'{working_dir}/../ersa_input_final_haploid.txt'

    df = pd.read_csv(output_file, sep='\t', header=None, names=['id1', 'id2', 'start', 'end', 'length', 'chrom'])

    df['id1_clean'] = df['id1'].apply(remove_haplotype)
    df['id2_clean'] = df['id2'].apply(remove_haplotype)
    df = df[df['id1_clean'] != df['id2_clean']] # remove segments that are shared in same ID on different haplotypes

    id_pairs = set()
    for _, row in df.iterrows():
        pair = sort_ids_numerically(row['id1_clean'], row['id2_clean'])
        id_pairs.add(pair)

    id_pairs = sorted(list(id_pairs), key=lambda x: (int(x[0][2:]), int(x[1][2:])))

    for pair in id_pairs:
        id1, id2 = pair
        pair_df = df[((df['id1_clean'] == id1) & (df['id2_clean'] == id2))]

        for chrom in range(1,23):
            pair_chrom_df = pair_df[pair_df['chrom'] == chrom]

            if not pair_chrom_df.empty:
                shared_segments = find_shared_segments_new(pair_chrom_df)
                result_df = pd.concat([result_df, shared_segments], ignore_index=True)

    # When done, filter to only segments >5cM, and to IBD2
    result_df = result_df[(result_df['length'] >= 5.0) & (result_df['ibd'] == 2)]
    result_df = result_df[['id1','id2','chrom','start','end','length','ibd']]


    '''
    NEW: 
    - Read the baseline ersa input final in as a DF and set ibd equal to 1
    - Concat these 2 dfs 
    - Save as a new one 

    id1	id2	chrom	start	end	length	ibd

    '''

    # Read in normal germline output file
    original_df = pd.read_csv(f'{working_dir}/../ersa_input_final.txt', sep='\t', header=None, names=['id1', 'id2', 'start', 'end', 'length', 'chrom'])
    original_df = original_df[original_df['length'] >= 5.0]
    original_df['ibd'] = 1

    # Concat
    original_df = pd.concat([original_df, result_df], ignore_index=True)

    ### ALSO, we need to take the SECOND occurrence of the segment, not the first in the case of duplicates 

    mask = original_df.duplicated(subset=['id1', 'id2', 'start', 'end', 'chrom'], keep='last')
    df_filtered = original_df[~mask]


    if write_output == True:
        df_filtered.to_csv(f'{working_dir}/../ersa_input_final_haploid_ibd_NEWEST.txt', sep='\t', index=None) # Formerly *_haploid_ibd.txt
    else:
        return df_filtered


if __name__ == "__main__":
    
    simdir = sys.argv[1] # /data100t1/home/grahame/projects/compadre/unified-simulations/analysis-output/simulations/EUR/uniform3_size20_sim200/uniform3_size20_sim200-4
    run_transform(simdir)