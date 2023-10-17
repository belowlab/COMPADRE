import pandas as pd
import os, math, sys

'''Goal: take a vector of probabilities from PRIMUS + ID pair, run ERSA on that pair (with the match file path), perform updates to the vector, and return a new vector to Perl'''

def calculate_ersa_props(model_file):

    #model_file = '/data100t1/home/grahame/projects/compadre/primus-ersa-v2/output/sim64/ersa/sim64_id2_id6_ersa.model'
    model_df = pd.read_csv(model_file, sep='\t', dtype={'maxlnl': float, 'degree_of_relatedness': int})
    model_df['maxlnl2'] = model_df['maxlnl'].apply(lambda x: math.exp(x))
    #model_df

    # every line is the pair we want, so we can skip the original step of checking that the IDs match

    model_df_2d = model_df[model_df['degree_of_relatedness'] == 2] 
    model_df_2d_sum = model_df_2d['maxlnl2'].sum()

    model_df_3d = model_df[model_df['degree_of_relatedness'] == 3] 
    model_df_3d_sum = model_df_3d['maxlnl2'].sum()

    model_df_4d = model_df[model_df['degree_of_relatedness'] == 4] 
    model_df_4d_sum = model_df_4d['maxlnl2'].sum()

    model_df_un = model_df[(model_df['degree_of_relatedness'] >= 5) & (model_df['degree_of_relatedness'] <= 9)] 
    model_df_un_sum = model_df_un['maxlnl2'].sum()

    # re-proportion them 
    total_prop = model_df_2d_sum + model_df_3d_sum + model_df_4d_sum + model_df_un_sum
    print (f'Total proportion found: {total_prop}')
    model_df_2d_prop = model_df_2d_sum/total_prop
    model_df_3d_prop = model_df_3d_sum/total_prop
    model_df_4d_prop = model_df_4d_sum/total_prop
    model_df_un_prop = model_df_un_sum/total_prop

    print (model_df_2d_prop, model_df_3d_prop, model_df_4d_prop, model_df_un_prop)

    return (model_df_2d_prop, model_df_3d_prop, model_df_4d_prop, model_df_un_prop)



ersa = '/data100t1/home/grahame/projects/compadre/ersa/ersa2_py3/ersa.py'

id1 = sys.argv[1]
id2 = sys.argv[2]
vector_str = sys.argv[3].strip()
vector_arr = [float(x) for x in vector_str.split(',')]
print (vector_arr)

prop2plus = float(1 - vector_arr[0] + vector_arr[1])
print (f'Potential ERSA prop: {prop2plus}')
ersa_match_file = sys.argv[4].strip()

ersa_dir = ersa_match_file.split('/')[:-1]
ersa_dir = '/'.join(ersa_dir)

if prop2plus > 0.2:

    ersa_outfile = f'{ersa_dir}/{id1}_{id2}_ersa'
    ersa_command = f'python3 {ersa} --single_pair={id1}:{id2} --segment_files={ersa_match_file} --model_output_file={ersa_outfile}.model --output_file={ersa_outfile}.out'
    os.system(ersa_command)

    output_model_file = ersa_outfile + '.model'

    # there is a chance that ERSA was not able to build results off a small # of segments, so we need to check if there's even results to use
    omf_df = pd.read_csv(output_model_file, sep='\t')
    if len(omf_df) < 5:
        print (vector_str) # return what was provided since it can't be improved 

    else:
        ersa_props = calculate_ersa_props(output_model_file)
        ersa_props_new = [x*prop2plus for x in ersa_props]
        return_output = f'{vector_arr[0]},{vector_arr[1]},{ersa_props_new[0]},{ersa_props_new[1]},{ersa_props_new[2]},{ersa_props_new[3]}'
        print (return_output)

else:

    print (vector_str)