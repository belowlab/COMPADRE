import pandas as pd
import os, sys, json, math, time, re
from dotenv import dotenv_values
import numpy as np
config = dotenv_values('.env')

# Config needs to be a full path because it has to be irrespective of where PRIMUS is being run 

python3 = config['PYTHON3_EXE']
ersa = config['ERSA_EXE']

###############################################################

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

    model_df_un = model_df[(model_df['degree_of_relatedness'] >= 5) & (model_df['degree_of_relatedness'] <= 40)] 
    model_df_un_sum = model_df_un['maxlnl2'].sum()

    # re-proportion them 
    total_prop = model_df_2d_sum + model_df_3d_sum + model_df_4d_sum + model_df_un_sum
    model_df_2d_prop = model_df_2d_sum/total_prop
    model_df_3d_prop = model_df_3d_sum/total_prop
    model_df_4d_prop = model_df_4d_sum/total_prop
    model_df_un_prop = model_df_un_sum/total_prop

    return (model_df_2d_prop, model_df_3d_prop, model_df_4d_prop, model_df_un_prop)


'''
PRIMUS supplies the pairwise IDs, the probability vector, and the ERSA INPUT***, so we can just assume that we need to run ERSA every time 

IMPORTANT -- THE IDS NEED TO MATCH HOW THEY SHOW UP IN THE .MATCH FILE 
'''

id1 = sys.argv[1]
id2 = sys.argv[2]
vector_str = sys.argv[3].strip()
ersa_input_file = sys.argv[4].strip()

# 2/23/24 -- ersa match file should be good enough to grab directory to look for 

if 'Missing' not in id1 and 'Missing' not in id2:   # CHANGE THIS???? 

    
    vector_arr = [float(x) for x in vector_str.split(',')]
    # prop2plus = float(1 - (vector_arr[0] + vector_arr[1]))

    #print (f'\n{vector_str}\n')
    #print (f'Potential ERSA probability space : {prop2plus}\n')

    #if prop2plus > 0.4:

    ersa_dir = ersa_input_file.split('/')[:-1]
    ersa_dir = '/'.join(ersa_dir) + '/ersa2'
    #print(ersa_dir)
    # make an output dir for ersa pairwise results
    if not os.path.exists(ersa_dir):
        os.makedirs(ersa_dir, exist_ok=True)

    # match the id string structure to what's in the ersa match input
    # pattern = r'^(.*?_)id\d+'
    # prefix = ''
    # try:
    #     with open(ersa_input_file, 'r') as file:
    #         line = file.readline()
    #         columns = line.split('\t')
    #         match = re.match(pattern, columns[0])
    #         prefix = match.group(1)
    # except:
    #     pass

    # id1_temp = f'{prefix}{id1}'
    # id2_temp = f'{prefix}{id2}'

    id1_temp = id1.split('_')[-1]
    id2_temp = id2.split('_')[-1]

    ### WE CAN UPDATE THIS ^ NOW THAT THE ERSA IDS ARE ALL CORRECT 

    ersa_outfile = f'{ersa_dir}/output_new_{id1_temp}_{id2_temp}'
    ersa_command = f'{python3} {ersa} --single_pair={id1_temp}:{id2_temp} --segment_files={ersa_input_file} --model_output_file={ersa_outfile}.model --output_file={ersa_outfile}.out'
    print (ersa_command + '\n')

    #sys.exit()

    os.system(ersa_command)

    output_model_file = f'{ersa_outfile}.model'
    omf_df = pd.read_csv(output_model_file, sep='\t')

    if len(omf_df) == 0:
        print (vector_str) # return what was provided since it can't be improved 

    else:
        # function call 
        ersa_props = calculate_ersa_props(output_model_file)

        # multiply the ERSA proportions by p0-2 to re-fit them back into the PRIMUS probability space
        prop02 = 1 - (vector_arr[0] + vector_arr[1])
        ersa_props_updated = tuple(x * prop02 for x in ersa_props) # 2,3,4,UN

        updated_vector = f'{vector_arr[0]},{vector_arr[1]},{ersa_props_updated[0]},{ersa_props_updated[1]},{ersa_props_updated[2]},{ersa_props_updated[3]}'

        print (updated_vector)

else:
    print (vector_str)