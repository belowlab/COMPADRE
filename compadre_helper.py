import pandas as pd
import os, sys, json, math, time, re
import ersa 

'''
TODO:

- incorporate the file preprocessing as an initial step of the socket connection 
    - convert germline file into leaner dictionary
    - keep in system memory
    
'''

###############################################################

def calculate_ersa_props(model_df):

    model_df['degree_of_relatedness'] = model_df['degree_of_relatedness'].astype(int)
    model_df['maxlnl2'] = model_df['maxlnl'].apply(lambda x: math.exp(x))
    
    model_df_2d = model_df[model_df['degree_of_relatedness'] == 2] 
    model_df_2d_sum = model_df_2d['maxlnl2'].sum()
    model_df_3d = model_df[model_df['degree_of_relatedness'] == 3] 
    model_df_3d_sum = model_df_3d['maxlnl2'].sum()
    # model_df_4d = model_df[model_df['degree_of_relatedness'] == 4] 
    # model_df_4d_sum = model_df_4d['maxlnl2'].sum()
    # model_df_un = model_df[(model_df['degree_of_relatedness'] >= 5) & (model_df['degree_of_relatedness'] <= 40)] 
    # model_df_un_sum = model_df_un['maxlnl2'].sum()
    model_df_4d = model_df[(model_df['degree_of_relatedness'] >= 4) & (model_df['degree_of_relatedness'] < 10)] # updated 7.2.24
    model_df_4d_sum = model_df_4d['maxlnl2'].sum()
    model_df_un = model_df[(model_df['degree_of_relatedness'] >= 10) & (model_df['degree_of_relatedness'] <= 40)] 
    model_df_un_sum = model_df_un['maxlnl2'].sum()

    # re-proportion them 
    total_prop = model_df_2d_sum + model_df_3d_sum + model_df_4d_sum + model_df_un_sum
    model_df_2d_prop = model_df_2d_sum/total_prop
    model_df_3d_prop = model_df_3d_sum/total_prop
    model_df_4d_prop = model_df_4d_sum/total_prop
    model_df_un_prop = model_df_un_sum/total_prop

    return model_df_2d_prop, model_df_3d_prop, model_df_4d_prop, model_df_un_prop

###############################################################

# Figure out how to put in the socket listener here (?)

id1 = sys.argv[1]
id2 = sys.argv[2]
vector_str = sys.argv[3].strip()
ersa_input_file = sys.argv[4].strip()

id1_temp = id1.split('_')[-1] # just in case
id2_temp = id2.split('_')[-1]

segment_dict = {}
with open (ersa_input_file, 'r') as f:
    for line in f:
        ls = line.split('\t')
        iid1, iid2, start, end, cmlen, chrom = ls[0], ls[1], int(ls[2]), int(ls[3]), round(float(ls[4]), 2), int(ls[5].strip())
        key = f"{iid1}:{iid2}"
        value = (chrom, start, end, cmlen)
        if cmlen >= 5.0:
            if key not in segment_dict:
                segment_dict[key] = [value,]
            else:
                segment_dict[key] += [value,]

idk = f"{id1_temp}:{id2_temp}"

if idk in segment_dict.keys(): # Pairs with one genetically unrelated individual won't be in here because of how we filter 

    segment_obj = {idk : segment_dict[idk]}
    segment_obj = json.dumps(segment_obj)

    ersa_dir = ersa_input_file.split('/')[:-1]
    ersa_dir = '/'.join(ersa_dir) + '/ersa'
    if not os.path.exists(ersa_dir):
        os.makedirs(ersa_dir, exist_ok=True)

    ersa_outfile = f'{ersa_dir}/output_new_{id1_temp}_{id2_temp}'

    # OLD
    #ersa_command = f'{python3} {ersa} --single_pair={id1_temp}:{id2_temp} --segment_files={ersa_input_file} --model_output_file={ersa_outfile}.model --output_file={ersa_outfile}.out'
    #print (ersa_command + '\n')
    #os.system(ersa_command)

    # NEW
    ersa_options = {
        "single_pair": f"{id1_temp}:{id2_temp}",
        "segment_dict": segment_obj,
        "segment_files": ersa_input_file,
        "model_output_file": f"{ersa_outfile}.model",
        "output_file": f"{ersa_outfile}.out",
        "return_output": True,
        "write_output": False
    }
    output_model_df = ersa.runner(ersa_options)

    if len(output_model_df) == 0:
        print (vector_str) # return what was provided since it can't be improved 

    else:
        ersa_props = calculate_ersa_props(output_model_df)
        vector_arr = [float(x) for x in vector_str.split(',')]
        prop02 = 1 - (vector_arr[0] + vector_arr[1])
        ersa_props_updated = tuple(x * prop02 for x in ersa_props) # multiply the ERSA proportions by p0-2 to re-fit them back into the PRIMUS probability space
        updated_vector = f'{vector_arr[0]},{vector_arr[1]},{ersa_props_updated[0]},{ersa_props_updated[1]},{ersa_props_updated[2]},{ersa_props_updated[3]}'
        print (updated_vector) # returned to PRIMUS

else:
    print (vector_str)