import pandas as pd
import os, sys, json, math, time, re
import ersa 


'''
TODO:

- incorporate the file preprocessing as an initial step of the socket connection 
    - convert germline file into leaner dictionary
    - keep in system memory

- update ersa.py to iterate through dictionary, not regular g2 file 

'''


###############################################################

def calculate_ersa_props(model_df):

    model_df['degree_of_relatedness'] = model_df['degree_of_relatedness'].astype(int)
    model_df['maxlnl2'] = model_df['maxlnl'].apply(lambda x: math.exp(x))
    
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

    return model_df_2d_prop, model_df_3d_prop, model_df_4d_prop, model_df_un_prop

###############################################################

# Figure out how to put in the socket listener here (?)

id1 = sys.argv[1]
id2 = sys.argv[2]
vector_str = sys.argv[3].strip()
ersa_input_file = sys.argv[4].strip()

if 'Missing' not in id1 and 'Missing' not in id2:   # CHANGE THIS???? 

    ersa_dir = ersa_input_file.split('/')[:-1]
    ersa_dir = '/'.join(ersa_dir) + '/ersa2'
    if not os.path.exists(ersa_dir):
        os.makedirs(ersa_dir, exist_ok=True)

    id1_temp = id1.split('_')[-1]
    id2_temp = id2.split('_')[-1]
    ersa_outfile = f'{ersa_dir}/output_new_{id1_temp}_{id2_temp}'

    # OLD
    #ersa_command = f'{python3} {ersa} --single_pair={id1_temp}:{id2_temp} --segment_files={ersa_input_file} --model_output_file={ersa_outfile}.model --output_file={ersa_outfile}.out'
    #print (ersa_command + '\n')
    #os.system(ersa_command)

    # NEW
    ersa_options = {
        "single_pair": f"{id1_temp}:{id2_temp}",
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