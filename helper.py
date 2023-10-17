import pandas as pd
import os, math, sys

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

    #model_df_un = model_df[(model_df['degree_of_relatedness'] >= 5) & (model_df['degree_of_relatedness'] <= 9)]
    model_df_un = model_df[model_df['degree_of_relatedness'] >= 5] 
    model_df_un_sum = model_df_un['maxlnl2'].sum()

    # re-proportion them 
    total_prop = model_df_2d_sum + model_df_3d_sum + model_df_4d_sum + model_df_un_sum
    model_df_2d_prop = model_df_2d_sum/total_prop
    model_df_3d_prop = model_df_3d_sum/total_prop
    model_df_4d_prop = model_df_4d_sum/total_prop
    model_df_un_prop = model_df_un_sum/total_prop

    return (model_df_2d_prop, model_df_3d_prop, model_df_4d_prop, model_df_un_prop)



ersa = '/data100t1/home/grahame/projects/compadre/ersa/ersa2_py3/ersa.py'


vector_file = sys.argv[1].strip()
ersa_input_file = sys.argv[2].strip() # .match simulated file 


ersa_dir = ersa_input_file.split('/')[:-1]
ersa_dir = '/'.join(ersa_dir)


df = pd.read_csv(vector_file, sep='\t')
df['values'] = df['LIKELIHOOD_VECTOR'].str.split(',')
df['values'] = df['values'].apply(lambda x: [float(val) for val in x])
df = pd.concat([df, df['values'].apply(lambda x: pd.Series(x))], axis=1)
df = df.rename(columns={0: 'p0', 1: 'p1', 2: "p2", 3: "p3", 4: "p4", 5: "p5"})
df = df.drop(columns=['values'])


#########################

mdf = pd.read_csv(ersa_input_file, sep='\t', header=None)

#########################
# Isolate all potentially improve-able lines into their own df and iterate over them

df['p0-1'] = df['p0'] + df['p1']
subset = df[df['p0-1'] <= 0.8].reset_index(drop=True)
the_rest = df[df['p0-1'] > 0.8].reset_index(drop=True)

counter = 0

for x in range(len(subset)):

    id1 = subset.loc[x, 'IID1']
    id2 = subset.loc[x, 'IID2']

    # NEED TO CHECK IF THERE'S ANY PAIRWISE SHARING IN THE FIRST PLACE -- If not, ERSA won't help us, and we can skip it 
    match_subset = mdf[(mdf[0] == id1) & (mdf[1] == id2)]
    if len(match_subset) != 0:

        #print (f'{id1}, {id2}')

        ersa_outfile = f'{ersa_dir}/{id1}_{id2}_ersa'
        ersa_command = f'python3 {ersa} --single_pair={id1}:{id2} --segment_files={ersa_input_file} --model_output_file={ersa_outfile}.model --output_file={ersa_outfile}.out'
        os.system(ersa_command)

        # after ersa is run on that pair, do the same steps as in ersa_helper 
        output_model_file = ersa_outfile + '.model'

        # there is a chance that ERSA was not able to build results off a small # of segments, so we need to check if there's even results to use
        omf_df = pd.read_csv(output_model_file, sep='\t')
        if len(omf_df) != 0:
        
            ersa_props = calculate_ersa_props(output_model_file)

            # multiply the ERSA proportions by p0-2 to re-fit them back into the PRIMUS probability space
            prop01 = subset.loc[x, 'p0-1']
            prop2plus = 1-prop01 
            ersa_props_updated = tuple(x * prop2plus for x in ersa_props) # 2,3,4,UN

            # now we can update the subset values and ultimately LIKELIHOOD_VECTOR
            p0_old = subset.loc[x, 'p0']
            p1_old = subset.loc[x, 'p1']
            subset.loc[x, 'p2'] = ersa_props_updated[0]
            subset.loc[x, 'p3'] = ersa_props_updated[1]
            subset.loc[x, 'p4'] = ersa_props_updated[2]
            subset.loc[x, 'p5'] = ersa_props_updated[3]

            # update LIKELIHOOD_VECTOR
            subset.loc[x, 'LIKELIHOOD_VECTOR'] = f'{p0_old},{p1_old},{ersa_props_updated[0]},{ersa_props_updated[1]},{ersa_props_updated[2]},{ersa_props_updated[3]}'

            counter += 1


new_total_df = pd.concat([subset, the_rest]).reset_index(drop=True)
new_total_df = new_total_df.iloc[:, :-7]
new_total_df.to_csv(vector_file + '.updated', sep='\t', index=False)
