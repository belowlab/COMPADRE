
import sys, os, json, math, argparse, time

parser = argparse.ArgumentParser()
parser.add_argument('-f1', '--fid1', help = 'fid1', dest = 'fid1')
parser.add_argument('-f2', '--fid2', help = 'fid2', dest = 'fid2')
parser.add_argument('-p', '--primus_dir', help = 'primus_dir', dest = 'primus_dir')
parser.add_argument('-s', '--primus_sum', help = 'primus_sum', dest = 'primus_sum')
args = parser.parse_args()

# --------------------------------------

# not needed anymore
def do_we_need_ersa(fid1, fid2): 

    with open('../../../example_data/eur_small_ibdfile_KDE_likelihood_vectors_small', 'r') as f:
        for line in f:
            if line.split('\t')[0] == fid1 and line.split('\t')[2] == fid2:
                vectorlist = []
                for item in line.split('\t')[5].split(','):
                    vectorlist.append(float(item))
                ersa_sum_prob = 1 - sum(vectorlist[0:3])
                if sum(vectorlist[0:3]) >= 0.8:
                    return (False)
                else:
                    return (1 - sum(vectorlist[0:3]))

# not currently being used since we wanna query all degrees
def get_init_ersa_info(fid1, fid2, primus_dir):

    with open (f'{primus_dir}/example_data/ersa_eur.ersa', 'r') as ersafile:
        for line in ersafile:
            if line.split('\t')[0] == fid1 and line.split('\t')[1] == fid2:
                ersa_top_degree = line.split('\t')[3]
                return (ersa_top_degree) 

# --------------------------------------------------------------------------------------------

def sum_unrelateds(fid1, fid2, primus_dir):

    # open ersa file 2, get all relationship listings higher than 3rd degree, convert and sum 
    with open (f'{primus_dir}/example_data/ersa_eur_2.ersa', 'r') as f:

        unrelated_sum = 0
        counter = 0

        for line in f:
            if line.split('\t')[0] == fid1 and line.split('\t')[1] == fid2 and int(line.split('\t')[3]) > 3:

                maxlnl = float(line.split('\t')[4].strip('\n'))
                inverted_likelihood = math.exp(maxlnl)
                unrelated_sum += inverted_likelihood
                counter += 1

    f.close()

    return (unrelated_sum, counter)
    
# ---------------------------------------------------------------

def sum_upto9(fid1, fid2, primus_dir):

    # open ersa file 2, get all relationship listings higher than 3rd degree, convert and sum 
    with open (f'{primus_dir}/example_data/ersa_eur_2.ersa', 'r') as f:

        upto9_sum = 0
        counter = 0

        for line in f:
            if line.split('\t')[0] == fid1 and line.split('\t')[1] == fid2 and (3 < int(line.split('\t')[3]) <= 9):

                maxlnl = float(line.split('\t')[4].strip('\n'))
                inverted_likelihood = math.exp(maxlnl)
                upto9_sum += inverted_likelihood
                counter += 1

    f.close()

    return (upto9_sum, counter)

# ---------------------------------------------------------------

def sum_1st_degree(fid1, fid2, primus_dir):

    # open ersa file 2, get all relationship listings higher than 3rd degree, convert and sum 
    with open (f'{primus_dir}/example_data/ersa_eur_2.ersa', 'r') as f:

        likelihood_sum = 0
        counter = 0

        for line in f:
            if line.split('\t')[0] == fid1 and line.split('\t')[1] == fid2 and int(line.split('\t')[3]) == 1:

                maxlnl = float(line.split('\t')[4].strip('\n'))
                inverted_likelihood = math.exp(maxlnl)
                likelihood_sum += inverted_likelihood
                counter += 1

    f.close()

    return (likelihood_sum, counter)
    
# ---------------------------------------------------------------

def sum_2nd_degree(fid1, fid2, primus_dir):

    with open (f'{primus_dir}/example_data/ersa_eur_2.ersa', 'r') as f:

        likelihood_sum = 0
        counter = 0
        lines = 0

        for line in f:
            if line.split('\t')[0] == fid1 and line.split('\t')[1] == fid2 and line.split('\t')[3] == '2':
                maxlnl = float(line.split('\t')[4].strip('\n'))
                # maxlnl = maxlnl * -1
                inverted_likelihood = math.exp(maxlnl)
                likelihood_sum += inverted_likelihood
                counter += 1

    f.close()

    return (likelihood_sum, counter)

# ---------------------------------------------------------------

def sum_3rd_degree(fid1, fid2, primus_dir):

    with open (f'{primus_dir}/example_data/ersa_eur_2.ersa', 'r') as f:

        likelihood_sum = 0
        counter = 0

        for line in f:
            if line.split('\t')[0] == fid1 and line.split('\t')[1] == fid2 and line.split('\t')[3] == '3':
                maxlnl = float(line.split('\t')[4].strip('\n'))
                # maxlnl = maxlnl * -1
                inverted_likelihood = math.exp(maxlnl)
                likelihood_sum += inverted_likelihood
                counter += 1

    f.close()

    return (likelihood_sum, counter)

# ---------------------------------------------------------------

def runner(fid1, fid2, primus_sum, primus_dir):

    total_ersa_sum = float(1 - float(primus_sum))
    #print (f'Probability to be divided among ERSA values: {str(total_ersa_sum)}')
    top_degree = get_init_ersa_info(fid1, fid2, primus_dir) # not currently being used

    # these give us the summed component for each of the last 3 VECTOR pieces

    #start_time = time.time()

    first_info = sum_1st_degree(fid1, fid2, primus_dir)
    #print (f'Sum of 1st degree probabilities: {first_info[0]}\nNumber of relationships: {first_info[1]}')

    second_info = sum_2nd_degree(fid1, fid2, primus_dir)
    #print (f'Sum of 2nd degree probabilities: {second_info[0]}\nNumber of relationships: {second_info[1]}')

    third_info = sum_3rd_degree(fid1, fid2, primus_dir)
    #print (f'Sum of 3rd degree probabilities: {third_info[0]}\nNumber of relationships: {third_info[1]}')

    upto9_info = sum_upto9(fid1, fid2, primus_dir)

    unrelated_info = sum_unrelateds(fid1, fid2, primus_dir)
    #print (f'Sum of UN degree probabilities: {unrelated_info[0]}\nNumber of relationships: {unrelated_info[1]}\n')
    #print ("Elapsed time: " + str(round(float(time.time() - start_time))) + ' seconds')


    ##################################################################################################################
    # now, weight them accordingly based on their makeup relative to the rest of the ERSA weights
    # currently done by weighing likelihood sum of the degree over total summed likelihoods

    #first_prop = float(first_info[0] / (first_info[0] + second_info[0] + third_info[0] + unrelated_info[0]))
    #print (f"\n1st Degree Proportion: {first_prop}")

    #denominator = float(second_info[0] + third_info[0] + unrelated_info[0])
    denominator = float(second_info[0] + third_info[0] + upto9_info[0])

    second_prop = float(second_info[0] / denominator)
    #print (f"2nd Degree Proportion: {second_prop}")

    third_prop = float(third_info[0] / denominator)
    #print (f"3rd Degree Proportion: {third_prop}")

    upto9_prop = float(upto9_info[0] / denominator)

    unrelated_prop = float(unrelated_info[0] / denominator)
    #print (f"UN Degree Proportion: {unrelated_prop}")


    ##################################################################################################################
    # need to also take into account the amount of probability left over from PRIMUS indices 0,1,2

    second_vector_value = float(total_ersa_sum * second_prop)
    third_vector_value = float(total_ersa_sum * third_prop)
    upto9_vector_value = float(total_ersa_sum * upto9_prop)
    unrelated_vector_value = float(total_ersa_sum * unrelated_prop)

    # 3 print outputs that get captured by PRIMUS in PERL
    print (second_vector_value)
    print (third_vector_value)
    print (upto9_vector_value)
    #print (unrelated_vector_value)

# -----------------------------------------------

if __name__ == '__main__':

    #fid1, fid2 = 'R201111371', 'R239357796'
    if 'fid1' in args and "fid2" in args and "primus_sum" in args and "primus_dir" in args:
        runner(args.fid1, args.fid2, args.primus_sum, args.primus_dir)