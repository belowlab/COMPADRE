import csv, argparse, os

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vector', help = 'vector', dest = 'vector')
parser.add_argument('-f', '--famfile', help = 'famfile', dest = 'famfile')
args = parser.parse_args()

# -------------------------------------------------------------------------

def isolate_subset_from_vectorfile(vectorFile, idList):

    # initialize smaller vector file 2
    if vectorFile[-1] == '/':
        newfile = vectorFile.strip('/') + '_small'
    else:
        newfile = vectorFile + '_small'

    # make this file if not exist yet
    if os.path.exists(newfile) == False:   
        open(newfile, "w").close

    # initialize new file
    outfile = csv.writer(open(newfile, 'w'), delimiter = '\t')
    outfile.writerow(['FID1'] + ['IID1'] + ['FID2'] + ['IID2'] + ['MOST_LIKELY_REL'] + ['LIKELIHOOD_VECTOR'] + ['IBD0'] + ['IBD1'] + ['IBD2'] + ['PI_HAT'] + ['-1'] + ['POSSIBLE_RELS'] + ['LIKELIHOOD_CUTOFF'])     # headers
    #outfile.close()

    # save list of IDS from famfile to memory
    idlist = []
    with open(idList, 'r') as famfile:
        for line in famfile:
            idlist.append(line.split(' ')[0])
    #log(f'ID list populated, length {len(idlist)}', 'success')
    
    # iterate through BIG vector file and only add lines whose IDs correspond to ones that are in the idlist
    with open(vectorFile, 'r') as vectorfile, open(newfile, 'a') as outfile2:
        for line in vectorfile:
            if line.split('\t')[0] != "FID1":
                if line.split('\t')[0] in idlist and line.split('\t')[2] in idlist:
                    # print (line.split('\t'))
                    outfile2.write(line)
                    # sys.exit()
    #log('Done!', 'success')

    print (newfile) # what gets captured by PRIMUS PERL



if __name__ == "__main__":

    if "vector" in args and "famfile" in args:
        isolate_subset_from_vectorfile(args.vector, args.famfile)