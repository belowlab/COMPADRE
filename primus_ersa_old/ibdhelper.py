import csv, argparse, os

from logger import logger
log = logger().log

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--ibdfile', help = 'ibdfile', dest = 'ibdfile')
parser.add_argument('-f', '--famfile', help = 'famfile', dest = 'famfile')
args = parser.parse_args()

# -------------------------------------------------------------------------

def make_smaller_ibdfile(ibdFile, idList):

    # initialize smaller vector file 2
    if ibdFile[-1] == '/':
        newfile = ibdFile.strip('/') + '_small'
    else:
        newfile = ibdFile + '_small'

    # make this file if not exist yet
    if os.path.exists(newfile) == False:   
        open(newfile, "w").close

    outfile = csv.writer(open(newfile, 'w'), delimiter = '\t')
    outfile.writerow(['FID1'] + ['IID1'] + ['FID2'] + ['IID2'] + ['RT'] + ['EZ'] + ['Z0'] + ['Z1'] + ['Z2'] + ['PI_HAT'] + ['PHE'] + ['DST'] + ['PPC'] + ['RATIO'])     # headers
    #log('Headers written', 'success')

    # save list of IDS from famfile to memory
    idlist = []
    with open(idList, 'r') as famfile:
        for line in famfile:
            idlist.append(line.split(' ')[0])
    #log(f'ID list populated, length {len(idlist)}', 'success')

    with open(ibdFile, 'r') as ibdfile, open(newfile, 'a') as outfile2:
        for line in ibdfile:
            if 'FID1' not in line:
                f1 = line.split('  ')[1]
                f2 = line.split('  ')[3]
                if f1 in idlist and f2 in idlist:
                    outfile2.write(line)

    print (newfile) # what gets captured by PRIMUS PERL



if __name__ == "__main__":

    if "ibdfile" in args and "famfile" in args:
        make_smaller_ibdfile(args.ibdfile, args.famfile)