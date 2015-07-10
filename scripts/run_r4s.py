'''
rate4site -s phy/unique_test.phy -t RAxML_bestTree.test.txt -o r4s_test.res -Ml
'''

##############################################################################
##  Script to generate sitewise evolutionary rates from phylogenetic tree
##  using Rate4Site.
##
##  Written by Benjamin Jack (bjack913@gmail.com)
##############################################################################

import sys
import os
import subprocess

def run_r4s(rate4site, seqfile, treefile, outfile, outfile_raw):
    call_r4s = rate4site+' -s '+seqfile+' -t '+treefile+' -o '+outfile+' -y '+outfile_raw+' -Ml'
    print call_r4s
    with open(outfile+'.log', 'w') as f:
        r4s = subprocess.call(call_r4s, shell=True, stderr=f, stdout=f)
    assert(r4s == 0), "Call to Rate4Site failed: "+seqfile


def main(argv):
    usage = "Usage: python run_r4s.py <r4s_path> <sequence_file> <tree_file> <outfile> <outfile_raw>. Requires both a tree file and sequence file. "
    if len(argv) != 5:
        print(usage)
        sys.exit()
    rate4site = argv[0]
    seqfile = argv[1]
    treefile = argv[2]
    outfile = argv[3]
    outfile_raw = argv[4]
    assert(os.path.isfile(seqfile) and os.path.isfile(treefile)), usage
    run_r4s(rate4site, seqfile, treefile, outfile, outfile_raw)


if __name__ == '__main__':
    main(sys.argv[1:])