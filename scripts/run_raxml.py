#! /usr/bin/env python

##############################################################################
##  Script to build phylogenetic tree from an aligned multi-fasta file. Uses
##  RAxML.
##
##  Written by Benjamin Jack (bjack913@gmail.com)
##############################################################################

import sys
import os
import subprocess
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

########### Options ##############
# raxml_version = 'raxmlHPC-PTHREADS-SSE3'
phy_dir = 'phy'
# threads = '12'


def run_raxml(version, threads, phy_file, infile):

    outfile_ext = infile+'.txt'
    call_raxml = version+" -T "+str(threads)+" -s "+phy_dir+"/"+str(phy_file)+" -n "+outfile_ext+" -m PROTCATLG -p 12345"
    print call_raxml
    raxml = subprocess.call(call_raxml, shell=True)
    assert(raxml == 0), "Call to RAxML failed."


def main():
    # Where is the fasta file containing aligned sequences?
    usage = "Usage: python run_raxml.py <fasta_alignment_file> <pdb_name> <raxml_version or raxml_path> <threads>. fasta_file must contain at least one sequence. "

    if len(sys.argv) != 5:
        print(usage)
        sys.exit()

    rep = str(sys.argv[1])
    assert(os.path.isfile(rep)), usage

    pdb_name = sys.argv[2]
    raxml_version = sys.argv[3]
    threads = sys.argv[4]

    if not os.path.exists(phy_dir):
        os.makedirs(phy_dir)

    # Convert fasta to phy for RAxML
    phy_name = pdb_name+'.phy'
    alignments = AlignIO.read(rep, "fasta")
    with open(phy_dir+'/'+phy_name, 'w') as f:
        AlignIO.write(alignments, f, "phylip-relaxed")

    run_raxml(raxml_version, threads, phy_name, pdb_name)


main()