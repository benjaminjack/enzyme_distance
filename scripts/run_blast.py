#! /usr/bin/env python

################################################################################
##  Script to collect Huang, et al homologous protein sequences using PSI-BLAST.
##
##  Adapted from script by Stephanie J. Spielman (stephanie.spielman@gmail.com)
################################################################################
import sys
import os
import subprocess

######### Input argument ##########


############ Blast Options ############
# Note: These parameters were set for the UniRef90 database. If using another database, you should probably play with these.
#######################################
e_value = 1e-10
num_iterations = 2
outfmt = '"6 sseqid pident qlen slen"'
max_target_seqs = 150000
length_error = 0.2
percent_identity = 30.0



def grab_blast_output(infile, outfile, length_error, percent_identity):
    ''' From the blast returns, create a unique set of ids, and write this set to outfile.'''
    id_set=set()
    f = open(infile, 'rU')
    for line in f:
        line = line.strip()
        acc = ''
        pident = float(line.split()[1])
        qlen = float(line.split()[2])
        slen = float(line.split()[3])
        diff = abs(slen - qlen)/(qlen) # error computation for length differences
        print diff, slen, qlen
        if pident >= percent_identity and diff <= length_error:
            acc = line.split()[0][9:]
            if len(acc) > 5:
                id_set.add(acc)
    outf = open(outfile, 'w')
    for id in id_set:
        outf.write(str(id) + '\n')
    outf.close()        
    return id_set


def run_blast(blast_in, database):

    pdb_id = os.path.splitext(os.path.basename(blast_in))[0]

    # Define output files
    blastout = "seed_"+str(pdb_id)+".out"
    idout = "blast_ids_"+str(pdb_id)+".txt"

    # Run the blast search.
    call_blast = "psiblast -db "+database+" -out "+str(blastout)+" -query "+str(blast_in)+" -num_iterations "+str(num_iterations)+" -evalue "+str(e_value)+" -outfmt "+str(outfmt)
    print call_blast
    blast = subprocess.call(call_blast, shell=True)
    assert(blast == 0), "Call to BLAST failed."

    # Collect NCBI ids that blast returns and save to file.
    print "getting id set"
    grab_blast_output(blastout, idout, length_error, percent_identity)

def main(argv):
    
    # Where is the PDB ID text file?
    usage = "Usage: python run_blast.py <fasta_file> <db_path>. fasta_file must contain protein sequence. "

    if len(argv) != 2:
        print(usage)
        sys.exit()

    input_file = argv[0]
    database = argv[1]

    assert(os.path.isfile(input_file)), usage

    run_blast(input_file, database)

if __name__ == "__main__":
    main(sys.argv[1:])

