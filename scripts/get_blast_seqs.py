#! /usr/bin/env python

##############################################################################
##  Script to convert list of sequence IDs into a multi-FASTA sequence file by
##  by querying local database and aligning with mafft.
##
##  Written by Benjamin Jack (bjack913@gmail.com)
##############################################################################
import sys
import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

########### blastdbcmd options ##############
# database = 'UniRef90'
# database = "$WORK/db/uniref90"
outfmt = '%f'


def get_seq_list(filename):
    '''
    Return a list of IDs read from a text file
    '''
    seq_id_list = []
    with open(filename, 'rU') as f:
        for line in f:
            seq_id_list.append(line.strip())

    return seq_id_list


def run_mafft(input_file, output_file):
    call_mafft = "mafft --auto --inputorder --anysymbol " + input_file + " > " + output_file

    # Must redirect stderr because of mafft/tacc wonkiness. Mafft will not work otherwise.
    FNULL = open(os.devnull, 'w')
    mafft = subprocess.call(call_mafft, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    assert (mafft == 0), "Call to mafft failed."


def call_blastdbcmd(filename, ref_seq_id, ref_seq_file, database, db_name):
    '''
    Query local database with list of IDs and generate a multi-FASTA sequence file
    '''
    seq_id_list = get_seq_list(filename)
    # entry_list = "'"+ref_seq_id+"'"
    entry_list = ''
    for seq_id in seq_id_list:
        entry_list += ",'"+db_name+"_"+seq_id+"'"
    # Write to a file
    call_blastdbcmd = "blastdbcmd -db "+database+" -dbtype prot -out "+ref_seq_id+"_seqs.fasta -outfmt "+outfmt+" -entry "+entry_list[1:]
    print call_blastdbcmd

    blastdbcmd = subprocess.call(call_blastdbcmd, shell=True)
    assert(blastdbcmd == 0), "Call to BLASTDBCMD failed."

    my_records = []

    with open(ref_seq_id+"_seqs.fasta", 'rU') as f:
        for record in SeqIO.parse(f, 'fasta'):
            # Switch selenocysteine (U) with cysteine (C)
            record.seq = Seq(str(record.seq).replace('U', 'C'), IUPAC.protein)
            my_records.append(record)

    with open(ref_seq_file) as ref:
        my_records.insert(0, SeqIO.read(ref, 'fasta'))

    SeqIO.write(my_records, ref_seq_id+"_seqs.fasta", "fasta")


def main(argv):
    # Where is the text file of sequence IDs?
    usage = "Usage: python get_blast_seqs.py <seq_id_file> <ref_seq_file> <ref_seq_id> <db_name> <db_path>\n\n" \
            "seq_id_file must contain at least one sequence ID and fasta ref sequence must be available. "
    if len(argv) < 4 or len(argv) > 5:
        print(usage)
        sys.exit()
    blast_seqs = argv[0]
    ref_seq_file = argv[1]
    ref_seq_id = argv[2]
    db_name = argv[3]

    if len(argv) == 5:
        db_path = argv[4]
    else:
        db_path = db_name


    # Pull the PDB ID from the filename and grab fasta file from /raw_pdbs directory
    assert(os.path.isfile(ref_seq_file) and os.path.isfile(blast_seqs)), usage

    # Query database
    call_blastdbcmd(blast_seqs, ref_seq_id, ref_seq_file, db_path, db_name)

    # Align the returned sequences
    run_mafft(ref_seq_id+"_seqs.fasta", ref_seq_id+"_aln.fasta" )


if __name__ == "__main__":
    main(sys.argv[1:])
