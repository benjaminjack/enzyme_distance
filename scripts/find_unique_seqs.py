#! /usr/bin/env python

##############################################################################
##  Script to find unique sequences from a set of aligned sequences. Accepts
##  an aligned fasta file as input, and returns a fasta file as
##  output. Must be re-aligned separately.
##
##  Written by Benjamin Jack (bjack913@gmail.com)
##############################################################################

import sys
import os
import numpy
import subprocess
from tempfile import NamedTemporaryFile
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

########### Options ##############
percent_diff = 10


def get_unique_seqs(alignment, percent_diff):
    '''
    Filter out sequences that are below a percentage different threshold.
    '''

    # Convert biopython alignment object to numpy array
    # First column is ID
    align_array = numpy.array([[rec.id]+list(rec) for rec in alignment], numpy.character)
    # Grab the number of rows
    align_array_len = len(align_array[0:,0])
    # Grab the reference sequence which will be the first row in the array
    ref_seq = align_array[0,0:]
    # Initialize output array with ref_seq as first entry
    final_array = numpy.array([ref_seq], numpy.character, ndmin=2)
    # Initialize list of output ids
    final_ids = []

    diff = float(percent_diff)/100.0

    # For each sequence in array...
    for i in range(align_array_len):
        add_seq = True

        for j in range(len(final_array[0:,0])):
            # For each position check if two sequences are different and also not '-'
            diff_logic = numpy.all([(align_array[i,1:] != final_array[j,1:]),
                                    (align_array[i,1:] != '-'),
                                    (final_array[j,1:] != '-')], axis=0)
            n_diff = float(numpy.sum(diff_logic))
            # Calculate all non '-' positions
            comp_logic = numpy.all([(align_array[i,1:] != '-'), (final_array[j,1:] != '-')], axis=0)
            n_comp = float(numpy.sum(comp_logic))

            # Check that the number of camparable sites of greater than 0 to avoid divide by 0 errors
            if (n_comp > 0):
                # Check if difference is below threshold
                if ((n_diff/n_comp) < diff) or (align_array[i,0] in final_ids):
                    add_seq = False
                    print align_array[i,0]
                    break

        if add_seq:
            final_ids.append(align_array[i,0])
            final_array = numpy.vstack([final_array, align_array[i,0:]])

    # Calculate number of sequences removed
    removed = align_array_len - len(final_array[0:,0])
    print str(removed) + " sequence(s) removed."

    # Initialize list for conversion back to biopython object
    fasta_aln = []

    for i in range(len(final_array[0:,0])):
        fasta_aln.append(SeqRecord(Seq(''.join(final_array[i,1:])), id=final_array[i,0], description=''))

    return fasta_aln


def run_mafft(input_file, output_file):
    call_mafft = "mafft --auto --inputorder --anysymbol " + input_file + " > " + output_file

    # Must redirect stderr because of mafft/tacc wonkiness. Mafft will not work otherwise.
    FNULL = open(os.devnull, 'w')
    mafft = subprocess.call(call_mafft, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    assert (mafft == 0), "Call to mafft failed."


def main():
    # Where is the fasta file containing aligned sequences?
    usage = "Usage: python find_unique_seqs.py <fasta_file>. fasta_file must contain at least one sequence. "
    rep = str(sys.argv[1])
    assert(os.path.isfile(rep)), usage

    alignment = AlignIO.read(rep, "fasta")
    unique_seqs = get_unique_seqs(alignment, percent_diff)

    # Write out unique seqs as temporary FASTA files
    outfile = os.path.splitext(os.path.basename(rep))[0][0:6]+'_uniq.fasta'
    with NamedTemporaryFile() as tempfile:
        SeqIO.write(unique_seqs, tempfile, "fasta")
        tempfile.flush()
        run_mafft(tempfile.name, outfile)


main()
