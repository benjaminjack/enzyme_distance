__author__ = 'ben'

import sys
import os
import subprocess
from random import sample
from Bio import AlignIO, SeqIO
from tempfile import NamedTemporaryFile

def run_mafft(input_file, output_file):
    call_mafft = "mafft --auto --inputorder --anysymbol " + input_file + " > " + output_file

    # Must redirect stderr because of mafft/tacc wonkiness. Mafft will not work otherwise.
    FNULL = open(os.devnull, 'w')
    mafft = subprocess.call(call_mafft, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    assert (mafft == 0), "Call to mafft failed."


def main():
    # Where is the fasta file containing aligned sequences?
    usage = "Usage: python downsample_seqs.py <fasta_alignment_file> <fasta_pdb_name> <sample_size>. fasta_file must contain at least one sequence. "
    assert(sys.argv[1] and sys.argv[2] and os.path.isfile(sys.argv[1])), usage
    rep = str(sys.argv[1])
    fasta_name = sys.argv[2]
    sample_size = int(sys.argv[3])

    outfile = fasta_name +'_sample.fasta'

    alignments = AlignIO.read(rep, "fasta")

    if len(alignments) > 300:
        downsample = sample(alignments[1:], sample_size)
        # Make sure ref seq is first in fasta file
        downsample.insert(0, alignments[0])
        # Write out unique seqs as temporary FASTA files
        with NamedTemporaryFile() as tempfile:
            SeqIO.write(downsample, tempfile, "fasta")
            tempfile.flush()
            run_mafft(tempfile.name, outfile)
    else:
        with open(outfile, 'w') as f:
            SeqIO.write(alignments, f, "fasta")


main()