
# Import from wilke lab toolbox
import sys
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.Polypeptide import *
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO

def main():

    parser = argparse.ArgumentParser(description="Extract amino acid sequence of specific chain from PDB file.")
    parser.add_argument("input_pdb", help="Input PDB file or list of PDB files for batch processing with PDB filename and one chain per line, separated by a space.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-c", "--chain", help="PDB chain to extract from PDB file.")
    group.add_argument("-b", "--batch", help="Batch process files supplied as PDB list.", action="store_true")
    parser.add_argument("-i", "--in_dir", default="./", help="Input directory for batch processing PDBs.")
    parser.add_argument("-o", "--out_dir", default="./fasta", help="Directory to output cleaned PDBs and chains. Defaults to './clean'")
    args = parser.parse_args()

    pdb_file = args.input_pdb
    chain_name = args.chain
    out_dir = args.out_dir
    input_dir = args.in_dir

    if not args.batch:
        # Run single PDB file
        pdb_name = os.path.basename(pdb_file).split('.')[0].upper()
        outfile = pdb_name + "_" + chain_name + ".fasta"
        extract_aa(input_dir + pdb_file, pdb_name, chain_name, outfile)
    else:
        # Check to make sure PDB list file exists
        if not os.path.isfile(pdb_file):
            raise argparse.ArgumentTypeError("PDB list file could not be found.")

        with open(pdb_file) as f:
            for line in f.readlines():
                split_line = line.split()
                pdb_file2 = split_line[0]
                chain_name = split_line[1]
                pdb_name = os.path.basename(pdb_file2).split('.')[0].upper()
                outfile = pdb_name + "_" + chain_name + ".fasta"
                extract_aa(input_dir + pdb_file2, pdb_name, chain_name, outfile)


def extract_aa(pdb_file, pdb_name, chain_name, outfile):

    if not os.path.isfile(pdb_file):
        raise argparse.ArgumentTypeError("PDB file could not be found.")

    structure = parsePDBStructure(pdb_file)
    chain = structure[0][chain_name]
    aa_seq = get_aa_fromPDB(chain)
    fasta_seq_rec = SeqRecord(Seq(str(aa_seq)), id=pdb_name+'_'+chain_name, description='')
    with open(outfile, 'w') as out:
        SeqIO.write(fasta_seq_rec, out, "fasta")

def parsePDBStructure( pdb_id ):
    parser = PDBParser()
    structure = parser.get_structure('test_rsa', pdb_id)
    return structure

def get_aa_fromPDB( chain ):
    polypeptides = ""
    ppb=PPBuilder()
    for pp in ppb.build_peptides(chain):
        polypeptides += pp.get_sequence()
    return polypeptides

if __name__ == "__main__":
    main()
