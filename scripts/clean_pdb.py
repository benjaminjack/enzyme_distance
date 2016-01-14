__author__ = 'ben'

#############################################################################
# This script was written by Benjamin Jack in the laboratory of Dr. Claus   #
# Wilke at the University of Texas at Austin. For any questions or concerns #
# email Ben at benjamin.r.jack@gmail.com.                                   #
#############################################################################

import sys
import argparse
import os
from Bio.PDB import PDBIO, Select
from renumber_pdb import parsePDBStructure, renumberResidues


# Define a custom chain selector. This is probably not the most 'pythonic' way of doing this

class ChainSelect(Select):

    def __init__(self, chain_name):
        self.chain_name = chain_name

    def accept_chain(self, chain):
        if chain.get_id() == self.chain_name:
            return 1
        else:
            return 0

    def accept_model(self, model):
        # We only pull out model 0. This is fine for our data set, but keep in mind
        # that a true biological assembly may include multiple models.
        if model.get_id() == 0:
            return 1
        else:
            return 0

def clean_pdb(pdb_file, pdb_chain, out_dir):

    out_dir_chain = out_dir + '/' + 'chain'

    if not os.path.isfile(pdb_file):
        raise argparse.ArgumentTypeError("PDB file could not be found.")

    # Create output directories if they do not already exist
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(out_dir_chain):
        os.makedirs(out_dir_chain)

    # Grab PDB name and make sure it's converted to uppercase
    pdb_name = os.path.basename(pdb_file).split('.')[0].upper()

    # Extract chain of interest
    structure = parsePDBStructure(pdb_file)
    # Make sure chain exists, otherwise throw an error
    try:
        chain = structure[0][pdb_chain]
    except KeyError:
        print("\nERROR:\n\n\t"+pdb_name+": chain "+pdb_chain+" could not be found.\n")
        return

    io = PDBIO()
    chain_select = ChainSelect(pdb_chain)
    io.set_structure(structure)
    pdb_chain_file = out_dir_chain+'/'+pdb_name+'_'+pdb_chain+'_temp.pdb'
    io.save(pdb_chain_file, chain_select)

    # Remove HetAtoms
    temp_file = out_dir + "/" + pdb_name + "_temp.pdb"
    removeHetAtoms(pdb_file, temp_file)

    temp_file_chain = out_dir_chain + "/" + pdb_name + '_' + pdb_chain + "_temp2.pdb"
    removeHetAtoms(pdb_chain_file, temp_file_chain)

    # Renumber PDB
    structure = parsePDBStructure(temp_file)
    (new_pdb, renumbered_pdb) = renumberResidues(structure, pdb_name)

    structure_chain = parsePDBStructure(temp_file_chain)
    (new_pdb_chain, renumbered_pdb) = renumberResidues(structure_chain, pdb_name + '_' + pdb_chain)

    # Remove waters
    removeWaters(new_pdb, out_dir + "/" + pdb_name + ".pdb")
    removeWaters(new_pdb_chain, out_dir_chain+'/'+pdb_name+'_'+pdb_chain+'.pdb')

    # Clean up temporary files
    os.remove(pdb_chain_file)
    os.remove(temp_file)
    os.remove(temp_file_chain)
    os.remove(new_pdb_chain)
    os.remove(new_pdb)

def main():

    parser = argparse.ArgumentParser(description="Clean PDB and split into individual PDBs for each chain.")
    parser.add_argument("input_pdb", help="Input PDB file or list of PDB files for batch processing with PDB filename and one chain per line, separated by a space.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-c", "--chain", help="PDB chain to extract from PDB file.")
    group.add_argument("-b", "--batch", help="Batch process files supplied as PDB list.", action="store_true")
    parser.add_argument("-i", "--in_dir", default="./", help="Input directory for batch processing PDBs.")
    parser.add_argument("-o", "--out_dir", default="./clean", help="Directory to output cleaned PDBs and chains. Defaults to './clean'")
    args = parser.parse_args()

    pdb_file = args.input_pdb
    pdb_chain = args.chain
    out_dir = args.out_dir
    input_dir = args.in_dir

    if not args.batch:
        # Run single PDB file
        clean_pdb(input_dir + pdb_file, pdb_chain, out_dir)
    else:
        # Check to make sure PDB list file exists
        if not os.path.isfile(pdb_file):
            raise argparse.ArgumentTypeError("PDB list file could not be found.")

        with open(pdb_file) as f:
            for line in f.readlines():
                split_line = line.split()
                pdb_file2 = split_line[0]
                pdb_chain = split_line[1]
                clean_pdb(input_dir + pdb_file2, pdb_chain, out_dir)



def removeHetAtoms(infile, outfile):
    with open(infile, "rU") as f:
        with open(outfile, "w") as fout:
            for line in f:
                if "HETATM" not in line:
                    fout.write(line)
                    # Record line

def removeWaters(infile, outfile):
    with open(infile, "rU") as f:
        with open(outfile, "w") as fout:
            for line in f:
                if "HOH" not in line and "ATOM" in line:
                    fout.write(line)
                    # Record line

## Execute the main function
if __name__ == "__main__":
    main()
