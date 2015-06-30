__author__ = 'ben'

#############################################################################
# This script was written by Benjamin Jack in the laboratory of Dr. Claus   #
# Wilke at the University of Texas at Austin. For any questions or concerns #
# Ben at benjamin.r.jack@gmail.com.                                         #
#############################################################################

import sys
import getopt
import os
from Bio.PDB import PDBIO
from renumber_pdb import parsePDBStructure, renumberResidues

def main(argv):

    usage = "\n\nUsage:\n\n\tclean_pdb.py <input file> <chain> <output directory>\n\n" \
            "Output: \n\n" \
            "\t<output_directory>/<PDB>.pdb - cleaned biological assembly\n" \
            "\t<output_directory>/chain/<PDB>_<chain>.pdb - cleaned single selected chain for rate and distance processing\n"

    if len(argv) != 3:
        print(usage)
        sys.exit()

    pdb_file = argv[0]
    pdb_chain = argv[1]
    out_dir = argv[2]

    out_dir_chain = out_dir + '/' + 'chain'

    assert os.path.exists(pdb_file), usage

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
        print("\nERROR:\n\n\t"+pdb_name+": chain "+pdb_chain+" could not be found.\n"+usage)
        sys.exit()

    io = PDBIO()
    io.set_structure(chain)
    pdb_chain_file = out_dir_chain+'/'+pdb_name+'_'+pdb_chain+'_temp.pdb'
    io.save(pdb_chain_file)

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
    main(sys.argv[1:])







