__author__ = 'ben'

#############################################################################
# This script was written by Benjamin Jack in the laboratory of Dr. Claus   #
# Wilke at the University of Texas at Austin. For any questions or concerns #
# Ben at benjamin.r.jack@gmail.com.                                         #
#############################################################################

import sys
import getopt
import os
from renumber_pdb import parsePDBStructure, renumberResidues

def main(argv):

    usage = "\nUsage:\n\n\tclean_pdb.py <input file> <output directory>\n"

    if len(argv) != 2:
        print(usage)
        sys.exit()

    pdb_file = argv[0]
    out_dir = argv[1]

    assert os.path.exists(pdb_file) and os.path.exists(out_dir), usage

    # Grab PDB name and make sure it's converted to uppercase
    pdb_name = os.path.basename(pdb_file).split('.')[0].upper()

    # Remove HetAtoms
    temp_file = out_dir + "/" + pdb_name + "_temp.pdb"
    removeHetAtoms(pdb_file, temp_file)

    # Renumber PDB
    structure = parsePDBStructure(temp_file)
    (new_pdb, chain_renumbered_pdb) = renumberResidues(structure, pdb_name)

    # Remove waters
    removeWaters(new_pdb, out_dir + "/" + pdb_name + ".pdb")

    # Clean up temporary files
    os.remove(temp_file)
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







