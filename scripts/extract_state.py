import sys, subprocess, os
from Bio.PDB.Polypeptide import *
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from renumber_pdb import parsePDBStructure

###########################################################################
# This script takes a full PDB file and determines if it is multimeric (1)
# or monomeric (0).
###########################################################################

def main():
    args = sys.argv
    pdb = args[1]
    pdb_name = os.path.basename(pdb).split('.')[0].upper()

    usage = "\nUsage:\n\n\textract_state.py <pdb>\n\n" \
            "\tThis script prints whether a given PDB is monomeric or multimeric. Must be used with RAW biological assembly for accurate results."

    assert os.path.exists(pdb), usage

    structure = parsePDBStructure(pdb)
    print pdb_name,",",checkMultimer(structure)

    return 0


def checkMultimer(structure):
    i = 0
    for model in structure:
        for chain in model:
            i += 1

    if i > 1:
        return 1
    elif i == 1:
        return 0
    else:
        return "NA"

## Execute the main function
if __name__ == "__main__":
    main()

