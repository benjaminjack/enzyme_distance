__author__ = 'ben'

#############################################################################
# This script was written by Benjamin Jack in the laboratory of Dr. Claus   #
# Wilke at the University of Texas at Austin. For any questions or concerns #
# Ben at benjamin.r.jack@gmail.com.                                         #
#############################################################################

from Bio.PDB import *
import csv
import os
import sys
import protein_toolbox_helper as ph
import math

#This is a dictionary that relates the three letter amino acid abbreviation with its one letter abbreviation
resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }


def main(argv):

    usage = "\nUsage:\n\n\tcalc_distances.py <input file> <output directory>\n"

    if len(argv) != 2:
        print(usage)
        sys.exit()

    pdb_file = argv[0]
    out_dir = argv[1]

    assert os.path.exists(pdb_file), usage

    # Make output directory if it does not exist
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Grab PDB name and make sure it's converted to uppercase
    pdb_name = os.path.basename(pdb_file).split('.')[0].upper()

    calcDistances(pdb_name, pdb_file, out_dir)


def calcDistances(pdb_name, input_file, output_dir):
    p = PDBParser()
    structure = p.get_structure('X', input_file)

    resnam   = []     # A list containing all Residue Names
    resnum   = []     # A list containing all Residue Numbers
    reschain = []     # A list containing the chain name in which the residues lie
    crdCA   = []
    crdSC   = []
    sizeSC  = []

    for residue in structure.get_residues():
        resnam.append(residue.resname)
        resnum.append(residue.get_full_id()[3][1])
        reschain.append(residue.get_full_id()[2])
        noSC = True
        noCA = True
        rescrd_SC = []  # A list containing the coordinates of all side chain atoms of the current residue. Will be used to calculate the COM of the side chain.
        for atom in structure.get_atoms():
            # atom.name is equivalent to atom.get_id()
            if atom.parent.id == residue.id and atom.name == 'CA':
                noCA = False
                crdCA.append(atom.get_coord())

            elif atom.parent.id == residue.id and atom.name not in ['C','CA','O','N']:
                noSC = False
                rescrd_SC.append(atom.get_coord())

        if noCA:
            print('\nERROR: ', pdb_name, 'has a missing alpha carbon at residue ', resnum[-1],
                  '. Please manually check PDB file for extraneous structural information (e.g. partially resolved'
                  'residues, bound DNA or ligands, etc.). Distances cannot be calculated.\n')
            sys.exit()

        if noSC:
            print('Missing side chain in residue: ', resnum[-1], resnam[-1], ', possibly a GLY:', pdb_name)
            crdSC.append(crdCA[-1])
            sizeSC.append(0)
        else:
            # Calculate side chain properties:
            sizeSC.append(len(rescrd_SC))
            crdSC.append(sum(rescrd_SC)/float(sizeSC[-1]))

    single_letter_aa = ph.convert_to_one_letter_code(resnam)
    distSC = [['PDB_AA'] + [str(res.get_id()[1])+resdict[res.get_resname()] for res in structure.get_residues()]]
    for i in range(len(resnam)):
        dist_sc_per_aa = [single_letter_aa[i]]
        for j in range(len(single_letter_aa)):
            distSCi = 0.
            if crdSC[i][0] == 'NA' or crdSC[j][0] == 'NA':
                dist_sc_per_aa.append('NA')
                continue
            if i != j:
                distSCi = math.sqrt((crdSC[i][0]-crdSC[j][0])**2 +
                                    (crdSC[i][1]-crdSC[j][1])**2 +
                                    (crdSC[i][2]-crdSC[j][2])**2)
            dist_sc_per_aa.append(distSCi)
        distSC.append(dist_sc_per_aa)

    with open(output_dir+'/'+pdb_name+'_dist.csv', 'wb') as output_file:
        csv_writer = csv.writer(output_file, delimiter=',')
        csv_writer.writerows(distSC)


if __name__ == "__main__":
    main(sys.argv[1:])



