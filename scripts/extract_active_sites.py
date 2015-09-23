#############################################################################
# This script was written by Benjamin Jack in the laboratory of Dr. Claus   #
# Wilke at the University of Texas at Austin. For any questions or concerns #
# Ben at benjamin.r.jack@gmail.com.                                         #
#############################################################################

import pandas as pd
import sys
import os
from Bio.PDB import PDBParser

resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }


def main(argv):

    usage = "\n\nUsage: extract_active_sites.py <active site db> <pdb name> <pdb file>\n\n"

    if len(argv) != 3:
        print(usage)
        sys.exit()

    active_sites_file = argv[0]
    active_site = pd.read_csv(active_sites_file)
    pdb = argv[1]
    pdb_file = argv[2]

    assert os.path.exists(pdb_file) and os.path.exists(active_sites_file), usage

    active_list = list()
    p = PDBParser()
    structure = p.get_structure('X', pdb_file)

    for model in structure:
        for chain in model:
            active2 = active_site[(active_site['PDB ID'] == pdb.lower()) &
                                  (active_site['CHAIN ID'] == chain.get_id()) &
                                  (active_site['EVIDENCE TYPE'] == 'LIT')]
            for residue in chain.get_residues():
                if residue.get_resname() not in resdict.keys():
                    print pdb.upper()+'_'+chain.get_id()+': '+str(residue.get_id()[1])+' '+residue.get_resname()+' skipped.'
                    continue
                if residue.get_id()[1] in active2['RESIDUE NUMBER'].values:
                    # Sanity check
                    active_res = active2.loc[active2['RESIDUE NUMBER'] == residue.get_id()[1], 'RESIDUE TYPE'].values[0]
                    if residue.get_resname().upper() == active_res.upper():
                        active_list.append({'PDB': pdb,
                                            'RES': resdict[residue.get_resname()],
                                            'SITE_NUMBER': residue.get_id()[1],
                                            'ACTIVE_SITE': 1})
                    else:
                        print residue.get_id()[1], residue.get_resname().upper(), active_res.upper()
                        raise Exception, "Residues do not match, numbering is not consistent with PDB."
                else:
                    active_list.append({'PDB': pdb.upper(),
                                        'RES': resdict[residue.get_resname()],
                                        'SITE_NUMBER': residue.get_id()[1],
                                        'ACTIVE_SITE': 0})

    output = pd.DataFrame(data=active_list)
    output.to_csv(pdb+'_act.csv', index=False)

if __name__ == '__main__':
    main(sys.argv[1:])
