
# Import from wilke lab toolbox
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.Polypeptide import *
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO

def main(argv):
    pdb_file = argv[0]
    pdb_name = argv[1]
    chain_name = argv[2]
    outfile = argv[3]

    usage = "\nUsage:\n\n\textract_aa.py <pdb_file> <pdb_name> <chain> <output>\n\n" \
            "\tExtracts an amino acid sequence from a given PDB file and chain and outputs it as a fasta file."

    if len(argv) != 4:
        print(usage)
        sys.exit()

    assert os.path.exists(pdb_file), usage

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
    main(sys.argv[1:])


