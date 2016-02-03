import os
import argparse
import subprocess
import merge_df_by_alignment as merge_df

def main():
    parser = argparse.ArgumentParser(description="Merge structure information tables for a set of structures.")
    parser.add_argument("pdb_list", help="List of PDBs and chains, separated by a space, one PDB and chain per line.")
    parser.add_argument("out_dir", help="Output directory for merged files.")
    args = parser.parse_args()

    if not os.path.isfile(args.pdb_list):
        raise argparse.ArgumentTypeError("PDB list file could not be found.")

    with open(args.pdb_list) as f:
        for line in f.readlines():
            pdb = line.strip().split()
            pdb_name = pdb[0]
            pdb_chain = pdb[1]
            command = "Processing... " + pdb_name + " " + pdb_chain
            print(command)
            try:
                merge_df.main(['', pdb_name, pdb_chain, args.out_dir])
            except:
                print("PDB " + pdb[0] + " failed.")

if __name__ == "__main__":
    main()
