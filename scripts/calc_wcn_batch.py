#!/usr/bin/python

import os
import argparse
import subprocess
import calc_wcn as calc_wcn

def main():
    parser = argparse.ArgumentParser(description="Calculate WCN values for a set of PDB files.")
    parser.add_argument("pdb_list", help="List of PDB files for batch processing, one file per line.")
    parser.add_argument("input_dir", help="Input directory containing PDB files.")
    parser.add_argument("out_dir", help="Output directory for WCN files.")
    args = parser.parse_args()

    if not os.path.isfile(args.pdb_list):
        raise argparse.ArgumentTypeError("PDB list file could not be found.")

    with open(args.pdb_list) as f:
        for line in f.readlines():
            pdb_file = line.strip()
            pdb_name = os.path.splitext(pdb_file)[0]
            command = "calc_wcn.py " + args.input_dir + "/" + pdb_file + " " + args.out_dir + "/" + pdb_name + "_wcn.csv"
            print(command)
            calc_wcn.main(['', args.input_dir + "/" + pdb_file, args.out_dir + "/" + pdb_name + "_wcn.csv"])

if __name__ == "__main__":
    main()
