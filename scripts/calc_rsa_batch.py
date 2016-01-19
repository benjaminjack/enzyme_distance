#!/usr/bin/python

import os
import argparse
import subprocess
import calc_rsa as calc_rsa

def main():
    parser = argparse.ArgumentParser(description="Calculate RSA values for a set of PDB files.")
    parser.add_argument("pdb_list", help="List of PDB files for batch processing, one file per line.")
    parser.add_argument("input_dir", help="Input directory containing PDB files.")
    parser.add_argument("asa_out", help="Output directory for ASA files.")
    parser.add_argument("rsa_out", help="Output directory for RSA files.")
    args = parser.parse_args()

    if not os.path.isfile(args.pdb_list):
        raise argparse.ArgumentTypeError("PDB list file could not be found.")

    with open(args.pdb_list) as f:
        for line in f.readlines():
            pdb_file = line.strip()
            pdb_name = os.path.splitext(pdb_file)[0]
            command = "calc_rsa.py " + args.input_dir + "/" + pdb_file + " " + args.asa_out + "/" + pdb_name + "_asa.csv" + args.rsa_out + "/" + pdb_name + "_rsa.csv"
            print(command)
            calc_rsa.main(['', args.input_dir + "/" + pdb_file, args.asa_out + "/" + pdb_name + "_asa.csv", args.rsa_out + "/" + pdb_name + "_rsa.csv"])

if __name__ == "__main__":
    main()
