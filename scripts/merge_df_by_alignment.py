#############################################################################
# This script was written by Benjamin Jack in the laboratory of Dr. Claus   #
# Wilke at the University of Texas at Austin. For any questions or concerns #
# Ben at benjamin.r.jack@gmail.com.                                         #
#############################################################################



import sys
import os
import subprocess
import pandas as pd
import protein_toolbox_helper as ph


def align_seqs_mafft(seqs, headers, out_fasta):
    in_file = "temp_seqs_12345.txt"
    ph.write_seq_to_file(seqs, headers, in_file)
    # Linsi is slower but produces shorter, less gappy alignments. This speeds up the processing of the resulting
    # dataframes in python so it's a worthwhile trade-off.
    align_command = "linsi " + in_file + " > " + out_fasta
    subprocess.call(align_command, shell=True)
    os.unlink(in_file)


def merge_df_by_aa(all_dfs, col_list, retain_gaps=False):
    align_outfile = "temp.txt"
    all_seqs = []
    all_headers = []

    # We'll use this later for looping through all the dataframes/files
    df_range = xrange(len(all_dfs))

    for i in df_range:
        # Read in CSV files as dataframes
        seq_list = all_dfs[i][col_list[i]].values
        print seq_list
        # Grab AA sequence from column
        if len(seq_list[0]) == 3:
            seq = ''.join(ph.convert_to_one_letter_code(seq_list))
        elif len(seq_list[0]) == 1:
            seq = ''.join(seq_list)
        else:
            raise "Could not recognized amino acid column in dataframe."

        # Append sequences and filenames to list for maaft processing
        all_seqs.append(seq)
        all_headers.append('>AA_' + str(i))

    # Align sequences with mafft
    align_seqs_mafft(all_seqs, all_headers, align_outfile)

    # Extract the sequences from the aligned file
    (seqs, headers) = ph.get_sequences(align_outfile)

    # Now go backwards and build a dataframe with the sequences aligned
    # Column headers are file names and rows are each AA in sequence.
    aligned_seqs = dict()
    for i in df_range:
        aligned_seqs[headers[i][1:]] = list(seqs[i])
    aligned_df = pd.DataFrame(data=aligned_seqs)

    # Initialize list for storing newly constructed dataframes which will eventually be merged
    # with aligned_df.
    partial_dfs = []

    # Loop through file list and append rows to a new dataframe corresponding to the aligned
    # sequences. This is where the core functionality of this script happens.

    for s in df_range:
        j = 0
        # How many rows in the original CSV?
        row_length = len(all_dfs[s].index)
        current_col = all_headers[s][1:]
        # Initialize dict to store new rows. We will convert this back to a dataframe later.
        # Building a dictionary of lists first and then converting to a DF, rather than appending rows
        # to an existing DF is MUCH, MUCH faster.
        new_df_dict = dict()
        for index, row in aligned_df[current_col].iteritems():
            # If the sequence in the alignment is not blank, and we haven't reached the end of
            # the original CSV sequence length...
            if j >= row_length:
                break
            if row != '-':
                # Append rows from original CSV to our new dictionary
                new_df_dict[index] = all_dfs[s].iloc[j, :].values
                j += 1

        # Append newly generated dataframe to a list to concatenate later
        new_df = pd.DataFrame.from_dict(new_df_dict, orient='index')
        # Add back in appropriate column names
        new_df.columns = all_dfs[s].columns
        partial_dfs.append(new_df)

    # Insert aligned_df to the beginning of the list of new dataframes
    partial_dfs.insert(0, aligned_df)
    # Concatenate dataframes together by column (i.e., like cbind in R)
    aligned_df = pd.concat(partial_dfs, axis=1)
    # Remove any rows containing NaN
    if retain_gaps is False:
        aligned_df = aligned_df.dropna()

    return aligned_df


def main():
    # List of input CSV files
    if len(sys.argv) != 4:
        print("\n\nUsage:\n\n     "+sys.argv[0]+" <input pdb id> <chain> <output directory>")
        return
    else:
        pdb = [sys.argv[1], sys.argv[2]]
        output_dir = sys.argv[3]
        print pdb

    # calculate_distances('2ENG_A_renumbered', "raw_pdbs/chains/", "distances_test/")

    print "Processing: " + pdb[0]
    rsa_mono = pd.read_csv('rsa_mono/' + pdb[0] + '_rsa.csv')
    # Rename RSA column
    rsa_mono = rsa_mono.rename(columns={'RSA': 'RSA_mono'})
    rsa_multi = pd.read_csv('rsa/' + pdb[0] + '_rsa.csv')
    rates = pd.read_csv('extracted_rates/' + pdb[0] + '_' + pdb[1] + '_r4s.csv')
    raw_rates = pd.read_csv('extracted_rates/' + pdb[0] + '_' + pdb[1] + '_r4s_raw.csv')
    # Rename raw rate column
    raw_rates = raw_rates.rename(columns={'rate': 'raw_rate'})
    wcn_multi = pd.read_csv('wcn/' + pdb[0] + '_wcn.csv')
    wcn_mono = pd.read_csv('wcn_mono/' + pdb[0] + '_wcn.csv')
    active_sites = pd.read_csv('active_sites/' + pdb[0] + '_act.csv')
    distances = pd.read_csv('distances/' + pdb[0] + '_' + pdb[1] + '_dist.csv')
    df = merge_df_by_aa([rsa_mono, rsa_multi, rates, raw_rates, wcn_multi, wcn_mono, active_sites, distances],
                         ['Amino_Acid', 'Amino_Acid', 'pdb_aa', 'pdb_aa', 'resnam', 'resnam', 'RES', 'PDB_AA'],
                         retain_gaps=False)


    with open(output_dir + '/' + pdb[0] + '_merged.csv', 'w') as f:
        df.to_csv(f, index=False)

if __name__ == "__main__":
    main()
