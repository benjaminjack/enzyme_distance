# Enzyme Analysis Pipeline
This repository contains all scripts and data required to reproduce the analysis in the following paper:


I have tried to automate as much of the pipeline as possible. However, since many portions of the pipeline would need to be run on a cluster to be practical, the pipeline is broken into several steps that need to be manually completed. See "run_pipeline.sh" for an example of how all of the scripts from parts I, II, and III can be stitched together for a single protein.

## Part 0: Dependencies

The following software tools are required and must be available in your $PATH for the pipeline:

- MAFFT
- PSI-BLAST
- RAxML
- Rate4Site

## Part I: Data preprocessing

1.  Clean the data
    
    `scripts/clean_pdb.py <input file> <chain> <output directory>`
         
    Begin with biological assemblies for all proteins in the directory `~/data/structures/raw`. Note, you will need to choose a single chain for calculating rates and some structural metrics. Generally speaking, this is the chain for which active site information is available. If you choose a chain with no active site information, it will be impossible to later calculate distance to the a catalytic residue and the entire protein will be excluded from analysis later in the pipeline. Here is an example command:
         
    `scripts/clean_pdb.py data/structures/raw/1a2t.pdb A data/structures/clean/`

2.  Manually clean the data (optional)
    
    PDBs often have bound DNA and ligands in the structure, or other components that are not part of the protein itself. It is very difficult to remove these in an automated fashion. If you run into errors later in the pipeline, you should first check the structure visually to make sure there are no nucleic acids, ligands, or partial residues in the structure.

3.  Extract multimeric state

    `scripts/extract_state.py <input file>`
    
    Prints the name of the PDB and its multimeric state (0 = single subunit, 1 = multiple subunits). Input file must be _raw_ PDB file, otherwise the state returned will not be accurate. Here is an example command:
    
    `scripts/extract_state.py data/structures/raw/1h7o.pdb > states.csv`

4. Extract catalytic residues
   
   `scripts/extract_active_sites.py <active site db> <pdb name> <pdb file>`
   
   Extracts the catalytic residues from a text file of the CSA database and places them in a separate CSV. The database file is available to download [here](https://www.ebi.ac.uk/thornton-srv/databases/CSA/Downloads.php). Make sure you download the 'flat file' for the script to work correctly. 
    
   `scripts/extract_active_sites.py data/CSA_2_0_121113.txt 1A2T data/structures/raw/1a2t.pdb`
   
5. Extract polypeptide sequence
   
   `scripts/extract_aa.py <pdb_file> <pdb_name> <chain> <output_file>`
   
   Extract the amino-acid sequence as a fasta file for calculation of evolutionary rates. Input should be a _raw_ PDB file.
   
   `scripts/extract_aa.py data/structures/raw/12as.pdb 12AS A data/fasta/12AS_A.fasta`

## Part II: Site-wise Evolutionary Rates

1. Run PSI-BLAST

   `scripts/run_blast.py <fasta_file> <db_path or just db_name if using local PSI-BLAST set up>`
   
   Takes a fasta amino-acid sequence and runs a PSI-BLAST query on a local database. The parameters are optimized for the UniRef90 database. If you are using another database, then you may want to edit this script and change some of the parameters. The output is a list of homologous sequence IDs.
   
   `scripts/run_blast.py data/fasta/12AS_A.fasta UniRef90`
   
2. Construct alignments
   
   `scripts/get_blast_seqs.py <seq_id_file> <ref_seq_file> <ref_seq_id> <db_name> <db_path (optional)>`
   
   Takes a list of sequence IDs from Step 1 and pulls the actual sequences from the local database to create an alignment. You must use the same database as in Step 1!
   
   `scripts/get_blast_seqs.py data/blast/blast_ids_12AS_A.txt data/fasta/12AS_A.fasta 12AS_A UniRef90`
   
3. Downsample alignments

   `scripts/downsample_seqs.py <fasta_alignment_file> <fasta_pdb_name> <sample_size>`
   
   Downsample alignments to a maximum of 300 sequences. This is very important! Rate4Site (step 5) will not process alignments larger than 300 sequences without memory errors.
   
   `scripts/downsample_seqs.py data/full_alignments/12AS_A_aln.fasta 12AS_A 300`

4. Build trees with RAxML
   
   `scripts/run_raxml.py <fasta_alignment_file> <pdb_name> <raxml_version or raxml_path> <threads>`
   
   Build phylogenetic trees from sequence alignment. This step will be very slow, so I suggest using a multi-threaded version of RAxML.
   
   `scripts/run_raxml.py data/alignments_300/12AS_A_sample.fasta 12AS_A raxmlHPC-PTHREADS-SSE3 2`
   
5. Run Rate4Site
   
   `scripts/run_r4s.py <r4s_path> <sequence_file> <tree_file> <outfile> <outfile_raw>`
   
   Compute rates with Rate4Site software. The sequence file must be in PHY format or Rate4Site will give errors.
   
   `scripts/run_r4s.py rate4site trees/phy/12AS_A.phy trees/RAxML_bestTree.12AS_A.txt rates/12AS_A_r4s.txt rates/12AS_A_r4s_raw.txt`
   
6. Extract rates from Rate4Site output
   
   `scripts/extract_unmapped_rate4site_rates.py <input rate4site file> <output file>`
   
   Extract rates from Rate4Site output files and convert to CSV. Note that you'll want to run this for both the normalized and raw rates individually.
   
   `scripts/extract_unmapped_rate4site_rates.py rates/12AS_A_r4s.txt extracted_rates/12AS_A_r4s.csv`
   
   
## Part III: Structural Metrics

1.  Calculate distance matrices
    
    `scripts/calc_distances.py <input file> <output directory>`
    
    These distances matrices will be used to calculate distance to the active site later on. The input file should be a _single cleaned chain_ from Step 1 of Data preprocessing. Here is an example command:
    
    `scripts/calc_distances.py data/structures/clean/chain/1A2T_A.pdb data/distances/`
     
2.  Calculate relative solvent accessibility (RSA)
    
    `scripts/calc_rsa.py <input file> <output ASA> <output RSA>`
    
    Absolute solvent accessibilities (ASA) are output as well as RSA. Both output files are in the CSV format. We calculate both RSA with respect to the bioligical assembly and RSA with respect to a single subunit. This means that the above script will need to be run twice, for example:
    ```
    scripts/calc_rsa.py structures/clean/1A2T.pdb rsa/12AT_asa.csv rsa/12AT_rsa.csv
    scripts/calc_rsa.py structures/clean/chain/1A2T_A.pdb rsa_mono/12AT_asa.csv rsa_mono/12AT_rsa.csv
    ```
    
3.  Calculate weighted contact number (WCN)
    
    `scripts/calc_wcn.py <input file> <output file>`
    
    Again, we must calculate WCN for both the biological assembly and a single subunit. Here is an example command:
    ```
    scripts/calc_wcn.py structures/clean/1A2T.pdb wcn/12AT_wcn.csv
    scripts/calc_wcn.py structures/clean/chain/1A2T_A.pdb wcn_mono/12AT_wcn.csv
    ```
     
## Part IV: Analyses and Figures

