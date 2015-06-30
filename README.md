# Enzyme Analysis Pipeline
This repository contains all scripts and data required to reproduce the analysis in the following paper:


I have tried to automate as much of the pipeline as possible. However, since many portions of the pipeline would need to be run on a cluster to be practical, the pipeline is broken into several steps that need to be manually completed. 

## Part I: Data preprocessing
1.  Clean the data
    `scripts/clean_pdb.py <input file> <chain> <output directory>`
         
    Begin with biological assemblies for all proteins in the directory ~/data/structures/raw. Note, you will need to choose a single chain for calculating rates and some structural metrics. Generally speaking, this is the chain for which active site information is available. If you choose a chain with no active site information, it will be impossible to later calculate distance to the a catalytic residue and the entire protein will be excluded from analysis later in the pipeline. Here is an example command:
         
    `scripts/clean_pdb.py data/structures/raw/1a2t.pdb A data/structures/clean/`

2.  Manually clean the data (optional)
    
    PDBs often have bound DNA and ligands in the structure, or other components that are not part of the protein itself. It is very difficult to remove these in an automated fashion. If you run into errors later in the pipeline, you should first check the structure visually to make sure there are no nucleic acids, ligands, or partial residues in the structure.

3.  Extract multimeric state
    `scripts/extract_states.py <input file>`
    
    Prints the name of the PDB and its multimeric state (0 = single subunit, 1 = multiple subunits). Input file must be _raw_ PDB file, otherwise the state returned will not be accurate. Here is an example command:
    
    `scripts/extract_states.py data/structures/raw/1h7o.pdb > states.csv`

4. Extract catalytic residues

5. Extract polypeptide sequence

## Part II: Site-wise Evolutionary Rates

## Part III: Structural Metrics

1. Calculate distance matrices
    -__scripts/calc_distances.py \<input file\> \<output directory\>__
    
    These distances matrices will be used to calculate distance to the active site later on. The input file should be a _single cleaned chain_ from Step 1 of Data preprocessing. Here is an example command:
    
     _scripts/calc_distances.py data/structures/clean/chain/1A2T_A.pdb data/distances/_
     
2. Calculate relative solvent accessibility (RSA)
    -__scripts/calc_rsa.py \<input file\> \<output ASA\> \<output RSA\>
    
    Absolute solvent accessibilities (ASA) are output as well as RSA. Both output files are in the CSV format. We calculate both RSA with respect to the bioligical assembly and RSA with respect to a single subunit. This means that the above script will need to be run twice, for example:
    
    _scripts/calc_rsa.py structures/clean/1A2T.pdb rsa/12AT_asa.csv rsa/12AT_rsa.csv_
    _scripts/calc_rsa.py structures/clean/chain/1A2T_A.pdb rsa_mono/12AT_asa.csv rsa_mono/12AT_rsa.csv_
    
3. Calculate weighted contact number (WCN)
    -__scripts/calc_rsa.py \<input file\> \<output ASA\> \<output RSA\>__
    
    Again, we must calculate WCN for both the biological assembly and a single subunit. Here is an example command:
    
    _scripts/calc_wcn.py structures/clean/1A2T.pdb wcn/12AT_wcn.csv_
    _scripts/calc_wcn.py structures/clean/chain/1A2T_A.pdb wcn_mono/12AT_wcn.csv_

     
## Part IV: Analyses and Figures

