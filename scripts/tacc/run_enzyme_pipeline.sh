#!/bin/bash

echo "Positional Parameters"
echo '$0 = ' $0
echo '$1 = ' $1
echo '$2 = ' $2
echo '$3 = ' $3

PDB_FILE=$1
PDB=$2
CHAIN=$3

SCRIPTS="$HOME/projects/enzyme_distance/scripts"
BLAST_DB_PATH="/work/03284/bjack/blast_db/uniref90"
BLAST_DB_NAME="UniRef90"
BLAST_THREADS=6
RAXML_PATH="raxmlHPC-PTHREADS-SSE3"
RAXML_THREADS=6
R4S_PATH="rate4site"
ACTIVE_SITE_DB="CSA_2_0_121113.txt"

# Let's make all of our output directories
echo "Making output directories."
mkdir -p structures
mkdir -p structures/clean
mkdir -p structures/states
mkdir -p fasta
mkdir -p blast
mkdir -p full_alignments
mkdir -p alignments_uniq
mkdir -p alignments_300
mkdir -p trees
mkdir -p rates
mkdir -p extracted_rates
mkdir -p rsa_mono
mkdir -p rsa_multi
mkdir -p wcn_mono
mkdir -p wcn_multi
mkdir -p distances
mkdir -p active_sites

abort()
{
    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

trap 'abort' 0

set -e


echo "Processing... $PDB_FILE"

# CLEAN PDB
# echo "Cleaning PDB."
# python $SCRIPTS/clean_pdb.py structures/$PDB_FILE -c $CHAIN -o structures/clean/
# echo "PDB cleaned."

# EXTRACT MULTIMERIC STATE
# echo "Extracting multimeric state."
# python $SCRIPTS/extract_state.py structures/$PDB_FILE > structures/states/${PDB}_states.csv
# echo "Multimeric state determined."

# EXTRACT AMINO ACID SEQUENCE
# echo "Extracting amino acid sequence from PDB."
# python $SCRIPTS/extract_aa.py structures/$PDB_FILE -c $CHAIN -o ./fasta
# echo "Amino acid sequence extracted."

# EXTRACT ACTIVE SITE
cd active_sites/
python $SCRIPTS/extract_active_sites.py ../$ACTIVE_SITE_DB $PDB ../structures/$PDB_FILE
cd ../

# RUN BLAST
# Change directory for blast output
# echo "Collecting homologous sequences."
# cd blast/
python $SCRIPTS/run_blast.py ../fasta/${PDB}_${CHAIN}.fasta $BLAST_DB_PATH $BLAST_THREADS
# And change the directory back...
# cd ../
# echo "Homologous sequences collected."

# CONSTRUCT AND DOWNSAMPLE ALIGNMENTS
# echo "Building alignments."
# cd full_alignments
# python $SCRIPTS/get_blast_seqs.py ../blast/blast_ids_${PDB}_${CHAIN}.txt ../fasta/${PDB}_${CHAIN}.fasta ${PDB}_${CHAIN} $BLAST_DB_NAME $BLAST_DB_PATH
# cd ../alignments_uniq
# python $SCRIPTS/find_unique_seqs.py ../full_alignments/${PDB}_${CHAIN}_aln.fasta ${PDB}_${CHAIN}_uniq.fasta
# cd ../alignments_300
# python $SCRIPTS/downsample_seqs.py ../full_alignments/${PDB}_${CHAIN}_aln.fasta ${PDB}_${CHAIN} 300
# cd ../
# echo "Sequences aligned and downsampled."

# RUN RAXML
# echo "Constructing tree."
# cd trees
# python $SCRIPTS/run_raxml.py ../alignments_300/${PDB}_${CHAIN}_sample.fasta ${PDB}_${CHAIN} $RAXML_PATH $RAXML_THREADS
# cd ../
# echo "Tree constructed."

# RUN RATE4SITE AND EXTRACT RATES
# echo "Computing rates."
# python $SCRIPTS/run_r4s.py $R4S_PATH trees/phy/${PDB}_${CHAIN}.phy trees/RAxML_bestTree.${PDB}_${CHAIN}.txt rates/${PDB}_${CHAIN}_r4s.txt rates/${PDB}_${CHAIN}_r4s_raw.txt
# echo "Rates calculated. Extracting rate information."

# python $SCRIPTS/extract_unmapped_rate4site_rates.py rates/${PDB}_${CHAIN}_r4s.txt extracted_rates/${PDB}_${CHAIN}_r4s.csv
# python $SCRIPTS/extract_unmapped_rate4site_rates.py rates/${PDB}_${CHAIN}_r4s_raw.txt extracted_rates/${PDB}_${CHAIN}_r4s_raw.csv

# echo "Rates extracted. Script complete."

# Compute distance matrices
# python $SCRIPTS/calc_distances.py structures/clean/chain/${PDB}_${CHAIN}.pdb distances/

# Calculate RSAs
python $SCRIPTS/calc_rsa.py structures/clean/${PDB}.pdb rsa_multi/${PDB}_asa.csv rsa_multi/${PDB}_rsa.csv
python $SCRIPTS/calc_rsa.py structures/clean/chain/${PDB}_${CHAIN}.pdb rsa_mono/${PDB}_${CHAIN}_asa.csv rsa_mono/${PDB}_${CHAIN}_rsa.csv

# Calculate WCN
python $SCRIPTS/calc_wcn.py structures/clean/${PDB}.pdb wcn_multi/${PDB}_wcn.csv
python $SCRIPTS/calc_wcn.py structures/clean/chain/${PDB}_${CHAIN}.pdb wcn_mono/${PDB}_${CHAIN}_wcn.csv

# Done!
trap : 0

echo >&2 '
************
*** DONE ***
************
'
