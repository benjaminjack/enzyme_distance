#!/bin/bash

#SBATCH -J makeblastdb             # Job name
#SBATCH -o makeblastdb.%j.out       # Name of stdout output file (%j expands to jobId)
#SBATCH -e makeblastdb.%j.e       # Name of stdout output file (%j expands to jobId)
#SBATCH --mail-user=b.jack@utexas.edu
#SBATCH --mail-type=FAIL
#SBATCH -A A-bio7          # <-- Allocation name to charge job against

#development
#SBATCH -n 16                 # Total number of mpi tasks requested (16 cores/node) (change to 128 or 256 for production)
#SBATCH -t 48:00:00           # Run time (hh:mm:ss)
#SBATCH -p normal       # Queue name (change to normal for production)

TIME=/usr/bin/time

export MKL_MIC_ENABLE=1

module load blast

echo Starting job running
cd $WORK/blast_db/
makeblastdb -in uniref90.fasta -input_type fasta -dbtype prot -out uniref90 -parse_seqids

