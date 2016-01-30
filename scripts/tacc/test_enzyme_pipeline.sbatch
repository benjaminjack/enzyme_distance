#!/bin/bash

#SBATCH -J test_pipeline             # Job name
#SBATCH -o test_pipeline.%j.out       # Name of stdout output file (%j expands to jobId)
#SBATCH -e test_pipeline.%j.e       # Name of stdout output file (%j expands to jobId)
#SBATCH --mail-user=b.jack@utexas.edu
#SBATCH --mail-type=FAIL
#SBATCH -A A-bio7          # <-- Allocation name to charge job against

#development
#SBATCH -n 16                 # Total number of mpi tasks requested (16 cores/node) (change to 128 or 256 for production)
#SBATCH -t 2:00:00           # Run time (hh:mm:ss)
#SBATCH -p development       # Queue name (change to normal for production)

TIME=/usr/bin/time

export MKL_MIC_ENABLE=1

module load blast
module load python/2.7.3-epd-7.3.2
module load mafft

echo Starting job running
cd $WORK/enzyme_distance_data/enzymes/
$HOME/projects/enzyme_distance/scripts/tacc/run_enzyme_pipeline.sh structures/132l.pdb 132L A


