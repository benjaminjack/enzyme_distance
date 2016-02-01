#!/bin/csh
#
# Simple SGE script for submitting multiple serial
# jobs (e.g. parametric studies) using a script wrapper
# to launch the jobs.
#
# To use, build the launcher executable and your
# serial application(s) and place them in your WORKDIR
# directory.  Then, edit the CONTROL_FILE to specify 
# each executable per process.
#-------------------------------------------------------
#-------------------------------------------------------
# 
#         <------ Setup Parameters ------>
#
#SBATCH -J calc_rates          # Job name
#SBATCH -N 10                   # Total number of nodes (16 cores/node)
#SBATCH -n 40                  # Total number of tasks
#SBATCH -c 4		       # Number of cores per task
#SBATCH -p normal              # Queue name
#SBATCH -o calc_rates.%j.o      # Name of stdout output file (%j expands to jobid)
#SBATCH -t 2:00:00            # Run time (hh:mm:ss)
#SBATCH -e calc_rates.%j.e       # Name of stdout output file (%j expands to jobId)
#SBATCH --mail-user=b.jack@utexas.edu
#SBATCH --mail-type=ALL

#      <------------ Account String ------------>
# <--- (Use this ONLY if you have MULTIPLE accounts) --->
##SBATCH -A A-bio7
#------------------------------------------------------

module load launcher
module load blast
module load python/2.7.3-epd-7.3.2
module load mafft
module load boost

setenv EXECUTABLE     $TACC_LAUNCHER_DIR/init_launcher 
setenv CONTROL_FILE   enzyme_command_list.sh
setenv WORKDIR        $WORK/enzyme_distance_data/enzymes/
# 
# Variable description:
#
#  EXECUTABLE     = full path to the job launcher executable
#  CONTROL_FILE   = text input file which specifies
#                   executable for each process
#                   (should be located in WORKDIR)
#  WORKDIR        = location of working directory
#
#      <------ End Setup Parameters ------>
#--------------------------------------------------------
#--------------------------------------------------------

#----------------
# Error Checking
#----------------

if ( ! -e $WORKDIR ) then
        echo " "
	echo "Error: unable to change to working directory."
	echo "       $WORKDIR"
	echo " "
	echo "Job not submitted."
	exit
endif

if ( ! -f $EXECUTABLE ) then
	echo " "
	echo "Error: unable to find launcher executable $EXECUTABLE."
	echo " "
	echo "Job not submitted."
	exit
endif

if ( ! -f $WORKDIR/$CONTROL_FILE ) then
	echo " "
	echo "Error: unable to find input control file $CONTROL_FILE."
	echo " "
	echo "Job not submitted."
	exit
endif


#----------------
# Job Submission
#----------------

cd $WORKDIR/
echo " WORKING DIR:   $WORKDIR/"

$TACC_LAUNCHER_DIR/paramrun $EXECUTABLE $CONTROL_FILE

echo " "
echo " Parameteric Job Complete"
echo " "
