#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#  These options are:
#  job name: -N
#  use the current working directory: -cwd
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem
#$ -N fix-sti1
#$ -pe sharedmem 16
#$ -t 1-80
#$ -R y
#$ -P roslin_wilson
#$ -cwd
#$ -o out
#$ -e out
#$ -l h_rt=2:00:00     # <-- ** make sure this is correct!!! **

# Initialise the environment modules
. /etc/profile.d/modules.sh

module load R/4.6
module load openmpi/5.0.7

unset R_HOME

# Run the program
Rscript tmp_fix_sti1_files.R $SGE_TASK_ID

