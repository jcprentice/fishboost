#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#  These options are:
#  job name: -N
#  use the current working directory: -cwd
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem
#$ -N fb-test7
#$ -pe sharedmem 16
#$ -t 1-6
#$ -R y
#$ -P roslin_wilson
#$ -cwd
#$ -o out
#$ -e out
#$ -l h_rt=600:00:00     # <-- ** make sure this is correct!!! **

# Initialise the environment modules
. /etc/profile.d/modules.sh

module load R/4.5
module load openmpi/5.0.7

# Run the program
Rscript batch.R fb-test-1e7 $SGE_TASK_ID

