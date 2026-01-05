#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#  These options are:
#  job name: -N
#  use the current working directory: -cwd
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem
#$ -N flatten
#$ -pe sharedmem 1
#$ -t 1-25
#$ -R y
#$ -P roslin_wilson
#$ -cwd
#$ -o out
#$ -e out
#$ -l h_rt=2:00:00     # <-- ** make sure this is correct!!! **

# Initialise the environment modules
. /etc/profile.d/modules.sh

module load R/4.5

# Run the program
Rscript batch_flatten_bici_states.R fb-test $SGE_TASK_ID # 2

