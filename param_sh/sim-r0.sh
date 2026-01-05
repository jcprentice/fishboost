#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#  These options are:
#  job name: -N
#  use the current working directory: -cwd
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem
#$ -N r0
#$ -t 1-500
#$ -R y
#$ -P roslin_wilson
#$ -cwd
#$ -o out
#$ -e out
#$ -l h_rt=1:00:00     # <-- ** make sure this is correct!!! **

# Initialise the environment modules
. /etc/profile.d/modules.sh

module load R/4.5

# Run the program
Rscript batch_selection_on_R0.R 1 $SGE_TASK_ID
Rscript batch_selection_on_R0.R 2 $SGE_TASK_ID
# Rscript batch_selection_on_R0.R 3 $SGE_TASK_ID
# Rscript batch_selection_on_R0.R 4 $SGE_TASK_ID
# Rscript batch_selection_on_R0.R 5 $SGE_TASK_ID

