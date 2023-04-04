#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N fb-prst
#$ -cwd
#$ -o ostr
#$ -e estr
#$ -l h_rt=72:00:00     # <-- ** make sure this is correct!!! **
#$ -l h_vmem=2G
#  These options are:
#  job name: -N
#  use the current working directory: -cwd
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem

# Initialise the environment modules
. /etc/profile.d/modules.sh

# Load R
module load roslin/R/4.2.0
module load openmpi/4.1.1

# Run the program
# Rscript batch.R sim $SGE_TASK_ID # 11
# Rscript batch.R sim-events1-mpi $SGE_TASK_ID # 3
# Rscript batch.R sim-donor_links3-mpi $SGE_TASK_ID # 4
# Rscript batch.R sim-donor_links4-mpi $SGE_TASK_ID # 4
# Rscript batch.R sim-Gsir_cov0_Dlidr-1-mpi $SGE_TASK_ID # 2
# Rscript batch.R sim-censored-2 $SGE_TASK_ID # 3
# Rscript batch.R sim-patch-fb-2 $SGE_TASK_ID # 3



# Rscript batch.R fb $SGE_TASK_ID # 9
# Rscript batch.R fb-minmax2-mpi $SGE_TASK_ID # 12
# Rscript batch.R fb-parasites $SGE_TASK_ID # 12
Rscript batch.R fb-parasites4 $SGE_TASK_ID # 6
# Rscript batch.R fb-ge $SGE_TASK_ID # 4

