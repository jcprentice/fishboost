#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N sim-ev1mp
#$ -cwd
#$ -o ostr
#$ -e estr
#$ -l h_rt=47:00:00     # <-- ** make sure this is correct!!! **
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
# Rscript batch.R sim $SGE_TASK_ID # 1-220
Rscript batch.R sim-events1-mpi $SGE_TASK_ID # 3 / 10
# Rscript batch.R sim-events2-mpi $SGE_TASK_ID # 60 / 4
# Rscript batch.R sim-donor_links3-mpi $SGE_TASK_ID # 4 / 10
# Rscript batch.R sim-donor_links4-mpi $SGE_TASK_ID # 80 / 4
# Rscript batch.R sim-Gsir_cov0_Dlidr-1-mpi $SGE_TASK_ID # 2 / 10
# Rscript batch.R sim-Gsir_cov0_Dlidr-2-mpi $SGE_TASK_ID # 40 / 4
# Rscript batch.R fb-minmax-mpi $SGE_TASK_ID # 40 / 4


# Rscript batch.R fb-minmax1-mpi $SGE_TASK_ID # 1-6 / 10
# Rscript batch.R fb-minmax2-mpi $SGE_TASK_ID # 1-120 / 4
