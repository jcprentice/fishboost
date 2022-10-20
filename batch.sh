#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N fb-mm1
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


# Run the program
# Rscript batch.R sim $SGE_TASK_ID # 1-220
# Rscript batch.R sim-events1 $SGE_TASK_ID # 1-30
# Rscript batch.R sim-events2 $SGE_TASK_ID # 1-60
# Rscript batch.R sim-fixed_seed $SGE_TASK_ID # 1-180
# Rscript batch.R sim-test_cov $SGE_TASK_ID # 1-120
# Rscript batch.R sim-simple1 $SGE_TASK_ID # 1-30
# Rscript batch.R sim-simple2 $SGE_TASK_ID # 1-60
# Rscript batch.R sim-donor_links1 $SGE_TASK_ID # 1-30
# Rscript batch.R sim-donor_links2 $SGE_TASK_ID # 1-60
# Rscript batch.R sim-donor_links3 $SGE_TASK_ID # 1-40
# Rscript batch.R sim-donor_links4 $SGE_TASK_ID # 1-80
# Rscript batch.R sim-Gsir_cov0_Dlidr-1 $SGE_TASK_ID # 1-20
# Rscript batch.R sim-Gsir_cov0_Dlidr-2 $SGE_TASK_ID # 1-40
# Rscript batch.R sim-inf $SGE_TASK_ID # 1-140

Rscript batch.R fb-minmax1 $SGE_TASK_ID # 1-60
# Rscript batch.R fb-minmax2 $SGE_TASK_ID # 1-120
