#!/bin/bash 

#---------------------------------------------------------------------------------
# Account information

#SBATCH --account=faculty              # basic (default), staff, phd, faculty

#---------------------------------------------------------------------------------
# Resources requested

#SBATCH --partition=standard       # standard (default), long, gpu, mpi, highmem
#SBATCH --cpus-per-task=1          # number of CPUs requested (for parallel tasks)
#SBATCH --mem-per-cpu=15G           # requested memory
#SBATCH --time=0-01:00:00          # wall clock limit (d-hh:mm:ss)

#---------------------------------------------------------------------------------
# Job specific name (helps organize and track progress of jobs)

#SBATCH --job-name=simTransition_noEps_noshocks   # user-defined job name

#---------------------------------------------------------------------------------
# Print some useful variables

echo "Job ID: $SLURM_JOB_ID"
echo "Job User: $SLURM_JOB_USER"
echo "Num Cores: $SLURM_JOB_CPUS_PER_NODE"

#---------------------------------------------------------------------------------
# Load necessary modules for the job

module load matlab/R2021b

#---------------------------------------------------------------------------------
# Commands to execute below...

srun matlab -nodisplay < simTransition_noEps_noshocks.m