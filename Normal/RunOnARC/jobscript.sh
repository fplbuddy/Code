#!/bin/bash
#SBATCH --job-name=kf32_mu5e-4
#SBATCH --time=48:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mail-user=winchester@maths.ox.ac.uk
#SBATCH --mail-type=END

. enable_arcus-b_mpi.sh

###mpirun $MPI_HOSTS ./run.exe -in $SLURM_JOB_ID.in &> $SLURM_JOB_ID.out
mpirun $MPI_HOSTS ./run.exe
