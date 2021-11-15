#!/bin/bash
#SBATCH --error=out_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --partition=defq 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100m
#SBATCH --array=1-85

export OMP_NUM_THREADS=1 
module load julia-uon/gcc6.3.0/1.1.0
julia -L src/NetworkResilience.jl scripts/ndriver.jl $SLURM_ARRAY_TASK_ID 