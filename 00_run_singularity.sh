#!/usr/bin/env bash
#SBATCH --cpus-per-task=6
#SBATCH --mem=20GB
#SBATCH --time=3:00:00
#SBATCH --job-name=interative_singularity

# Adjust parameters above depending on your requirement 

srun --time=03:00:00 --mem-per-cpu=2G bash singularity shell my_image.sif