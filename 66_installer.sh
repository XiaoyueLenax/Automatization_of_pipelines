#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=12GB
#SBATCH --time=3:00:00
#SBATCH --job-name=meta
#SBATCH --mail-user=xiaoyue.deng@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/storage/scratch/users/xd22m086/90_output_error/output/output_installer_%j.o
#SBATCH --error=/storage/scratch/users/xd22m086/90_output_error/error/error_installer_%j.e
#SBATCH --partition=epyc2


mamba create -n SqueezeMeta -c conda-forge -c bioconda -c anaconda -c fpusan  squeezemeta --no-channel-priority