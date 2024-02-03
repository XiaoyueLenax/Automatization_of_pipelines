#!/usr/bin/env bash
#SBATCH --cpus-per-task=18
#SBATCH --mem=80GB
#SBATCH --time=7-00:00:00
#SBATCH --job-name=mamba_installer
#SBATCH --mail-user=xiaoyue.deng@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/storage/scratch/users/xd22m086/90_output_error/output/output_mamba_installer_%j.o
#SBATCH --error=/storage/scratch/users/xd22m086/90_output_error/error/error_mamba_installer_%j.e
#SBATCH --partition=epyc2


# for squeezemeta
mamba create -n SqueezeMeta -c conda-forge -c bioconda -c anaconda -c fpusan  squeezemeta=1.6 --no-channel-priority

#for metawrap
mamba install --only-deps -c ursky metawrap-mg