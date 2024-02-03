#!/usr/bin/env bash
#SBATCH --cpus-per-task=24
#SBATCH --mem=GB
#SBATCH --time=10:00:00
#SBATCH --job-name=KRAKEN
#SBATCH --mail-user=xiaoyue.deng@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/storage/scratch/users/xd22m086/90_output_error/output_installKrakendb_%j.o
#SBATCH --error=/storage/scratch/users/xd22m086/90_output_error/error_installkrakendb_%j.e
#SBATCH --partition=epyc2

kraken-build --standard --threads 24 --db MY_KRAKEN_DATABASE
kraken-build --db MY_KRAKEN_DATABASE --clean