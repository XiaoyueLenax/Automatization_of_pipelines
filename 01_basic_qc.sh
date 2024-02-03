#!/usr/bin/env bash
#SBATCH --cpus-per-task=6
#SBATCH --mem=64GB
#SBATCH --time=20:00:00
#SBATCH --job-name=basic_qc_hu1
#SBATCH --mail-user=xiaoyue.deng@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/storage/scratch/users/xd22m086/90_output_error/output_qc_%j.o
#SBATCH --error=/storage/scratch/users/xd22m086/90_output_error/error_qc_%j.e
#SBATCH --partition=epyc2

#call module
module load vital-it/7 
module load UHTS/Quality_control/fastqc/0.11.9

#define path
input=/storage/scratch/users/xd22m086/02_human_reads
output=/storage/scratch/users/xd22m086/02_human_reads/QC

#run qc - debugging of failure to create html, creating temporary file
fastqc -t 2 -o $output $input/Human1_*.fastq.gz
