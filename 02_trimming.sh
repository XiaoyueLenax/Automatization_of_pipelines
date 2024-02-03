#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --mem=64GB
#SBATCH --time=20:00:00
#SBATCH --job-name=trimmomatic
#SBATCH --mail-user=xiaoyue.deng@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/storage/scratch/users/xd22m086/90_output_error/output_trimmomatic_%j.o
#SBATCH --error=/storage/scratch/users/xd22m086/90_output_error/error_trimmomatic_%j.e
#SBATCH --partition=epyc2

# run trimmomatic
INPUT_DIR="/storage/scratch/users/xd22m086/02_human_reads"
OUT_DIR="/storage/scratch/users/xd22m086/02_human_reads/all_in_one"
trimmomatic PE -threads 12 $INPUT_DIR/*R1.fastq $INPUT_DIR/*R2.fastq $OUT_DIR/QC $OUT_DIR/QC <output_R2.trimmed.fastq> <output_R2.unpaired.fastq> ILLUMINACLIP:<adapter_file>:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36