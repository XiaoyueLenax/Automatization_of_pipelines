#!/usr/bin/env bash
#SBATCH --cpus-per-task=6
#SBATCH --mem=64GB
#SBATCH --time=20:00:00
#SBATCH --job-name=multi_qc_hu1
#SBATCH --mail-user=xiaoyue.deng@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/storage/scratch/users/xd22m086/90_output_error/output_multiqc_%j.o
#SBATCH --error=/storage/scratch/users/xd22m086/90_output_error/error_multiqc_%j.e
#SBATCH --partition=epyc2


module load vital-it/7 
module load UHTS/Analysis/MultiQC/1.8
ex=/storage/scratch/users/xd22m086/02_human_reads/Human1/Human1_R_1.fastq
ex=/storage/scratch/users/xd22m086/02_human_reads/Human1/Human1_R_2.fastq

ex2=/storage/scratch/users/xd22m086/02_human_reads/Human2/Human2_R_1.fastq
ex2=/storage/scratch/users/xd22m086/02_human_reads/Human2/Human2_R_2.fastq

#define path
input=/storage/scratch/users/xd22m086/02_human_reads/QC
output=/storage/scratch/users/xd22m086/02_human_reads/QC/multiqc

fastqc
multiqc $input -o $output