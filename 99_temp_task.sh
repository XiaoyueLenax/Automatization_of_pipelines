#!/usr/bin/env bash
#SBATCH --cpus-per-task=20
#SBATCH --mem=90GB
#SBATCH --time=1-00:00:00
#SBATCH --job-name=kraken2_db_installer
#SBATCH --mail-user=xiaoyue.deng@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/storage/scratch/users/xd22m086/90_output_error/output/output_kraken2_db_installer_%j.o
#SBATCH --error=/storage/scratch/users/xd22m086/90_output_error/error/error_kraken2_db_installer_%j.e
#SBATCH --partition=epyc2

conda activate metawrap
module load vital-it/7
module load UHTS/Quality_control/fastqc/0.11.9
module load UHTS/Analysis/trimmomatic/0.36
module load Blast/blast/2.2.26
module load Blast/ncbi-blast/2.9.0+
module load UHTS/Aligner/bowtie2/2.3.4.1
module load UHTS/Assembler/SPAdes/3.15.4
module load UHTS/Assembler/megahit/1.1.4
module load UHTS/Quality_control/quast/4.6.0
module load UHTS/Analysis/samtools/1.10
module load UHTS/Analysis/metabat/2.12.1
module load Anaconda3

WORKDIR=/storage/scratch/users/xd22m086/04_metawrap_testground
RAW_DATA_DIR=/storage/scratch/users/xd22m086/02_human_reads
OUT_DIR=$WORKDIR/OUTPUT
QC_DIR=$OUT_DIR/QC
ASSEMBLY_DIR=$OUT_DIR/ASSEMBLY
BINNING_DIR=$OUT_DIR/BINNING
BB_DIR=$OUT_DIR/BLOBOLOGY
KRAKEN2_DIR=$OUT_DIR/KRAKEN2
Blast_db=/storage/scratch/users/rj23k073/programs/BLAST/Database    #Here we uses Russ' database to perform blast later
echo "Input directory = $RAW_DATA_DIR, output = $OUT_DIR"
IN_DIR=/storage/scratch/users/xd22m086/04_metawrap_testground/OUTPUT/ASSEMBLY/megahit
sample_name=("Human1" "Human2" "Human3")
human_num="Human1"
metawrap binning -o ${BINNING_DIR}/concoct/${human_num} -t 50 -a ${ASSEMBLY_DIR}/MEGAHIT/${human_num}/final_assembly.fasta --concoct ${RAW_DATA_DIR}/${human_num}/${human_num}_R_1.fastq ${RAW_DATA_DIR}/${human_num}/${human_num}_R_2.fastq
        