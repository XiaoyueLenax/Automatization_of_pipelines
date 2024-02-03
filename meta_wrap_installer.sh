#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --mem=40GB
#SBATCH --time=20:00:00
#SBATCH --job-name=meta_tagger_human_index
#SBATCH --mail-user=xiaoyue.deng@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/storage/scratch/users/xd22m086/90_output_error/output/output_tagger_human_index_%j.o
#SBATCH --error=/storage/scratch/users/xd22m086/90_output_error/error/error_tagger_human_index_%j.e
#SBATCH --partition=epyc2

#conda install -y -c ursky metawrap-mg
#conda install -y blas=2.5=mkl
#mamba install -y metawrap

#wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz" -P /storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/NCBI_nt
# after above is done, run this code:
#for a in nt.*.tar.gz
#do tar xzf $a
#done


#wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P /storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/NCBI_tax
#cd /storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/NCBI_tax
#tar -xvf taxdump.tar.gz


cd /storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/BMTAGGER_INDEX
#wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/*fa.gz -P /storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/BMTAGGER_INDEX
#gunzip *fa.gz
#cat *fa > hg38.fa
#rm chr*.fa
#
bmtool -d hg38.fa -o hg38.bitmask
srprism mkindex -i hg38.fa -o hg38.srprism -M 100000