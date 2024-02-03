#!/usr/bin/env bash
#SBATCH --cpus-per-task=60
#SBATCH --mem=80GB
#SBATCH --time=3-0:00:00
#SBATCH --job-name=Human3_kraken2_test
#SBATCH --mail-user=xiaoyue.deng@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/storage/scratch/users/xd22m086/90_output_error/output/output_Human3_kraken2_test_%j.o
#SBATCH --error=/storage/scratch/users/xd22m086/90_output_error/error/error_Human3_kraken2_test_%j.e
#SBATCH --partition=epyc2
#SBATCH --array=0-6


# ===========================================================================

#           Module setup - Loeading all required packages

#============================================================================
#module load Anaconda3
#Load all vital-IT packages preinstalled
conda activate metawrap
module load vital-it/7
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

# ====================================================================================

#                       User Input: Sets working directory

# ====================================================================================
WORKDIR=""   # Master directory where you want everything to be in there
RAW_DATA_DIR=""   # Raw data directory where the raw fastq files are 
        #WARNING: make sure your data is in _R1.fastq format!

# Replace these with your sample names so it can be analysed automatically in a loop.        
sample_names=("Human1" "Human2" "Human3") # Recommend to manually type in all your sample names to loop through later..
REF_GENOME="/storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/UHGG_reps.fasta" 


# ====================================================================================

#          Setting up pipeline structure + install metawrap packages

# ====================================================================================
# Data Strutcure -  makes directories if they do not already exist.
cd $WORKDIR
    
    mkdir -p Metawrap_Pipeline
    mkdir -p Scripts
    mkdir -p db
    mkdir -p OUTPUT

        # Here, where this script and the other scripts should be
        scp /storage/scratch/users/xd22m086/04_metawrap_testground/1_scripts_meta/meta_all.sh .
        scp /storage/scratch/users/xd22m086/04_metawrap_testground/1_scripts_meta/ncbi.sh .
        

    cd -p Metawrap_Pipeline
        mkdir -p db ; mkdir -p OUTPUT; mkdir -p metaWRAP


    cd metaWRAP

        # Install the metawrap pipeline
        git clone https://github.com/bxlab/metaWRAP.git
        echo 'export PATH="${WORKDIR}/metaWRAP/bin:$PATH"' >> ~/.bashrc


        # Here, echo and check whether metawrap is correctly configured and you can see it in the PATH.
        echo $PATH


            # In case you need to set up the environment on your own -- check the metawrap documentation
            # In case the documentation is too confusing, use this:
            #installations
            your_env_name="Insert_here"
            
            conda install -y mamba
            mamba create -y -n ${your_env_name} python=2.7
            conda activate your_env_name

            # configure channels: may not be necessary
            conda config --add channels defaults
            conda config --add channels conda-forge
            conda config --add channels bioconda
            conda config --add channels ursky

            # Install metawrap dependencies - may not be necessary
            mamba install --only-deps -c ursky metawrap-mg

            # Install metawrap
            mamba install --metawrap

            # Important step here after installation: Manually locate WIP (00 file)
    cd ../
    cd db
        # Caution: If you already have these databases installed, soft link it for faster performance, 
        #          since these databases are very large.
        mkdir -p kraken2
            cd kraken2
            # WIP - It is undecided whether to install the kraken or not, since with other packages storage is already full.
            # kraken-build --standard --threads 24 --db MY_KRAKEN_DATABASE
            # kraken-build --db MY_KRAKEN_DATABASE --clean
            #           Update your db link to the metawrap config file.
            # KRAKEN_DB=/path/to/my/database/MY_KRAKEN_DATABASE
            cd ../
        mkdir -p NCBI_nt
            cd NCBI_nt
            
            # Softlink the database from the repository.   
            ln -s /storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/NCBI_nt/* .
            cd ../

        mkdir -p NCBI_tax
            cd NCBI_tax

            # softlink from repository
            ln -s /storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/NCBI_tax/* .
            cd ../

        mkdir -p BMTAGGER_INDEX
        mkdir -p PROGRAMS
            cd PROGRAMS

            # Smaller program scripts are copied dicrectly for convinence. 
            scp /storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/NCBI_nt/ncbi-blast-2.14.0+/bin/blastn .
            blastn_script='/storage/scratch/users/xd22m086/04_metawrap_testground/DATABASE/NCBI_nt/ncbi-blast-2.14.0+/bin/blastn'
            cd ../
    cd OUTPUT
        mkdir -p QC
            cd QC
                mkdir -p QUAST
            cd ../
        
        mkdir -p ASSEMBLY
            cd ASSEMBLY
                mkdir -p MetaSPAdes
                mkdir -p MEGAHIT
            cd ../
       
        mkdir -p BINNING
        mkdir -p BLOBOLOGY
        mkdir -p KRAKEN2
            cd BINNING
                mkdir -p Metabat2
                mkdir -p Maxbin2
                mkdir -p CONCOCT
            cd ../

        mkdir -p ANNOTATION
            mkdir -p blastn
            mkdir -p megablast
            cd ../
        
# Sub folders set up complete
echo "Finished Setting up folder structures."

# Output directories-------------------------------------------------------------------
OUT_DIR=$WORKDIR/OUTPUT
QC_DIR=$OUT_DIR/QC
ASSEMBLY_DIR=$OUT_DIR/ASSEMBLY
BINNING_DIR=$OUT_DIR/BINNING
BB_DIR=$OUT_DIR/BLOBOLOGY
KRAKEN2_DIR=$OUT_DIR/KRAKEN2
Blast_db=/storage/scratch/users/rj23k073/programs/BLAST/Database
BLAST_DIR=$OUT_DIR/ANNOTATION    #Here we uses Russ' database to perform blast later
echo "Input directory = $RAW_DATA_DIR, output = $OUT_DIR"


# -------------------------------------------------------------------------------------------

#                 If you would like to use everything with metawrap:
#                                   Here is the guide

#            Modules:
#                    read_qc         Raw read QC module (read trimming and contamination removal)
#                    assembly        Assembly module (metagenomic assembly)
#                    kraken          KRAKEN module (taxonomy annotation of reads and assemblies)
#                    kraken2         KRAKEN2 module (taxonomy annotation of reads and assemblies)
#                    blobology       Blobology module (GC vs Abund plots of contigs and bins)

#                    binning         Binning module (metabat, maxbin, or concoct)
#                    bin_refinement  Refinement of bins from binning module
#                    reassemble_bins Reassemble bins using metagenomic reads
#                    quant_bins      Quantify the abundance of each bin across samples
#                    classify_bins   Assign taxonomy to genomic bins
#                    annotate_bins   Functional annotation of draft genomes

#                    --help | -h             show this help message
#                    --version | -v  show metaWRAP version
#                    --show-config   show where the metawrap configuration files are stored



#                  The section below runs the modules in one go.
#                   Make sure everything is configured correctly!
#


#       Step 1: Fast QC + multiQC

#           Step 2: Assembly

#               Step 2.5: MetaQUAST

#                   Step 3: Binning

#                        Step 3.5: Checkm

#                           Step 3.6: Quant_bins

#                                Step 3.7: Bin_refinement

#                                       Step 4: MEGAblast

#                                              Step 4.5 : Blastn

# Step 1 - FastQC--------------------------------------------------------------------------------------------------------------

#       The raw reads must be named with suffix *_R_1.fastq and *_R_2.fastq !!
for sample_id in "${sample_names[@]}"; do
    metawrap read_qc -1 $RAW_DATA_DIR/${sample_id}_R_1.fastq -2 -1 $RAW_DATA_DIR/${sample_id}_R_2.fastq -t 24 -o ${OUT_DIR}
    #Options:

    #    -1 STR          forward fastq reads
    #    -2 STR          reverse fastq reads
    #    -o STR          output directory
    #    -t INT          number of threads (default=1)
    #    -x STR          prefix of host index in bmtagger database folder (default=hg38)

    #    --skip-bmtagger         dont remove human sequences with bmtagger
    #    --skip-trimming         dont trim sequences with trimgalore
    #    --skip-pre-qc-report    dont make FastQC report of input sequences
    #    --skip-post-qc-report   dont make FastQC report of final sequences

done
        # Feedback Check:
        if [ $? -eq 0 ]; then
        echo "Step 1: QC completed successfully, output files at $OUT_DIR/QC"
        else
        echo "Step 1: QC failed. We recommend you restart again."
        exit 1
        fi


# Step 2 - Assembly: MetaSPAdes& megahit -------------------------------------------------------------------

    # Warning: This is the most computationally intense stage.
    # Make sure to allocate enough storage!
    for sample_id in "${sample_names[@]}"; do

        metawrap assembly -1 ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_1.fastq -2 ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_2.fastq -o ${ASSEMBLY_DIR}/megahit2/${sample_id} -t 60 -m 80
        metawrap assembly --metaspades -1 ${RAW_DATA_DIR}/${sample_id}_R_1.fastq -2 ${RAW_DATA_DIR}//${sample_id}_R_2.fastq -o ${ASSEMBLY_DIR}/metaspades2/${sample_id} -t 60 -m 60 

        # Usage: metaWRAP assembly [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir
        # Options:

        # -1 STR          forward fastq reads
        # -2 STR          reverse fastq reads
        # -o STR          output directory
        # -m INT          memory in GB (default=24)
        # -t INT          number of threads (defualt=1)
        # -l INT          minimum length of assembled contigs (default=1000)

        # --megahit       assemble with megahit (default)
        # --metaspades    assemble with metaspades instead of megahit (better results but slower and higher memory requirement)    
    
    done
    
    # Wait for all parallel tasks to finish
    wait
    
        if [ $? -eq 0 ]; then
        echo "Step 2: assembly completed. Check output at $OUT_DIR/ASSEMBLY"
        else
        echo "Step 2: assembly has failed, please check the input files again or try changing allocated resources."
        exit 1
        fi

# Step 2.5  MetaQUAST -----------------------------------------------------------------------------------------------------------
    
    # Run MetaQUAST to assess quality of the assemblies, independent from metaWRAP
    # With reference
        for sample_id in "${sample_names[@]}"; do
            metaquast.py ${ASSEMBLY_DIR}megahit/$sample_id/final_assembly.fasta --output-dir ${QC_DIR}/QUAST/refs/$sample_id -R ${REF_GENOME} -t 50
            # WIP: copy the metaquast.py to here
        done
        if [ $? -eq 0 ]; then
            echo "Step 2: metaQUAST completed. Check output at ${QC_DIR}/QUAST/refs/$sample_id"
        else
            echo "Step 2: failed, check error report"
            exit 1
        fi

    

# Step 2c   Kraken-------------------------------------------------------------------------------------------------------

#    This steps requires a functional database. 
#    WIP - inplement a mechanism that checks whether the kraken db is installed. 
#    for sample_id in "${sample_name[@]}"; do
#        metawrap kraken2 --no-preload -o KRAKEN -t $NUM_THREADS /storage/scratch/users/xd22m086/02_human_reads/Human3_R*.fastq;
#        ASSEMBLY/final_assembly.fasta -o ${OUT_DIR}/KRAKEN2/Human3
#    done


# Step 3 : Binning ------------------------------------------------------------------------------------------------------------
    
    # Will automatically perform binning on both megahit and metaspades file. 
    # If you only used 1 assembler, hash out the one you're not using.

    MEGAHIT_DIR=${ASSEMBLY_DIR}/MEGAHIT
    METASPADES_DIR=${ASSEMBLY_DIR}/MetaSPAdes
    CONCOCT_DIR=${OUT_DIR}/CONCOCT
    
    # Pipeline Option 1 - Russ's Concoct, independent of metawrap
    #  For MEGAHIT
        perl /storage/scratch/users/rj23k073/programs/Maxbin2/MaxBin-2.2.7/run_MaxBin.pl
        cut_up_fasta.py /storage/scratch/users/rj23k073/04_DEER/06_Assembly/6_2_deer.asm/scaffolds_filtered.fasta -c 10000 -o 0 --merge_last -b 6_2_deer_10K.bed > 6_2_deer_10K.fa
        concoct_coverage_table.py 6_2_deer_10K.bed /storage/scratch/users/rj23k073/04_DEER/07_BAM/6_2_filter.sorted.bam > 6_2_coverage_table.tsv
        concoct --composition_file 6_2_deer_10K.fa --coverage_file 6_2_coverage_table.tsv -b concoct_output/
        merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv
        mkdir concoct_output/fasta_bins
        extract_fasta_bins.py /storage/scratch/users/rj23k073/04_DEER/06_Assembly/6_2_deer.asm/scaffolds_filtered.fasta concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins
    # For MetaSPAdes
        perl /storage/scratch/users/rj23k073/programs/Maxbin2/MaxBin-2.2.7/run_MaxBin.pl
        cut_up_fasta.py /storage/scratch/users/rj23k073/04_DEER/06_Assembly/6_2_deer.asm/scaffolds_filtered.fasta -c 10000 -o 0 --merge_last -b 6_2_deer_10K.bed > 6_2_deer_10K.fa
        concoct_coverage_table.py 6_2_deer_10K.bed /storage/scratch/users/rj23k073/04_DEER/07_BAM/6_2_filter.sorted.bam > 6_2_coverage_table.tsv
        concoct --composition_file 6_2_deer_10K.fa --coverage_file 6_2_coverage_table.tsv -b concoct_output/
        merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv
        mkdir concoct_output/fasta_bins
        extract_fasta_bins.py /storage/scratch/users/rj23k073/04_DEER/06_Assembly/6_2_deer.asm/scaffolds_filtered.fasta concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins

    for sample_id in "${sample_name[@]}"; do

        # ALL FOR METAHIT RESULTS
        # You may not want metabat 1 now there's metabat2, but this is an option
        metawrap binning -o ${BINNING_DIR}/metabat2/${sample_id} -t 50 -a ${ASSEMBLY_DIR}/MetaSPAdes/${sample_id}/final_assembly.fasta --metabat2 ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_1.fastq ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_2.fastq
        # metawrap binning -o ${BINNING_DIR}/metabat1/${sample_id} -t 50 -a ${ASSEMBLY_DIR}/MetaSPAdes/${sample_id}final_assembly.fasta --metabat1 ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_1.fastq ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_2.fastq
        metawrap binning -o ${BINNING_DIR}/maxbin2/${sample_id} -t 50 -a ${ASSEMBLY_DIR}/MetaSPAdes/${sample_id}final_assembly.fasta --maxbin2 ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_1.fastq ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_2.fastq
        metawrap binning -o ${BINNING_DIR}/concoct/${sample_id} -t 50 -a ${ASSEMBLY_DIR}/MetaSPAdes/${sample_id}final_assembly.fasta --concoct ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_1.fastq ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_2.fastq
        
        # ALL FOR METASPADES RESULTS
        metawrap binning -o ${BINNING_DIR}/metabat2/${sample_id} -t 50 -a ${ASSEMBLY_DIR}/MEGAHIT/${sample_id}/final_assembly.fasta --metabat2 ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_1.fastq ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_2.fastq
        # metawrap binning -o ${BINNING_DIR}/metabat1/${sample_id} -t 50 -a ${ASSEMBLY_DIR}/MEGAHIT/${sample_id}/final_assembly.fasta --metabat1 ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_1.fastq ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_2.fastq
        metawrap binning -o ${BINNING_DIR}/maxbin2/${sample_id} -t 50 -a ${ASSEMBLY_DIR}/MEGAHIT/${sample_id}/final_assembly.fasta --maxbin2 ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_1.fastq ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_2.fastq
        metawrap binning -o ${BINNING_DIR}/concoct/${sample_id} -t 50 -a ${ASSEMBLY_DIR}/MEGAHIT/${sample_id}/final_assembly.fasta --concoct ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_1.fastq ${RAW_DATA_DIR}/${sample_id}/${sample_id}_R_2.fastq
        #  Parameters used:
        #       -a STR          metagenomic assembly file
        #        -o STR          output directory
        #        -t INT          number of threads (default=1)
        #        -m INT          amount of RAM available (default=4)
        #        -l INT          minimum contig length to bin (default=1000bp). Note: metaBAT will default to 1500bp minimum

        #        --metabat2      bin contigs with metaBAT2
        #        --metabat1      bin contigs with the original metaBAT
        #        --maxbin2       bin contigs with MaxBin2
        #        --concoct       bin contigs with CONCOCT

        #        --universal     use universal marker genes instead of bacterial markers in MaxBin2 (improves Archaea binning)
        #        --run-checkm    immediately run CheckM on the bin results (requires 40GB+ of memory)
        #        --single-end    non-paired reads mode (provide *.fastq files)
        #        --interleaved   the input read files contain interleaved paired-end reads
    done
        # Performance check
        if [ $? -eq 0 ]; then
        echo "Step 3: Binning completed successfully, output files at ${BINNING_DIR}"
        else
        echo "Step 3: Binning has failed !! Check the error report for more details."
        exit 1
        fi



# Step 3.5 checkm and GC plot ---------------------------------------------------------------------------------
#                      To add or remove options, use this guide below of checkm.

#                               ...::: CheckM v1.2.2 :::...

#           Lineage-specific marker set:
#               tree         -> Place bins in the reference genome tree
#               tree_qa      -> Assess phylogenetic markers found in each bin
#               lineage_set  -> Infer lineage-specific marker sets for each bin

#           Taxonomic-specific marker set:
#               taxon_list   -> List available taxonomic-specific marker sets
#               taxon_set    -> Generate taxonomic-specific marker set

#           Apply marker set to genome bins:
#               analyze      -> Identify marker genes in bins
#               qa           -> Assess bins for contamination and completeness

#           Common workflows (combines above commands):
#               lineage_wf   -> Runs tree, lineage_set, analyze, qa
#               taxonomy_wf  -> Runs taxon_set, analyze, qa

#           Reference distribution plots:
#               gc_plot      -> Create GC histogram and delta-GC plot
#               coding_plot  -> Create coding density (CD) histogram and delta-CD plot
#               tetra_plot   -> Create tetranucleotide distance (TD) histogram and delta-TD plot
#               dist_plot    -> Create image with GC, CD, and TD distribution plots together

#           General plots:
#               nx_plot      -> Create Nx-plots
#               len_hist     -> Sequence length histogram
#               marker_plot  -> Plot position of marker genes on sequences
#               gc_bias_plot -> Plot bin coverage as a function of GC

#           Bin exploration and modification:
#               unique       -> Ensure no sequences are assigned to multiple bins
#               merge        -> Identify bins with complementary sets of marker genes
#                outliers     -> [Experimental] Identify outlier in bins relative to reference distributions
#               modify       -> [Experimental] Modify sequences in a bin

#           Utility functions:
#               unbinned     -> Identify unbinned sequences
#               coverage     -> Calculate coverage of sequences
#               tetra        -> Calculate tetranucleotide signature of sequences
#               profile      -> Calculate percentage of reads mapped to each bin
#               ssu_finder   -> Identify SSU (16S/18S) rRNAs in sequences


for sample_id in "${samples[@]}"; do
	
    checkm lineage_wf -t 50 -x fa ${BIN_DIR}/${sample_id}/metabat2_bins ${OUT_DIR}/${sample_id}/checkm
    #runs tree, lineage_set, analyze, qa
    # usage: checkm lineage_wf [-h] [-r] [--ali] [--nt] [-g] [-u UNIQUE] [-m MULTI] [--force_domain] [--no_refinement] 
    #                     [--individual_markers] [--skip_adj_correction] [--skip_pseudogene_correction] [--aai_strain AAI_STRAIN] 
    #                     [-a ALIGNMENT_FILE] [--ignore_thresholds]
    #                     [-e E_VALUE] [-l LENGTH] [-f FILE] [--tab_table] [-x EXTENSION] [-t THREADS] [--pplacer_threads PPLACER_THREADS] [-q]
    #                     [--tmpdir TMPDIR]
    #                     bin_input output_dir

    checkm gc_plot --dpi 1000 -x fa ${BIN_DIR}/${sample_id}/metabat2_bins ${OUT_DIR}/${sample_id}/checkm 1
    # gc_plot: create GC histogram and delta-GC plot
    # usage: checkm gc_plot [-h] [--image_type {eps,pdf,png,ps,svg}] [--dpi DPI] [--font_size FONT_SIZE] 
    #                   [-x EXTENSION] [--width WIDTH] [--height HEIGHT]
    #                  [-w GC_WINDOW_SIZE] [-b GC_BIN_WIDTH] [-q]
    #                  bin_input output_dir dist_value [dist_value ...]
done

    # Performance check
        if [ $? -eq 0 ]; then
        echo "Step 3.5: Checkm completed successfully, output files at ${OUT_DIR}/${sample_id}/checkm"
        else
        echo "Step 3.5: Checkm has failed !! Check the error report for more details."
        exit 1
        fi

# Step 3.6 QUANT_bins----------------------------------------------------------------------------------------------------
#       This is a build in metawrap module.
#       This module takes in a set of bins and any number of paired-end read sets from metagenomic samples, and estimates the abundance of each bin in each 
#       sample with salmon. It then uses Seaborn to make clustered heatmaps genome abundances.
#
#       Author of pipeline: German Uritskiy. I do not clain any authorship of the many programs this pipeline uses.
#       For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.


# Set here the bins that you would like to use, by default is the metabat2

    SEL_BIN_DIR=${BINNING_DIR}/metabat2
    # WIP: give an option to skil salmon?
#       For megahit
for sample_id in "${sample_name[@]}"; do
    #metawrap quant_bins -b ${BIN_DIR}/${sample_id}/original_bins -o ${OUT_DIR}/${sample_id} -a ${BIN_DIR}/${sample_id}/binned_assembly/assembly.fa ${RAW_DATA_DIR}/${sample_id}/*R_1.fastq ${RAW_DATA_DIR}/${sample_id}/*R_2.fastq -t 50
    metawrap quant_bins -b ${BIN_DIR}/${sample_id}/original_bins \
    -o ${OUT_DIR}/${sample_id} \
    -a ${MEGAHIT_DIR}/${sample_id}/final_assembly.fasta ${RAW_DATA_DIR}/${sample_id}/*R_1.fastq ${RAW_DATA_DIR}/${sample_id}/*R_2.fastq -t 50
    #   Usage: metaWRAP quant_bins [options] -b bins_folder -o output_dir -a assembly.fa readsA_1.fastq readsA_2.fastq ... [readsX_1.fastq readsX_2.fastq]
    #   Options:
    #        -b STR          folder containing draft genomes (bins) in fasta format
    #        -o STR          output directory
    #        -a STR          fasta file with entire metagenomic assembly (strongly recommended!)
    #        -t INT          number of threads
done

#        For metaSPADES
for sample_id in "${sample_name[@]}"; do
    #metawrap quant_bins -b ${BIN_DIR}/${sample_id}/original_bins -o ${OUT_DIR}/${sample_id} -a ${BIN_DIR}/${sample_id}/binned_assembly/assembly.fa ${RAW_DATA_DIR}/${sample_id}/*R_1.fastq ${RAW_DATA_DIR}/${sample_id}/*R_2.fastq -t 50
    metawrap quant_bins -b ${BIN_DIR}/${sample_id}/original_bins \
    -o ${OUT_DIR}/${sample_id} \
    -a ${METASPADES_DIR}/${sample_id}/final_assembly.fasta ${RAW_DATA_DIR}/${sample_id}/*R_1.fastq ${RAW_DATA_DIR}/${sample_id}/*R_2.fastq -t 50
    #   Usage: metaWRAP quant_bins [options] -b bins_folder -o output_dir -a assembly.fa readsA_1.fastq readsA_2.fastq ... [readsX_1.fastq readsX_2.fastq]
    #   Options:
    #        -b STR          folder containing draft genomes (bins) in fasta format
    #        -o STR          output directory
    #        -a STR          fasta file with entire metagenomic assembly (strongly recommended!)
    #        -t INT          number of threads
done

    # Performance check
        if [ $? -eq 0 ]; then
        echo "Step 3.6: Quant_bins completed successfully, output files at ${OUT_DIR}/${sample_id}/checkm"
        else
        echo "Step 3.6: Quant_bins has failed !! Check the error report for more details."
        exit 1
        fi


# Step 3.7 Bin refinement -------------------------------------------------------
#
#   This is a metawrap built in module that run on the outputs of binning.sh pipeline 
#   to analyze the metagenomic bins and arrive at the best possible putative genomes.
#   There are several options to give additional binning results for comparison. 

    for sample_id in "${sample_name[@]}"; do    
        metawrap bin_refinement -o ${SEL_BIN_DIR}/${sample_id}/refined_bins -A ${SEL_BIN_DIR}/${sample_id}/metabat2_bins -1 ${RAW_DATA_DIR}/${sample_id}/*_1.fastq -2 ${RAW_DATA_DIR}/${sample_id}/*_2.fastq -t 50

        #  Usage: metaWRAP bin_refinement [options] -o output_dir -A bin_folderA [-B bin_folderB -C bin_folderC]
        #       Note: the contig names in different bin folders must be consistant (must come from the same assembly).

        #   Options:

        #    -o STR          output directory
        #    -t INT          number of threads (default=1)
        #    -m INT          memory available (default=40)
        #    -c INT          minimum % completion of bins [should be >50%] (default=70)
        #    -x INT          maximum % contamination of bins that is acceptable (default=10)

        #    -A STR          folder with metagenomic bins (files must have .fa or .fasta extension)
        #    -B STR          another folder with metagenomic bins
        #    -C STR          another folder with metagenomic bins

        #    --skip-refinement       dont use binning_refiner to come up with refined bins based on combinations of binner outputs
        #    --skip-checkm           dont run CheckM to assess bins
        #    --skip-consolidation    choose the best version of each bin from all bin refinement iteration
        #    --keep-ambiguous        for contigs that end up in more than one bin, keep them in all bins (default: keeps them only in the best bin)
        #    --remove-ambiguous      for contigs that end up in more than one bin, remove them in all bins (default: keeps them only in the best bin)
        #    --quick                 adds --reduced_tree option to checkm, reducing runtime, especially with low memory

    done

    # Performance check
        if [ $? -eq 0 ]; then
        echo "Step 3.7: bin_refinement module completed successfully, output files at ${SEL_BIN_DIR}/${sample_id}/refined_bins"
        else
        echo "Step 3.7: bin_refinement module has failed !! Check the error report for more details."
        exit 1
        fi


# Step 4: MEGA blast -----------------------------------------------------------------------------------------
#       This is a built-in metawrap module, classify_bins.  
#       Classifies all contigs in a set of bins by aligning them to the NCBI database with MEGABLAST, 
#       pruning the resulting hits, and assigning final taxonomy with taxator-k. 
#       The consensus taxonomy of each bin is called by contructing aweighted consensus tree, 
#       and traversing the tree by the most likely path.

    CLASSIFY_BINS_DIR=${SEL_BIN_DIR}/Classified_bins
    for sample_id in "${sample_name[@]}"; do
        metawrap classify_bins -b ${SEL_BIN_DIR} -o ${CLASSIFY_BINS_DIR} -t 50
        #   Usage: metaWRAP classify_bins [options] -b bin_folder -o output_dir
        #   Options:

        #    -b STR          folder with the bins to be classified (in fasta format)
        #    -o STR          output directory
        #    -t INT          number of threads

    done

        #   Performance check
            if [ $? -eq 0 ]; then
            echo "Step 4: classify_bins module completed successfully, output files at ${CLASSIFY_BINS_DIR}"
            else
            echo "Step 4: classify_bins module has failed !! Check the error report for more details."
            exit 1
            fi

# Step 4.5: Regular blast ------------------------------------------------------------------------------------
#       This is a metawrap independent section, only utilizes blastn, instead of megablast.
#       Choose either one of these to perform. 

        # Input a bin of your choice here.
        # I have not figured out a way to perform blastn on all bins without them being shut down
        # due to limited computation cost.
        Selected_best_bin="" 

#       I cannot make them run in a loop. A tiny bin already takes ages.
#       So, make sure to select one bins per sample, and we can run them for each sample
#       For example, you can choose the longest bin. 
#       This is not a great way to deal with them, but we need a method to bypass the limits...

        for sample_id in "${sample_name[@]}"; do
            ${blastn_script} -num_threads 50\
            -db /storage/scratch/users/rj23k073/programs/BLAST/Database/nt\
            -outfmt 6\
            -query ${Selected_best_bin} > ${BLAST_DIR}/${sample_id}/test_bin1_Human1_out.raw.tab
        done 

                    #   Performance check
                    if [ $? -eq 0 ]; then
                    echo "Step 4。5: blastn module completed successfully, output files at ${CLASSIFY_BINS_DIR}"
                    else
                    echo "Step 4.5: blastn module has failed !! Check the error report for more details."
                    exit 1
                    fi

# WIP - bloblogy ---------------------------------------------------------------------------------
#       This is a built in metawrap module, which still needs to be configured.
#       modified pipeline from the program 'BLOBOLOGY', which produces a GC vs Coverage plot of the contigs, helping visualize bins and phylogenetic 
#       composition.The original script has been modified to be better suited for running on clusters.
#
#       Author of original pipeline: Sujai Kumar (https://github.com/blaxterlab). Author of modifications: German Uritskiy. I do not take any credit for the original pipeline.
#       For questions, bugs, and suggestions, contact German Uritskiy at guritsk1@jhu.edu.

    for sample_id in "${sample_name[@]}"; do
        metawrap blobology -t 50 -a ${OUT_DIR}/${sample_id} -a ${BIN_DIR}/${sample_id}/binned_assembly/assembly.fa -o ${OUT_DIR}/${sample_id} ${RAW_DATA_DIR}/${sample_id}/*R_1.fastq ${RAW_DATA_DIR}/${sample_id}/*R_2.fastq
        # Options:

        # -a STR          assembly fasta file
        # -o STR          output directory
        # -t INT          number of threads

        # --subsamble     INT     Number of contigs to run blobology on. Subsampling is randomized. (default=ALL)
        # --bins          STR     Folder containing bins. Contig names must match those of the assembly file. (default=None)
    done