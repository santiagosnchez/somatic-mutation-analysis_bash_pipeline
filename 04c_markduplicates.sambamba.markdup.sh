#!/bin/bash
#PBS -l nodes=1:ppn=10,vmem=30g,mem=30g
#PBS -e ${sample}.sambamba.markdup.log
#PBS -j eo
# scheduler settings

############### INFO #################
#                                    
# Integrated Bash Pipeline for       
# Somatic Mutation Discovery         
#                                    
# Author: Santiago Sanchez-Ramirez 
# Year: 2021  
# Email: santiago.snchez@gmail.com   
#                                    
# More info on README.md             
#                                    
# Notes:                             
# (1) Uses sambamba markdup approach rather than
#     Picards MarkDuplicates or samtools markdup
#
#####################################

#############
# script 4c #
#############


# load modules
module load sambamba/0.7.0
module load samtools/1.10
module load gatk/4.2.2.0

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# check if previous run finish
# after main if block

# load all paths
source /hpf/largeprojects/tabori/santiago/pipeline/export_paths_to_reference_files.sh

# create dir for preprocessed bam files
if [[ ! -e preprocessed_bam ]]; then
    mkdir preprocessed_bam
fi

# run sambamba to mark duplicates
if [[ -e aligned_bam/${sample}.merged.bam && $(samtools quickcheck aligned_bam/${sample}.merged.bam && echo 1) == 1 ]]; then
    sambamba markdup \
     -t 10 \
     aligned_bam/${sample}.merged.bam \
     preprocessed_bam/${sample}.markdup.bam
     # prev step already generates index
     # index bam
     # gatk BuildBamIndex -I preprocessed_bam/${sample}.markdup.bam 
else
    echo "resubmitting previous step and increase time by 2hrs"
    # add two more hours of walltime
    wt=$(( wt + 2 ))
    # resubmit previous script and exit
    qsub -l walltime=${wt}:00:00 -v sample=${sample},wt=${wt},mode=${mode} ${pipeline_dir}/03_merge_bams.sambamba.sh
    exit 0
fi

check_finish=$?

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

if [[ "$check_finish" == 0 ]]; then
     # remove unnecessary files
     rm preprocessed_bam/${sample}.merged.*
     qsub -l walltime=${wt}:00:00 -v sample=${sample},wt=${wt},mode=${mode} ${pipeline_dir}/05_run_bqsr.gatk.BaseRecalibrator.sh
     # move log files to dir
     mv ${sample}.sambamba.markdup.log all_logfiles
fi

