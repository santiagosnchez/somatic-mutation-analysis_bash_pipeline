#!/bin/bash
#PBS -l nodes=1:ppn=12,vmem=30g,mem=30g
#PBS -e ${sample}.samtools.markdup.log
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
# (1) Uses samtools markdup approach rather than
#     Picards MarkDuplicates
#
#####################################

#############
# script 4a #
#############


# load modules
module load samtools/1.10

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# load all paths
source /hpf/largeprojects/tabori/shared/software/somatic-mutation-discovery/export_paths_to_reference_files.sh

# create dir for preprocessed bam files
if [[ ! -e preprocessed_bam ]]; then
    mkdir preprocessed_bam
fi

# create tmp dir
if [[ ! -e tmp ]]; then
    mkdir tmp
fi

# run samtools to sort and mark duplicates
if [[ -e aligned_bam/${sample}.bam && $(samtools quickcheck aligned_bam/${sample}.bam && echo 1) == 1 ]]; then
    if [[ ! -e preprocessed_bam/${sample}.samtools.readname.bam ]]; then
        samtools sort -n -@ 12 -o preprocessed_bam/${sample}.samtools.readname.bam aligned_bam/${sample}.bam
    else
        if [[ $(samtools quickcheck preprocessed_bam/${sample}.samtools.readname.bam && echo 1) != 1 ]]; then
            samtools sort -n -@ 12 -o preprocessed_bam/${sample}.samtools.readname.bam aligned_bam/${sample}.bam
        fi
    fi
    if [[ ! -e preprocessed_bam/${sample}.samtools.fixmate.bam ]]; then
        samtools fixmate -@ 12 -m preprocessed_bam/${sample}.samtools.readname.bam preprocessed_bam/${sample}.samtools.fixmate.bam
    else
        if [[ $(samtools quickcheck preprocessed_bam/${sample}.samtools.fixmate.bam && echo 1) != 1 ]]; then
            samtools fixmate -@ 12 -m preprocessed_bam/${sample}.samtools.readname.bam preprocessed_bam/${sample}.samtools.fixmate.bam
        fi
    fi
    if [[ ! -e preprocessed_bam/${sample}.samtools.sorted.bam ]]; then
        samtools sort -@ 12 -o preprocessed_bam/${sample}.samtools.sorted.bam preprocessed_bam/${sample}.samtools.fixmate.bam
    else
        if [[ $(samtools quickcheck preprocessed_bam/${sample}.samtools.sorted.bam && echo 1) != 1 ]]; then
            samtools sort -@ 12 -o preprocessed_bam/${sample}.samtools.sorted.bam preprocessed_bam/${sample}.samtools.fixmate.bam
        fi
    fi
    if [[ ! -e preprocessed_bam/${sample}.samtools.markdup.bam ]]; then
        samtools markdup -@ 12 -d 2500 --reference ${reference} -s -f all_logfiles/${sample}.markdup.stats.log preprocessed_bam/${sample}.samtools.sorted.bam preprocessed_bam/${sample}.markdup.bam
        samtools index preprocessed_bam/${sample}.markdup.bam
    else
        if [[ $(samtools quickcheck preprocessed_bam/${sample}.samtools.markdup.bam && echo 1) != 1 ]]; then
            samtools markdup -@ 12 -d 2500 --reference ${reference} -s -f all_logfiles/${sample}.markdup.stats.log preprocessed_bam/${sample}.samtools.sorted.bam preprocessed_bam/${sample}.markdup.bam
            samtools index preprocessed_bam/${sample}.markdup.bam
        fi
    fi
else
    echo "resubmitting previous step and increase time by 2hrs"
    # add two more hours of walltime
    wt=$(( wt + 2 ))
    # resubmit previous script and exit
    qsub -l walltime=${wt}:00:00 -v sample=${sample},wt=${wt},mode=${mode} ${pipeline_dir}/03_align_bam_to_ref.ngm-core.sh
    exit 0
fi

check_finish=$?

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

if [[ "$check_finish" == 0 ]]; then
     # move log files to dir
     mv ${sample}.samtools.markdup.log all_logfiles
     # remove unnecessary files
     rm $(ls preprocessed_bam/${sample}.samtools.* | grep -v "markdup")
     qsub -l walltime=${wt}:00:00 -v sample=${sample},wt=${wt},mode=${mode} ${pipeline_dir}/05_run_bqsr.gatk.BaseRecalibrator.sh
     # log to main
     echo "04: duplicate reads have been marked for ${sample} and preceeding files have been deleted." | tee -a main.log
fi
