#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=30g,mem=30g,walltime=10:00:00
#PBS -e ${tumor}__${normal}.CalculateContamination.log
#PBS -j eo
# scheduler settings

# load modules
module load java/1.8
#module load gatk/4.2.2.0
module load samtools/1.10

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# create dir for contamination
if [[ ! -e contamination ]]; then
    mkdir contamination
fi
# set bam dir
if [[ ! -e bam ]]; then
    dir=BQSR
else
    dir=bam
fi

# load reference path and other reference files
# for details check script
source /hpf/largeprojects/tabori/shared/software/somatic-mutation-discovery/export_paths_to_reference_files.sh
# change intervals to null if not WES
if [[ "${mode}" != "wes" ]]; then
    intervals=null
fi

# run gatk's CalculateContamination
$gatk_path/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./tmp" CalculateContamination \
 -I contamination/${tumor}.getpileupsummaries.table \
 -matched contamination/${normal}.getpileupsummaries.table \
 -O contamination/${tumor}__${normal}.calculatecontamination.table \
 --tumor-segmentation contamination/${tumor}__${normal}.tumorsegmentation.table

# check if finished
check_finish=$?

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    # move logfile
    mv ${tumor}__${normal}.CalculateContamination.log all_logfiles
    # log to main
    echo "${tumor}__${normal} Calculate Contamination completed." | tee -a main.log
    # next round of jobs are submitted manually or not
fi
