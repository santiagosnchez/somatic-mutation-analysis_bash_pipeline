#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=10g,mem=10g,walltime=10:00:00
#PBS -e ${tumor}__${normal}.GetPileupSummaries.log
#PBS -j eo
# scheduler settings

# load modules
module load java/1.8
#module load gatk/4.0.1.2
module load samtools/1.10

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# create dir for contamination
if [[ ! -e contamination ]]; then
    mkdir contamination
fi

# load reference path and other reference files
# for details check script
source /hpf/largeprojects/tabori/santiago/pipeline/export_paths_to_reference_files.sh
# change intervals to null if not WES
if [[ "${mode}" != "wes" ]]; then
    intervals=${gnomad_resource}
fi

# run gatk's GetPileupSummaries
gatk-4.2.2.0 GetPileupSummaries \
-I BQSR/${tumor}.bqsr.bam \
-V ${gnomad_resource} \
-L ${intervals} \
-O contamination/${tumor}.getpileupsummaries.table

gatk-4.2.2.0 GetPileupSummaries \
-I BQSR/${normal}.bqsr.bam \
-V ${gnomad_resource} \
-L ${intervals} \
-O contamination/${normal}.getpileupsummaries.table

# check if finished
check_finish=$?

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    mv ${tumor}__${normal}.GetPileupSummaries.log all_logfiles
    # next round of jobs are submitted manually or not
    qsub -v normal=${normal},tumor=${tumor},mode=${mode} ${pipeline_dir}/06b_check_crosscontamination.gatk.CalculateContamination.sh
fi


