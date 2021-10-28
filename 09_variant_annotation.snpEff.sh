#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=30g,mem=30g,walltime=12:00:00
#PBS -e ${tumor}__${normal}.FilterMutectCalls.log
#PBS -j eo
# scheduler settings

# load modules
module load java/1.8
module load samtools/1.10
module load snpEff/4.11
module load bcftools/1.11

# set working dir
cd $PBS_O_WORKDIR

# create output dirs
if [[ ! -e vcf ]]; then
    mkdir -p vcf
fi


# load reference path and other reference files
# for details check script
source /hpf/largeprojects/tabori/santiago/pipeline/export_paths_to_reference_files.sh
# change intervals to null if not WES
if [[ "${mode}" != "wes" ]]; then
    intervals=null
fi

# normalize variants, reduce complex alleles
bcftools annotate -- mutect2/${tumor}__${normal}.mutect2.normalized.${mode}.vcf.gz 

# run gatk's mutect2
gatk FilterMutectCalls \
 -V mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf \
 --contamination-table contamination/${tumor}__${normal}.calculatecontamination.table \
 --ob-priors orientation/read-orientation-model.tar.gz \
 -O mutect2/${tumor}__${normal}.mutect2.filtered.${mode}.vcf 

# check if finished
check_finish=$?

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    mv ${tumor}__${normal}.FilterMutectCalls.log all_logfiles
    # next round of jobs are submitted manually or not
    # annotate VCF file
    # 
fi


