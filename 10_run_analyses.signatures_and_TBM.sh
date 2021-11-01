#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=10g,mem=10g,walltime=5:00:00
#PBS -e ${tumor}__${normal}.analyses.log
#PBS -j eo
# scheduler settings

# load modules
module load java/1.8
module load gatk/4.0.1.2
module load samtools/1.10
module load bcftools/1.11

# set working dir
cd $PBS_O_WORKDIR

# create output dirs
if [[ ! -e analyses ]]; then
    mkdir analyses
fi

# load reference path and other reference files
# for details check script
source /hpf/largeprojects/tabori/santiago/pipeline/export_paths_to_reference_files.sh
# change intervals to null if not WES
if [[ "${mode}" != "wes" ]]; then
    intervals=null
fi

# estimate coverage
coverage=$(samtools depth -b $intervals_bed -q20 -Q20 -d1000 BQSR/${tumor}.bqsr.bam BQSR/${normal}.bqsr.bam | awk '$3 >= 4 && $4 >= 4' | wc -l)
# expected coverage
expected=$(cat $intervals_bed | awk '{ count = count + ($3 - ($2 + 1)) } END { print count }')

# estimate tumor mutation burden (TMB)
# use prev coverage estimate

# total snvs
total_snvs=$(bcftools view --types snps vcf/${tumor}__${normal}.mutect2.annotated.${mode}.vcf.gz | grep -v "^#" | wc -l)
# total indels
total_indels=$(bcftools view --types indels vcf/${tumor}__${normal}.mutect2.annotated.${mode}.vcf.gz | grep -v "^#" | wc -l)
# calc TMB
TMB_snvs=$( echo "scale=4; ${total_snvs}/(${coverage}/1000000)" | bc )
TMB_indels=$( echo "scale=4; ${total_indels}/(${coverage}/1000000)" | bc )

# output
echo "${tumor},${normal},${coverage},${expected},${total_snvs},${total_indels},${TMB_snvs},${TMB_indels}" >> analyses/coverage_and_tmb.csv

# extract COSMIC signatures


# add more steps

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


