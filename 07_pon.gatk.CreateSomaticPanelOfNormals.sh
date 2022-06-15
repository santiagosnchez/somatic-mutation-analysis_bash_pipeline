#!/bin/bash
#PBS -l mem=10g,walltime=5:00:00
#PBS -e gatk.CreateSomaticPanelOfNormals.log
#PBS -j eo
# scheduler settings

# set date to calculate running time
start=$(date)

# load modules
module load java/1.8
#module load gatk/4.2.2.0
module load samtools/1.10
module load bcftools/1.11
module load tabix

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# create tmp dir
if [[ ! -e .tmp ]]; then
    mkdir .tmp
fi

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# create PoN dir
if [[ ! -e PoN ]]; then
    mkdir PoN
fi

# load reference path and other reference files
# for details check script
source ${pipeline_dir}/00_export_pipeline_environment.sh ${organism} ${genome} ${mode}

# fetch sample names
samples=$(cat file_list.csv | cut -d, -f1 | sort -u)
vcffiles=""
# check if there are VCFs for all samples
for sample in ${samples}; do
  if [[ ! -e mutect2/${sample}.mutect2.pon.merged.vcf.gz ]]; then
    echo "07_pon: File mutect2/${sample}.mutect2.pon.merged.vcf.gz not found" | tee -a main.log
  else
    vcffiles="${vcffiles}-V mutect2/${sample}.mutect2.pon.merged.vcf.gz "
  fi
done
# count total samples
total_samples=$(echo $vcffiles | tr ' ' '\n' | grep -c '^\-V$')
echo "07_pon: Building PoN with ${total_samples} samples." | tee -a main.log

# fetch date for timestamp
date_for_pon=$(date -I)

# run gatk's GenomicsDBImport
$gatk_path/gatk --java-options "-Djava.io.tmpdir=./.tmp" GenomicsDBImport \
 -R $reference \
 -L $intervals \
 --genomicsdb-workspace-path PoN/pon_db \
 $vcffiles

# run gatk's CreateSomaticPanelOfNormals
$gatk_path/gatk --java-options "-Djava.io.tmpdir=./.tmp" CreateSomaticPanelOfNormals \
 -R reference.fasta \
 -V gendb://PoN/pon_db \
 -O PoN/pon.${total_samples}_samples.${date_for_pon}.vcf.gz

# check if finished
check_finish=$?

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    # log to main
    echo "07_pon: CreateSomaticPanelOfNormals completed." | tee -a main.log
    # calc runtime
    runtime=$( how_long "${start}" h )
    echo "07_pon: Step gatk.CreateSomaticPanelOfNormals.log took ${runtime} hours" | tee -a main.log
    # echo delete dirs and clean-up
    if [[ -e aligned_bam ]]; then
      rm -rf aligned_bam
    elif [[ -e preprocessed_bam ]]; then
      rm -rf preprocessed_bam
    elif [[ -e BQSR ]]; then
      mv BQSR bam
    fi
    rm -rf .tmp
    # change permissions
    find all_logfiles -type f -name "*.log" -exec chmod 644 {} \;
    # move log
    mv gatk.CreateSomaticPanelOfNormals.log all_logfiles

fi
