#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=10g,mem=10g,walltime=5:00:00
#PBS -e ${tumor}__${normal}.bcftools.filter.log
#PBS -j eo
# scheduler settings

# set date to calculate running time
start=$(date)

# load modules
module load bcftools/1.11

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# create output dirs
if [[ ! -e haplotypecaller ]]; then
    mkdir -p haplotypecaller
fi

# create tmp dir
if [[ ! -e .tmp ]]; then
    mkdir .tmp
fi

# set bam dir
if [[ ! -e bam ]]; then
    dir=BQSR
else
    dir=bam
fi

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# load reference path and other reference files
# for details check script
source ${pipeline_dir}/00_export_pipeline_environment.sh ${organism} ${genome} ${mode}

#if [[ "${gvcf}" == 0 ]]; then

#else
# filter VCF with bcftools

if [[ -e haplotypecaller/${tumor}__${normal}.haplotypecaller.unfiltered.${mode}.merged.vcf ]]; then

# first normalize VCF to split multiple alleles into different lines
# and left align indels
# then filter by genotype quality (GQ), and discordant tumor-normal calls
# and annotate calls
bcftools norm -m- -f ${reference} haplotypecaller/${tumor}__${normal}.haplotypecaller.unfiltered.${mode}.merged.vcf | \
 bcftools filter -e 'GQ[0] < 20 || GQ[1] < 20' -s LowQual | \
 bcftools filter -e '(GT[0] = "AA" && GT[1] != "AA") | (GT[0] != "AA" && GT[1] = "AA")' -s Diff | \
 bcftools filter -e '(GT[0] = "RA" && GT[1] != "RA") | (GT[0] != "RA" && GT[1] = "RA")' -s Diff > \
 haplotypecaller/${tumor}__${normal}.haplotypecaller.filtered-norm.${mode}.vcf

# keep variants that PASS in a different VCF
bcftools view -f PASS haplotypecaller/${tumor}__${normal}.haplotypecaller.filtered-norm.${mode}.vcf > haplotypecaller/${tumor}__${normal}.haplotypecaller.selected.${mode}.vcf

  # index VCFs
  index-vcf haplotypecaller/${tumor}__${normal}.haplotypecaller.filtered-norm.${mode}.vcf
  index-vcf haplotypecaller/${tumor}__${normal}.haplotypecaller.selected.${mode}.vcf

  # check if finished
  check_finish=$?

else
   echo "06: No input VCF named haplotypecaller/${tumor}__${normal}.haplotypecaller.unfiltered.${mode}.merged.vcf" | tee -a main.log
   check_finish=1
fi


# check if command finished
if [[ "$check_finish" == 0 ]]; then
  # log to main
  echo "06: Germline (haplotypecaller) variant calling finished for ${tumor}__${normal}." | tee -a main.log
  # delete pileups
  #rm varscan/pileups/${tumor}.pileup
  # if [[ -e varscan/pileups/${normal}.pileup ]]; then
  #     # how many T are paired with N
  #     paired_tumors=$(grep -c ",${normal}$" tumors_and_normals.csv)
  #     if [[ $(ls varscan/*__${normal}.varscan.all.Somatic.hc.${mode}.vcf.gz | wc -l) == "${paired_tumors}" ]]; then
  #         rm varscan/pileups/${normal}.pileup
  #     fi
  # fi
  # submit annotation
  qsub -v \
tumor=${tumor},\
normal=${normal},\
tissue="Germline",\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/09b_variant_annotation.snpEff-funcotator.sh
  # calc runtime
  runtime=$( how_long "${start}" h )
  echo "06: Step ${tumor}__${normal}.bcftools.filter.log took ${runtime} hours" | tee -a main.log
  # move logfile
  mv ${tumor}__${normal}.bcftools.filter.log all_logfiles
fi
