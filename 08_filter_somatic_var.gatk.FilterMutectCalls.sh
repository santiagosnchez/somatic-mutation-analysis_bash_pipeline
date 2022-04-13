#!/bin/bash
#PBS -l mem=10g,walltime=5:00:00
#PBS -e ${tumor}__${normal}.FilterMutectCalls.log
#PBS -j eo
# scheduler settings

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

# load reference path and other reference files
# for details check script
source ${pipeline_dir}/export_paths_to_reference_files.sh ${organism} ${genome} ${mode}

if [[ -e contamination/${tumor}__${normal}.calculatecontamination.table && -e contamination/${tumor}__${normal}.tumorsegmentation.table ]]; then
# run gatk's FilterMutectCalls
$gatk_path/gatk --java-options "-Djava.io.tmpdir=./.tmp" FilterMutectCalls \
 -R $reference \
 -V mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf \
 --contamination-table contamination/${tumor}__${normal}.calculatecontamination.table \
 --tumor-segmentation contamination/${tumor}__${normal}.tumorsegmentation.table \
 --ob-priors mutect2/f1r2/${tumor}__${normal}.read-orientation-model.tar.gz \
 -O mutect2/${tumor}__${normal}.mutect2.filtered.${mode}.vcf
 # skipping read orientation filtering due to high numbers of false negatives
 $gatk_path/gatk --java-options "-Djava.io.tmpdir=./.tmp" FilterMutectCalls \
  -R $reference \
  -V mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf \
  --contamination-table contamination/${tumor}__${normal}.calculatecontamination.table \
  --tumor-segmentation contamination/${tumor}__${normal}.tumorsegmentation.table \
  -O mutect2/${tumor}__${normal}.mutect2.filtered_no-obpriors.${mode}.vcf

# select passed variants
# $gatk_path/gatk --java-options "-Djava.io.tmpdir=./.tmp" SelectVariants \
#  -V mutect2/${tumor}__${normal}.mutect2.filtered.${mode}.vcf \
#  --exclude-filtered \
#  -O mutect2/${tumor}__${normal}.mutect2.selected.${mode}.vcf

# compress and index
index-vcf mutect2/${tumor}__${normal}.mutect2.filtered.${mode}.vcf
index-vcf mutect2/${tumor}__${normal}.mutect2.filtered_no-obpriors.${mode}.vcf

# normalize vcf file, compress, and tabix
bcftools norm -m- -f ${reference} mutect2/${tumor}__${normal}.mutect2.filtered.${mode}.vcf.gz > mutect2/${tumor}__${normal}.mutect2.filtered-norm.${mode}.vcf
index-vcf mutect2/${tumor}__${normal}.mutect2.filtered-norm.${mode}.vcf

# selection done later

else
    # resubmit until file is available with dependency
    # first get the job id if the GetPileupSummaries command
    if [[ ${normal} == "PON" ]]; then
        # filter without calc contamination
        # run gatk's FilterMutectCalls
        $gatk_path/gatk --java-options "-Djava.io.tmpdir=./.tmp" FilterMutectCalls \
         -R $reference \
         -V mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf \
         --ob-priors mutect2/f1r2/${tumor}__${normal}.read-orientation-model.tar.gz \
         -O mutect2/${tumor}__${normal}.mutect2.filtered.${mode}.vcf
         # skipping read orientation filtering due to high numbers of false negatives
         $gatk_path/gatk --java-options "-Djava.io.tmpdir=./.tmp" FilterMutectCalls \
          -R $reference \
          -V mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf \
          -O mutect2/${tumor}__${normal}.mutect2.filtered_no-obpriors.${mode}.vcf

        # select passed variants
        # $gatk_path/gatk --java-options "-Djava.io.tmpdir=./.tmp" SelectVariants \
        #  -V mutect2/${tumor}__${normal}.mutect2.filtered.${mode}.vcf \
        #  --exclude-filtered \
        #  -O mutect2/${tumor}__${normal}.mutect2.selected.${mode}.vcf

        # compress and index
        index-vcf mutect2/${tumor}__${normal}.mutect2.filtered.${mode}.vcf
        index-vcf mutect2/${tumor}__${normal}.mutect2.filtered_no-obpriors.${mode}.vcf

        # normalize vcf file, compress, and tabix
        bcftools norm -m- -f ${reference} mutect2/${tumor}__${normal}.mutect2.filtered.${mode}.vcf.gz > mutect2/${tumor}__${normal}.mutect2.filtered-norm.${mode}.vcf
        index-vcf mutect2/${tumor}__${normal}.mutect2.filtered-norm.${mode}.vcf

        # selection done later

    else
        if [[ -e ${tumor}__${normal}.CalculateContamination.log ]]; then
            # get job ids
            jobid_cc=$(head -1 ${tumor}__${normal}.CalculateContamination.log | sed 's/Job Id: //' )
            # resubmit as dependency job
            # log
            echo "08: Waiting for Calculate Contamination to finish for ${tumor}__${normal}: ${jobid_cc}" | tee -a main.log
            qsub -W depend=afterok:${jobid_cc} -v \
tumor=${tumor},\
normal=${normal},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/08_filter_somatic_var.gatk.FilterMutectCalls.sh
        # wait for file
        else
            echo "08: Waiting for Calculate Contamination to start for ${tumor} and ${normal}" | tee -a main.log
            qsub -v \
file="all_logfiles/${tumor}__${normal}.CalculateContamination.log",\
sample=${sample},\
script=08_filter_somatic_var.gatk.FilterMutectCalls.sh,\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/wait_for_file.sh
        fi
    fi
fi

# check if finished
check_finish=$?

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    # log to main
    echo "08: FilterMutectCalls completed for ${tumor}__${normal}. Submitting VCF annotation SnpEff and Funcotator" | tee -a main.log
    # next round of jobs are submitted manually or not
    # annotate VCF file
    qsub -v \
tumor=${tumor},\
normal=${normal},\
tissue="Somatic",\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/09_variant_annotation.snpEff-funcotator.sh
    # log to main
    echo "08: FilterMutectCalls completed for ${tumor}__${normal}. Submitting VCF annotation Annovar" | tee -a main.log
    qsub -v \
tumor=${tumor},\
normal=${normal},\
tissue="Somatic",\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/09a_variant_annotation.annovar.sh
    # move log
    mv ${tumor}__${normal}.FilterMutectCalls.log all_logfiles
fi
