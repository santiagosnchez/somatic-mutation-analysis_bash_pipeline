#!/bin/bash
#PBS -l mem=10g,walltime=5:00:00
#PBS -e ${tumor}__${normal}.FilterMutectCalls.log
#PBS -j eo
# scheduler settings

# load modules
module load java/1.8
module load gatk/4.2.2.0
module load samtools/1.10
module load bcftools/1.11
module load tabix

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# load reference path and other reference files
# for details check script
source /hpf/largeprojects/tabori/santiago/pipeline/export_paths_to_reference_files.sh
# change intervals to null if not WES
if [[ "${mode}" != "wes" ]]; then
    intervals=null
fi

if [[ -e contamination/${tumor}__${normal}.calculatecontamination.table && -e contamination/${tumor}__${normal}.tumorsegmentation.table ]]; then
# run gatk's FilterMutectCalls
gatk FilterMutectCalls \
 -R $reference \
 -V mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf \
 --contamination-table contamination/${tumor}__${normal}.calculatecontamination.table \
 --ob-priors mutect2/f1r2/${tumor}__${normal}.read-orientation-model.tar.gz \
 --tumor-segmentation contamination/${tumor}__${normal}.tumorsegmentation.table \
 -O mutect2/${tumor}__${normal}.mutect2.filtered.${mode}.vcf 

# select passed variants
gatk SelectVariants \
 -V mutect2/${tumor}__${normal}.mutect2.filtered.${mode}.vcf \
 --exclude-filtered \
 -O mutect2/${tumor}__${normal}.mutect2.selected.${mode}.vcf 

# compress and tabix selected file
bgzip -f mutect2/${tumor}__${normal}.mutect2.selected.${mode}.vcf
tabix mutect2/${tumor}__${normal}.mutect2.selected.${mode}.vcf.gz

# normalize vcf file, compress, and tabix
bcftools norm -m- -f ${reference} -Oz mutect2/${tumor}__${normal}.mutect2.selected.${mode}.vcf.gz > mutect2/${tumor}__${normal}.mutect2.normalized.${mode}.vcf.gz
tabix mutect2/${tumor}__${normal}.mutect2.normalized.${mode}.vcf.gz

else
    # resubmit until file is available with dependency
    # first get the job id if the GetPileupSummaries command
    if [[ -e ${tumor}__${normal}.GetPileupSummaries.log ]]; then
        getpileupsum_job=$(head -1 ${tumor}__${normal}.GetPileupSummaries.log)
    elif [[ -e all_logfiles/${tumor}__${normal}.GetPileupSummaries.log ]]; then
        getpileupsum_job=$(head -1 all_logfiles/${tumor}__${normal}.GetPileupSummaries.log)
    fi
    # then search for the second job id for CalculateContamination
    running_jobid=$(qstat -f -u `whoami` ${getpileupsum_job} | grep "beforeok" | sed 's/.*://')
    #running_jobid=$( head -1 ${tumor}__${normal}.CalculateContamination.log  )
    echo "waiting for CalculateContamination to finish"
    qsub -W depend=afterok:${running_jobid} -v tumor=${tumor},normal=${normal},mode=${mode} ${pipeline_dir}/08_filter_somatic_var.gatk.FilterMutectCalls.sh
    exit 0
fi

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
    qsub -v tumor=${tumor},normal=${normal},mode=${mode} ${pipeline_dir}/09_variant_annotation.snpEff.sh
fi


