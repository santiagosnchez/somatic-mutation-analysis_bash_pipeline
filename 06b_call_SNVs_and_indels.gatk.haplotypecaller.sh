#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=30g,mem=30g,walltime=15:00:00
#PBS -e ${tumor}__${normal}.haplotypecaller.${index}.log
#PBS -j eo
# scheduler settings

# set date to calculate running time
start=$(date)

# load modules
module load java/1.8
#module load gatk/4.2.2.0
module load samtools/1.10

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# log date
echo $start

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
# run gatk's haplotypecaller
if [[ ! -e haplotypecaller/${tumor}__${normal}.haplotypecaller.unfiltered.${mode}.merged.vcf ]]; then
$gatk_path/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./.tmp" HaplotypeCaller \
 -I ${dir}/${tumor}.bqsr.bam \
 -I ${dir}/${normal}.bqsr.bam \
 -R ${reference} \
 -O haplotypecaller/${tumor}__${normal}.haplotypecaller.unfiltered.${mode}.${index}.vcf \
 --max-mnp-distance 0 \
 -L ${bed30intervals}/${bed}
 #  -bamout mutect2/${tumor}__${normal}.${index}.bam \
 #  --create-output-bam-index \
else
  # for $? = 0
  ls &> /dev/null
fi

# check if finished
check_finish=$?

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    # check if all HaplotypeCaller operations finished
    # first check for files
    ls all_logfiles/${tumor}__${normal}.haplotypecaller.[1-9]*.log &> /dev/null
    # if HC is still running
    if [[ "$?" == 0 ]]; then
        hc_logfiles=$(ls all_logfiles/${tumor}__${normal}.haplotypecaller.[1-9]*.log | wc -l)
        # try to wrap up in one go
        if [[ "${hc_logfiles}" == 29 ]]; then
            # gather vcffiles
            # generate list of files with their own -I flag
            vcffiles=$(ls haplotypecaller/${tumor}__${normal}.haplotypecaller.unfiltered.${mode}.*.vcf  | sort -V | sed 's/^/-I /')
            $gatk_path/gatk GatherVcfs $vcffiles -O haplotypecaller/${tumor}__${normal}.haplotypecaller.unfiltered.${mode}.merged.vcf
            # delete if finished
            if [[ "$?" == 0 ]]; then
                rm haplotypecaller/${tumor}__${normal}.haplotypecaller.unfiltered.${mode}.[1-9]*.vcf
                rm haplotypecaller/${tumor}__${normal}.haplotypecaller.unfiltered.${mode}.[1-9]*.vcf.idx
            fi
            # log to main
            echo "06: all scattered HC calls merged for ${tumor}__${normal}." | tee -a main.log
            # submit bcftools filtering
            qsub -v \
tumor=${tumor},\
normal=${normal},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/06c_call_SNVs_and_indels.bcftools.filter.sh
            # first scatter
            first_scatter_date=$(ls ${tumor}__${normal}.haplotypecaller.${index}.log ${tumor}__${normal}.haplotypecaller.[0-9]*.log | \
                   parallel 'head -2 {} | tail -1' | parallel date --date={} +%s | sort -n | parallel date --date=@{} | head -1)
            runtime=$( how_long "${first_scatter_date}" h )
            # log
            echo "06: ${tumor}__${normal} HaplotypeCaller took ${runtime} hours" | tee -a main.log
            # concat logfiles
            cat $(ls ${tumor}__${normal}.haplotypecaller.${index}.log all_logfiles/${tumor}__${normal}.haplotypecaller.[0-9]*.log | sort -V) > all_logfiles/${tumor}__${normal}.haplotypecaller.log
            rm $(ls all_logfiles/${tumor}__${normal}.haplotypecaller.[0-9]*.log )
            rm ${tumor}__${normal}.haplotypecaller.${index}.log
        else
            # log to main
            echo "06: ${tumor}__${normal} HaplotypeCaller variant calling completed for interval ${index}." | tee -a main.log
            # move logfile
            mv ${tumor}__${normal}.haplotypecaller.${index}.log all_logfiles
        fi
    # no scattered logfiles found
    else
        # check if HC is running
        ls ${tumor}__${normal}.haplotypecaller.[1-9]*.log &> /dev/null
        if [[ "$?" == 0 ]]; then
            # log to main
            echo "06: ${tumor}__${normal} HaplotypeCaller variant calling completed for interval ${index}." | tee -a main.log
            # move logfile
            mv ${tumor}__${normal}.haplotypecaller.${index}.log all_logfiles
        fi
    fi
fi
