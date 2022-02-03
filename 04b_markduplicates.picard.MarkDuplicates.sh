#!/bin/bash
#PBS -l nodes=1:ppn=12,vmem=30g,mem=30g
#PBS -e ${sample}.picard.MarkDuplicates.log
#PBS -j eo
# scheduler settings

# load modules
module load java/1.8
module load gatk/4.2.2.0
#module load picard-tools/2.18.0
module load samtools/1.10
module load parallel/20210322

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# load all paths
source /hpf/largeprojects/tabori/shared/software/somatic-mutation-discovery/export_paths_to_reference_files.sh

# create dir for preprocessed bam
if [[ ! -e preprocessed_bam ]]; then
    mkdir preprocessed_bam
fi

# create tmp dir
if [[ ! -e tmp ]]; then
    mkdir tmp
fi

# run picards FixMateInformation and MarkDuplicates
if [[ -e aligned_bam/${sample}.bam && $(samtools quickcheck aligned_bam/${sample}.bam && echo 1) == 1 ]]; then
    if [[ ! -e preprocessed_bam/${sample}.picard.sorted.mfixed.bam ]]; then
gatk FixMateInformation \
 -I aligned_bam/${sample}.bam \
 -O preprocessed_bam/${sample}.picard.sorted.mfixed.bam \
 -SO coordinate \
 -COMPRESSION_LEVEL 6 -MAX_RECORDS_IN_RAM 100000 -CREATE_INDEX true
    else
        if [[ $(samtools quickcheck preprocessed_bam/${sample}.picard.sorted.mfixed.bam && echo 1) != 1 ]]; then
gatk FixMateInformation \
 -I aligned_bam/${sample}.bam \
 -O preprocessed_bam/${sample}.picard.sorted.mfixed.bam \
 -SO coordinate \
 -COMPRESSION_LEVEL 6 -MAX_RECORDS_IN_RAM 100000 -CREATE_INDEX true
        fi
    fi
    if [[ ! -e preprocessed_bam/${sample}.picard.sorted.markdup.bam ]]; then
gatk MarkDuplicates \
 -I preprocessed_bam/${sample}.picard.sorted.mfixed.bam \
 -O preprocessed_bam/${sample}..markdup.bam \
 -ASO coordinate \
 -TAGGING_POLICY All \
 -OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
 -METRICS_FILE all_logfiles/${sample}.picard.MarkDuplicates.txt \
 -COMPRESSION_LEVEL 6 -MAX_RECORDS_IN_RAM 100000 -CREATE_INDEX true
    else
        if [[ $(samtools quickcheck preprocessed_bam/${sample}.picard.sorted.markdup.bam && echo 1) != 1 ]]; then
gatk MarkDuplicates \
 -I preprocessed_bam/${sample}.picard.sorted.mfixed.bam \
 -O preprocessed_bam/${sample}.markdup.bam \
 -ASO coordinate \
 -TAGGING_POLICY All \
 -OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
 -METRICS_FILE all_logfiles/${sample}.picard.MarkDuplicates.txt \
 -COMPRESSION_LEVEL 6 -MAX_RECORDS_IN_RAM 100000 -CREATE_INDEX true
        fi
    fi
else
    echo "resubmitting previous step and increase time by 2hrs"
    # add two more hours of walltime
    wt=$(( wt + 2 ))
    # resubmit previous script and exit
    qsub -l walltime=${wt}:00:00 -v sample=${sample},wt=${wt},mode=${mode} 03_align_bam_to_ref.ngm-core.sh
    exit 0
fi

check_finish=$?

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

if [[ "$check_finish" == 0 ]]; then
     # move log files to dir
     mv ${sample}.picard.MarkDuplicates.log all_logfiles
     # remove unnecessary files
     rm $(ls preprocessed_bam/${sample}.picard.sorted.* | grep -v "markdup")
     # submit next job
     qsub -l walltime=${wt}:00:00 -v sample=${sample},wt=${wt},mode=${mode} ${pipeline_dir}/05_run_bqsr.gatk.BaseRecalibrator.sh
fi
