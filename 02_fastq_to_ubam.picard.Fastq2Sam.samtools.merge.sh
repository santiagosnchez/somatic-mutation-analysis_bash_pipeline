#!/bin/bash
#PBS -l nodes=1:ppn=8,vmem=100g,mem=100g,walltime=10:00:00
#PBS -e ${sample}.FastqToSam.merge.log
#PBS -j eo
# scheduler settings

# load modules
module load java/1.8
module load picard-tools/2.18.0
#module load samtools/1.10
module load parallel/20210322

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# picard
# don't change! or uncomment
export picard_jar_file=/hpf/tools/centos6/picard-tools/2.18.0/picard.jar

# load all paths
source /hpf/largeprojects/tabori/santiago/pipeline/export_paths_to_reference_files.sh

# create dir for ubam
if [[ ! -e unmapped_bam ]]; then
    mkdir unmapped_bam
fi

# run picard in parallel
cat file_list.csv | grep "^${sample}," | parallel --colsep="," java -jar $picard_jar_file FastqToSam F1={4} F2={5} O=unmapped_bam/{1}_{2}.unaligned.bam SM={1} RG={2} LB={3} PL=Illumina MAX_RECORDS_IN_RAM=100000

# check that command finished
check_finish1=$?

# if only one lane change name
if [[ "$check_finish1" == 0 ]]; then 
    # make input files into single string-var
    input=$(for i in unmapped_bam/${sample}_*unaligned.bam; do echo -n "I=$i "; done)
    # run picard to merge bam files
    eval "java -jar $picard_jar_file MergeSamFiles SO=unsorted $input O=unmapped_bam/${sample}.unaligned.merged.bam  MAX_RECORDS_IN_RAM=100000"
    check_finish2=$?
    if [[ "$check_finish2" == 0 ]]; then
        rm unmapped_bam/${sample}_*.unaligned.bam
    fi
fi

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# submit next job if previous commands succeeded
if [[ "$check_finish1" == 0 && "$check_finish2" == 0 ]]; then
    # and move logfiles to folder
    mv ${sample}.FastqToSam.merge.log all_logfiles
    # calculate the amount of time needed to run ngm
    size=$(du -s unmapped_bam/${sample}.unaligned.merged.bam | cut -f1)
    walltime=$(echo "((($size / 150000)*10)/60)" | bc)
    # submit job
    qsub -l walltime="${walltime}:00:00" -v sample=${sample},wt=${walltime},mode=${mode} ${pipeline_dir}/03_align_bam_to_ref.ngm-core.sh
    echo "submited next job"
else
    echo "check log. Step failed to complete"
fi

