#!/bin/bash
#PBS -l nodes=1:ppn=12,vmem=30g,mem=30g
#PBS -e ${sample}.BQSR.log
#PBS -j eo
# scheduler settings

# load modules
module load java/1.8
#module load gatk/4.0.1.2
module load picard-tools/2.18.0
module load samtools/1.10
module load parallel/20210322

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# picard jar file
export picard_jar_file=/hpf/tools/centos6/picard-tools/2.18.0/picard.jar

# check if prev step finished correctly
#if [[ ! -e BQSR/${sample}.bqsr.bam ]]; then
#    echo "resubmitting previous step and increase time by 2hrs"
#    wt=$(( wt + 2 ))
#    # can switch this to picards MarkDuplicates method
#    qsub -l walltime=${wt}:00:00 -v sample=${sample},wt=${wt} 04a_markduplicates.samtools.markdup.sh 
#    exit 0
#fi

# create dir for BQSR
if [[ ! -e BQSR ]]; then
    mkdir BQSR
fi

# load reference path and other reference files
# for details check script
source /hpf/largeprojects/tabori/santiago/pipeline/export_paths_to_reference_files.sh
# change intervals to null if not WES
if [[ "${mode}" != "wes" ]]; then
    intervals=null
fi

# check files of previous run
if [[ -e BQSR/${sample}.bqsr.bam ]]; then
    if [[ $(samtools quickcheck BQSR/${sample}.bqsr.bam && echo 1) == 1 ]]; then
        if [[ ! -e BQSR/${sample}.bqsr.bai ]]; then
            java -jar $picard_jar_file BuildBamIndex I=BQSR/${sample}.bqsr.bam
            #if [[ "$?" == 0 ]]; then
            #    rm preprocessed_bam/${sample}.*
            #    rm aligned_bam/${sample}.*
            #    rm unmapped_bam/${sample}.*
            #fi
            mv ${sample}.BQSR.log all_logfiles
            # submit next job
            qsub -v sample=${sample},wt=${wt},mode=${mode} ${pipeline_dir}/05_run_bqsr.gatk.BaseRecalibrator.sh 
            exit 0
        else
            ls &> /dev/null
        fi
    else
        rm BQSR/${sample}.bqsr.bam*
        echo "resubmitting previous step and increase time by 2 hrs"
        wt=$(( wt + 2 ))
        qsub -l walltime=${wt}:00:00 -v sample=${sample},wt=${wt},mode=${mode} 04a_markduplicates.samtools.markdup.sh
        exit 0
    fi
else

# run gatl's BaseRecalibrator
gatk-4.2.2.0 --java-options "-Xmx10G -XX:+UseParallelGC -XX:ParallelGCThreads=12" BaseRecalibrator \
 -R $reference \
 -L $intervals \
 --known-sites $knownsites \
 -I preprocessed_bam/${sample}.samtools.sorted.markdup.bam \
 -O ${sample}.baserecalibrator.txt

# use samtools to increase the compression level in the read-corrected output bam
gatk-4.2.2.0 --java-options "-Xmx10G -XX:+UseParallelGC -XX:ParallelGCThreads=12 -Dsamjdk.compression_level=6" ApplyBQSR \
 --bqsr ${sample}.baserecalibrator.txt \
 -I preprocessed_bam/${sample}.samtools.sorted.markdup.bam \
 -O BQSR/${sample}.bqsr.bam \
 --create-output-bam-index

#/dev/stdout | samtools view -@ 12 -hb --write-index --no-PG -o BQSR/${sample}.final.bam##idx##BQSR/${sample}.final.bam.bai -

fi

# collect last command status
check_finish=$?

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    if [[ -e BQSR/${sample}.bqsr.bam ]]; then
        #if [[ $(samtools quickcheck BQSR/${sample}.bqsr.bam && echo 1) == 1 ]]; then
        #    rm preprocessed_bam/${sample}.*bam*
        #    rm aligned_bam/${sample}.*bam*
        #    rm unmapped_bam/${sample}.*bam*
        #fi
        mv ${sample}.BQSR.log ${sample}.baserecalibrator.txt all_logfiles
    fi
    # next round of jobs are submitted manually or not
    if [[ -e tumors_and_normals.csv ]]; then
        for line in `cat tumors_and_normals.csv | grep "${sample}"`; do 
            tumor=$(echo $line | sed 's/,.*//')
            normal=$(echo $line | sed 's/^.*,//')
            # do all file existance and integrity checks
            check_tumor=$(samtools quickcheck BQSR/${tumor}.bqsr.bam && echo 1)
            check_normal=$(samtools quickcheck BQSR/${normal}.bqsr.bam && echo 1)
            if [[ "${check_normal}" == 1 && "${check_tumor}" == 1 ]]; then
                if [[ ! -e all_logfiles/${tumor}__${normal}.mutect2.0.log ]]; then
                    # submit all mutect2 jobs
                    export normal
                    export tumor
                    # start crosscontamination analyses
                    qsub -v normal=${normal},tumor=${tumor},mode=${mode} ${pipeline_dir}/06a_check_crosscontamination.gatk.GetPileupSummaries.sh 
                    # save a dry run of commands
                    ls $bed30intervals | grep ".bed" | parallel --dry-run "qsub -v normal=${normal},tumor=${tumor},bed={},mode=${mode},index={#} ${pipeline_dir}/06c_call_SNVs_and_indels.gatk.mutect2.sh" > all_logfiles/${tumor}__${normal}.mutect2.0.log
                    # submit mutect2 jobs on 30 intervals
                    ls $bed30intervals | grep ".bed" | parallel "qsub -v normal=${normal},tumor=${tumor},bed={},mode=${mode},index={#} ${pipeline_dir}/06c_call_SNVs_and_indels.gatk.mutect2.sh"
                fi
            elif [[ "${check_normal}" == 1 && "${check_tumor}" != 1 ]]; then
                # resubmit with dependency
                running_jobid=$( head -1 ${tumor}.BQSR.log )
                qsub -W afterok:${running_jobid} -v sample=${normal},mode=${mode} ${pipeline_dir}/05_run_bqsr.gatk.BaseRecalibrator.sh
            elif [[ "${check_normal}" != 1 && "${check_tumor}" == 1 ]]; then
                # resubmit with dependency
                running_jobid=$( head -1 ${normal}.BQSR.log )
                qsub -W afterok:${running_jobid} -v sample=${tumor},mode=${mode} ${pipeline_dir}/05_run_bqsr.gatk.BaseRecalibrator.sh
            fi
        done
    else
        echo "could not find \"tumors_and_normals.csv\" file"
        exit 1
    fi
fi

