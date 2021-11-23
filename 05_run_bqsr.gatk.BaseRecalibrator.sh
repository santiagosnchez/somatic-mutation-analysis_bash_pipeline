#!/bin/bash
#PBS -l nodes=1:ppn=12,vmem=30g,mem=30g
#PBS -e ${sample}.BQSR.log
#PBS -j eo
# scheduler settings

# load modules
module load java/1.8
module load gatk/4.2.2.0
module load samtools/1.10
module load parallel/20210322

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# create dir for BQSR
# and check if bam dir exists
# main dir with be changed to dir variable
if [[ ! -e bam ]]; then
    if [[ ! -e BQSR ]]; then
        mkdir BQSR
    fi
    dir=BQSR
else
    dir=bam
fi

# load reference path and other reference files
# for details check script
source /hpf/largeprojects/tabori/santiago/pipeline/export_paths_to_reference_files.sh
# change intervals to null if not WES
if [[ "${mode}" != "wes" ]]; then
    intervals=null
fi

# check files of previous run
if [[ -e ${dir}/${sample}.bqsr.bam ]]; then
    if [[ $(samtools quickcheck ${dir}/${sample}.bqsr.bam && echo 1) == 1 ]]; then
        if [[ ! -e ${dir}/${sample}.bqsr.bai ]]; then
            java -jar $picard_jar_file BuildBamIndex I=${dir}/${sample}.bqsr.bam
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
        rm ${dir}/${sample}.bqsr.*
        echo "resubmitting... increasing time by 2 hrs"
        wt=$(( wt + 2 ))
        qsub -l walltime=${wt}:00:00 -v sample=${sample},wt=${wt},mode=${mode} 05_run_bqsr.gatk.BaseRecalibrator.sh
        exit 0
    fi
else

# run gatl's BaseRecalibrator
gatk --java-options "-Xmx10G -XX:+UseParallelGC -XX:ParallelGCThreads=12 -Djava.io.tmpdir=./tmp" BaseRecalibrator \
 -R $reference \
 -L $intervals \
 --known-sites $knownsites \
 -I preprocessed_bam/${sample}.markdup.bam \
 -O ${sample}.baserecalibrator.txt

# use samtools to increase the compression level in the read-corrected output bam
gatk --java-options "-Xmx10G -XX:+UseParallelGC -XX:ParallelGCThreads=12 -Dsamjdk.compression_level=6 -Djava.io.tmpdir=./tmp" ApplyBQSR \
 --bqsr ${sample}.baserecalibrator.txt \
 -I preprocessed_bam/${sample}.markdup.bam \
 -O ${dir}/${sample}.bqsr.bam \
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
    #if [[ -e BQSR/${sample}.bqsr.bam ]]; then
        #if [[ $(samtools quickcheck BQSR/${sample}.bqsr.bam && echo 1) == 1 ]]; then
        #    rm preprocessed_bam/${sample}.*bam*
        #    rm aligned_bam/${sample}.*bam*
        #    rm unmapped_bam/${sample}.*bam*
        #fi
    #fi
    # next round of jobs are submitted manually or not
    # log to main
    echo "BQSR has been completed for sample ${sample}." | tee -a main.log
    # check if file exists and continue
    if [[ -e tumors_and_normals.csv ]]; then
        cat tumors_and_normals.csv | grep "^${sample},"
        if [[ "$?" == 0 ]]; then
            for line in `cat tumors_and_normals.csv | grep "^${sample},"`; do
                # first element is tumor, second is normal
                tumor=$(echo $line | sed 's/,.*//')
                normal=$(echo $line | sed 's/^.*,//')
                # do all file existance and integrity checks
                check_tumor=$(samtools quickcheck ${dir}/${tumor}.bqsr.bam && echo 1)
                check_normal=$(samtools quickcheck ${dir}/${normal}.bqsr.bam && echo 1)
                if [[ "${check_normal}" == 1 && "${check_tumor}" == 1 ]]; then
                    # submit all mutect2 jobs
                    export normal
                    export tumor
                    # start crosscontamination analyses
                    first_jobid=$(qsub -v normal=${normal},tumor=${tumor},mode=${mode} ${pipeline_dir}/06a_check_crosscontamination.gatk.GetPileupSummaries.sh)
                    # submit second crosscheck as dependency
                    qsub -W depend=afterok:${first_jobid} -v normal=${normal},tumor=${tumor},mode=${mode} ${pipeline_dir}/06b_check_crosscontamination.gatk.CalculateContamination.sh
                    # save a dry run of commands
                    ls $bed30intervals | grep ".bed" | parallel --tmpdir ./tmp --dry-run "qsub -v normal=${normal},tumor=${tumor},bed={},mode=${mode},index={#} ${pipeline_dir}/06c_call_SNVs_and_indels.gatk.mutect2.sh" > all_logfiles/${tumor}__${normal}.mutect2.0.log
                    # submit mutect2 jobs on 30 intervals
                    ls $bed30intervals | grep ".bed" | parallel --tmpdir ./tmp "qsub -v normal=${normal},tumor=${tumor},bed={},mode=${mode},index={#} ${pipeline_dir}/06c_call_SNVs_and_indels.gatk.mutect2.sh" | tee -a main.log
                    # move logfiles
                    if [[ "$?" == 0 ]]; then
                        mv ${tumor}.BQSR.log ${tumor}.baserecalibrator.txt all_logfiles
                        # log to main
                        echo "Mutect2 has started for ${tumor}__${normal} successfully." | tee -a main.log
                    fi
                elif [[ "${check_normal}" != 1 && "${check_tumor}" == 1 ]]; then
                    # resubmit with dependency
                    if [[ -e ${normal}.BQSR.log ]]; then
                        # wait until BQSR finishes
                        echo "${tumor} (tumor) waiting for BQSR ${normal} (normal) to finish: ${running_jobid}" | tee -a main.log
                        # get jobid from first line of log
                        running_jobid=$( head -1 ${normal}.BQSR.log )
                        qsub -W depend=afterok:${running_jobid} -v sample=${tumor},mode=${mode} ${pipeline_dir}/05_run_bqsr.gatk.BaseRecalibrator.sh
                        exit 0
                    elif [[ -e all_logfiles/${normal}.BQSR.log ]]; then
                        qsub -v sample=${tumor},mode=${mode} ${pipeline_dir}/05_run_bqsr.gatk.BaseRecalibrator.sh
                    else
                        # wait for the BQSR script to start
                        echo "${tumor} (tumor) waiting for ${normal} (normal) BQSR to start." | tee -a main.log
                        qsub -v file="${normal}.BQSR.log",sample=${sample},mode=${mode},script=05_run_bqsr.gatk.BaseRecalibrator.sh ${pipeline_dir}/wait_for_file.sh
                        exit 0
                    fi
                fi
            done
        else
            echo "sample $sample not in tumor column (1st)" | tee -a main.log
            mv ${sample}.BQSR.log ${sample}.baserecalibrator.txt all_logfiles
            exit 0
        fi
    else
        echo -e "could not find \"tumors_and_normals.csv\" file" | tee -a main.log
        exit 1
    fi
fi
