#!/bin/bash
#PBS -l nodes=1:ppn=12,vmem=30g,mem=30g
#PBS -e ${sample}.BQSR.log
#PBS -j eo
# scheduler settings

# load modules
module load java/1.8
#module load gatk/4.2.2.0
module load samtools/1.10
module load parallel/20210322

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# write job details to log
qstat -f $PBS_JOBID >> ${sample}.BQSR.log

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

# create .tmp dir
if [[ ! -e .tmp ]]; then
    mkdir .tmp
fi
# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# debug
echo "alignment-only mode: ${aln_only}"
# if aln_only variable is not set, do full analysis.
if [[ -z $aln_only ]]; then
    aln_only=0
fi

# load reference path and other reference files
# for details check script
source ${pipeline_dir}/export_paths_to_reference_files.sh ${organism} ${genome} ${mode}

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
        echo "05: Resubmitting 05 and increasing time by 2 hrs (${sample})" | tee -a main.log
        wt=$(( wt + 2 ))
        qsub -l walltime=${wt}:00:00 -v \
sample=${sample},\
wt=${wt},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome},\
aln_only=${aln_only} \
${pipeline_dir}/05_run_bqsr.gatk.BaseRecalibrator.sh
        exit 0
    fi
else

# run gatl's BaseRecalibrator
$gatk_path/gatk --java-options "-Xmx10G -XX:+UseParallelGC -XX:ParallelGCThreads=12 -Djava.io.tmpdir=./.tmp" BaseRecalibrator \
 -R $reference \
 -L $intervals \
 --known-sites $knownsites_snps \
 --known-sites $knownsites_indels \
 -I preprocessed_bam/${sample}.markdup.bam \
 -O ${sample}.baserecalibrator.txt

# use samtools to increase the compression level in the read-corrected output bam
$gatk_path/gatk --java-options "-Xmx10G -XX:+UseParallelGC -XX:ParallelGCThreads=12 -Dsamjdk.compression_level=6 -Djava.io.tmpdir=./.tmp" ApplyBQSR \
 --bqsr ${sample}.baserecalibrator.txt \
 -I preprocessed_bam/${sample}.markdup.bam \
 -O ${dir}/${sample}.bqsr.bam \
 --create-output-bam-index

#/dev/stdout | samtools view -@ 12 -hb --write-index --no-PG -o BQSR/${sample}.final.bam##idx##BQSR/${sample}.final.bam.bai -

fi

# collect last command status
check_finish=$?

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
    echo "05: BQSR has been completed for sample ${sample}." | tee -a main.log
    if [[ ${aln_only} == 0 ]]; then
        # submit varscan
        echo "05: submitting pileups for Varscan ${sample}." | tee -a main.log
        ls $bed30intervals | grep ".bed" | parallel --tmpdir ./.tmp "qsub -v \
sample=${sample},\
bed={},\
index={#},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/06d_call_SNVs_and_indels.samtools.pileup.sh" | tee -a main.log
        # submit GetPileupSummaries
        if [[ ! -e contamination/${sample}.getpileupsummaries.table ]]; then
            echo "05: submitting GetPileupSummaries for CalculateContamination ${sample}." | tee -a main.log
            sample_pid_gps=$(qsub -v \
sample=${sample},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/06a_check_crosscontamination.gatk.GetPileupSummaries.sh)
            echo "05: GPS jobid for ${sample}: ${sample_pid_gps}" | tee -a main.log
        else
            echo "05: GetPileupSummaries table found for ${sample}: contamination/${sample}.getpileupsummaries.table" | tee -a main.log
        fi
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
                    check_normal=2
                    if [[ "${normal}" == "" || "${normal}" == "NA" || "${normal}" == "PON" || "${normal}" == "pon" ]]; then
                        normal="PON"
                    else
                        check_normal=$(samtools quickcheck ${dir}/${normal}.bqsr.bam && echo 1)
                    fi
                    # check if both bams are ready
                    if [[ "${check_normal}" == 1 && "${check_tumor}" == 1 ]]; then
                        # submit all mutect2 jobs
                        export normal
                        export tumor

                        # find if CalculateContamination or GetPileupSummaries is needed

                        if [[ -e contamination/${tumor}__${normal}.calculatecontamination.table ]]; then
                            echo "05: CalculateContamination table found for ${sample}: contamination/${sample}.calculatecontamination.table" | tee -a main.log
                        elif [[ -e contamination/${tumor}.getpileupsummaries.table && -e contamination/${normal}.getpileupsummaries.table ]]; then
                            echo "05: GetPileupSummaries table found for ${tumor}: contamination/${tumor}.getpileupsummaries.table" | tee -a main.log
                            echo "05: GetPileupSummaries table found for ${normal}: contamination/${normal}.getpileupsummaries.table" | tee -a main.log
                            # submit CalculateContamination
                            qsub -v \
normal=${normal},\
tumor=${tumor},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/06b_check_crosscontamination.gatk.CalculateContamination.sh | tee -a main.log
                            echo "05: Submitting CalculateContamination step" | tee -a main.log
                        # get pids of tumor and normal
                        elif [[ -e ${normal}.GetPileupSummaries.log && -e contamination/${tumor}.getpileupsummaries.table ]]; then
                            normal_pid_gps=$( head -1 ${normal}.GetPileupSummaries.log )
                            # submit second crosscheck as dependency for both normal and tumor
                            qsub -W depend=afterok:${normal_pid_gps} -v \
normal=${normal},\
tumor=${tumor},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/06b_check_crosscontamination.gatk.CalculateContamination.sh
                            echo "05: submitting dependency on normal:${normal} GPS for CalculateContamination ${sample}." | tee -a main.log
                        # if GPS finished already submit as only tumor dependency
                        elif [[ -e contamination/${normal}.getpileupsummaries.table ]]; then
                            qsub -W depend=afterok:${sample_pid_gps} -v \
normal=${normal},\
tumor=${tumor},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/06b_check_crosscontamination.gatk.CalculateContamination.sh
                            echo "05: submitting dependency on tumor:${tumor} GPS for CalculateContamination ${sample}." | tee -a main.log

                         else
                            # get jobid info from main.log
                            tumor_pid_gps=$(cat main.log | grep -a -o " ${tumor}: .*$" | grep -o "[0-9]*")
                            normal_pid_gps=$(cat main.log | grep -a -o " ${normal}: .*$" | grep -o "[0-9]*")
                            # submit job dep
                            qsub -W depend=afterok:${tumor_pid_gps}:${tumor_pid_gps} -v \
normal=${normal},\
tumor=${tumor},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/06b_check_crosscontamination.gatk.CalculateContamination.sh
                            echo "05: Submitting CalculateContamination step with dependency on both tumor:${tumor} and normal:${normal} GPS." | tee -a main.log
                        fi

                        # submit Mutect2 scattered runs
                        # save a dry run of commands first
                        ls $bed30intervals | grep ".bed" | parallel --tmpdir ./.tmp --dry-run "qsub -v \
normal=${normal},\
tumor=${tumor},\
bed={},\
index={#},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/06c_call_SNVs_and_indels.gatk.mutect2.sh" > all_logfiles/${tumor}__${normal}.mutect2.0.log
                        # submit mutect2 jobs on 30 intervals
                        ls $bed30intervals | grep ".bed" | parallel --tmpdir ./.tmp "qsub -v \
normal=${normal},\
tumor=${tumor},\
bed={},\
index={#},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/06c_call_SNVs_and_indels.gatk.mutect2.sh" | tee -a main.log
                    # if no normal (i.e., normal is PON)
                    elif [[ "${normal}" == "PON" ]]; then
                        # submit dependency job for GPS in tumor only
                        qsub -W depend=afterok:${sample_pid_gps} -v \
normal=${normal},\
tumor=${tumor},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/06b_check_crosscontamination.gatk.CalculateContamination.sh
                        # dry run before submitting mutect2 runs
                        ls $bed30intervals | grep ".bed" | parallel --tmpdir ./.tmp --dry-run "qsub -v \
normal=${normal},\
tumor=${tumor},\
bed={},\
index={#},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/06c_call_SNVs_and_indels.gatk.mutect2.sh" > all_logfiles/${tumor}__${normal}.mutect2.0.log
                        # submit mutect2 jobs on 30 intervals
                        ls $bed30intervals | grep ".bed" | parallel --tmpdir ./.tmp "qsub -v \
normal=${normal},\
tumor=${tumor},\
bed={},\
index={#},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/06c_call_SNVs_and_indels.gatk.mutect2.sh" | tee -a main.log
                    else
                        # GPS has not started for normal, wait for file
                        # wait for the BQSR script to finish
                        echo "05: ${sample} (tumor) waiting for ${normal} (normal) GPS to start / BQSR to finish." | tee -a main.log
                        qsub -v \
file="${normal}.GetPileupSummaries.log",\
sample=${sample},\
script=05_run_bqsr.gatk.BaseRecalibrator.sh,\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/wait_for_file.sh
                        exit 0
                    fi
                    # if no errors move logfile
                    if [[ "$?" == 0 ]]; then
                        mv ${tumor}.BQSR.log ${tumor}.baserecalibrator.txt all_logfiles
                        # log to main
                        echo "05: Mutect2 has started successfully for ${tumor}__${normal}." | tee -a main.log
                    fi
                done
            else
                echo "05: sample $sample not in tumor column (1st)" | tee -a main.log
                if [[ -e ${sample}.BQSR.log ]]; then
                    mv ${sample}.BQSR.log ${sample}.baserecalibrator.txt all_logfiles
                fi
                exit 0
            fi
        else
            echo -e "05: could not find \"tumors_and_normals.csv\" file" | tee -a main.log
            exit 1
        fi
    else
        echo "05: cleaning up..." | tee -a main.log
        if [[ $(samtools quickcheck ${dir}/${sample}.bqsr.bam && echo 1) == 1 ]]; then
            rm aligned_bam/${sample}.merged.ba*
            rm preprocessed_bam/${sample}.markdup.ba*
        fi
        echo "05: alignment-only mode finished for ${sample}." | tee -a main.log
        mv ${sample}.BQSR.log ${sample}.baserecalibrator.txt all_logfiles
    fi
fi
