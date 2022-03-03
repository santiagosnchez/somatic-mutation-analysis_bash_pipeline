#!/bin/bash
#PBS -l nodes=1:ppn=10,vmem=30g,mem=30g,walltime=12:00:00
#PBS -e ${sample}.VarScan.pileup.${index}.log
#PBS -j eo
# scheduler settings

# load modules
module load parallel/20210322
#module load sambamba/0.7.0
module load samtools/1.10

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# create output dirs
if [[ ! -e varscan ]]; then
    mkdir -p varscan/pileups
fi

# create tmp dir
if [[ ! -e tmp ]]; then
    mkdir tmp
fi

# set bam dir
if [[ ! -e bam ]]; then
    dir=BQSR
else
    dir=bam
fi

# load reference path and other reference files
# for details check script
source /hpf/largeprojects/tabori/shared/software/somatic-mutation-discovery/export_paths_to_reference_files.sh

if [[ $(samtools quickcheck ${dir}/${sample}.bqsr.bam && echo 1) == 1 ]]; then
# run mpileup in parallel for 30 intervals
    samtools mpileup --reference ${reference} -l $bed30intervals/$bed ${dir}/${sample}.bqsr.bam > varscan/pileups/${tumor}.${index}.pileup
    # check if finished
    check_finish=$?
else
    # log
    echo "06: ${sample} resubmitting to previous." | tee -a main.log
    # resubmit previous
    qsub -v sample=${sample},mode=${mode} ${pipeline_dir}/05_run_bqsr.gatk.BaseRecalibrator.sh
fi

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    # log to main
    echo "06: ${sample} Pileup for ${index} completed." | tee -a main.log
    # check if any pileup finished
    ls all_logfiles/${sample}.VarScan.pileup.*.log &> /dev/null
    if [[ "$?" == 0 ]]; then
        # check how many pileups finished
        pileups=$(ls all_logfiles/${sample}.VarScan.pileup.*.log | wc -l)
        if [[ "${pileups}" == 29 ]]; then
            cat $(ls varscan/pileups/${tumor}.*.pileup | sort -V) > varscan/pileups/${sample}.pileup && \
            rm varscan/pileups/${tumor}.*.pileup
            # check if both tumor and normal pileups are found
            cat tumors_and_normals.csv | grep "^${sample},"
            if [[ "$?" == 0 ]]; then
                for line in `cat tumors_and_normals.csv | grep "^${sample},"`; do
                    # first element is tumor, second is normal
                    tumor=$(echo $line | sed 's/,.*//')
                    normal=$(echo $line | sed 's/^.*,//')
                    if [[ -e varscan/pileups/${tumor}.pileup && -e varscan/pileups/${normal}.pileup ]]; then
                        # submit calling step
                        qsub -v tumor=${tumor},normal=${normal},mode=${mode} ${pipeline_dir}/06e_call_SNVs_and_indels.varscan.sh
                        # mv log and merge pileup logs
                        mv ${sample}.VarScan.pileup.${index}.log all_logfiles
                        cat $(ls all_logfiles/${sample}.VarScan.pileup.*.log | sort -V) > all_logfiles/${sample}.VarScan.pileup.log
                        # delete
                        rm all_logfiles/${sample}.VarScan.pileup.*.log
                    elif [[ -e varscan/pileups/${tumor}.pileup && ! -e varscan/pileups/${normal}.pileup ]]; then
                        # log
                        echo "06: waiting for normal ${normal} pileup to finish." | tee -a main.log
                        # wait for file
                        qsub -v file="varscan/pileups/${normal}.pileup",tumor=${tumor},normal=${normal},mode=${mode},script=06e_call_SNVs_and_indels.varscan.sh ${pipeline_dir}/wait_for_file.sh
                        # mv log and merge pileup logs
                        mv ${sample}.VarScan.pileup.${index}.log all_logfiles
                        cat $(ls all_logfiles/${sample}.VarScan.pileup.*.log | sort -V) > all_logfiles/${sample}.VarScan.pileup.log
                        # delete
                        rm all_logfiles/${sample}.VarScan.pileup.*.log
                    fi
                done
            else
                # if sample is normal, just tidyup
                normal=$(grep ",${sample}$" tumors_and_normals.csv && echo ${sample})
                if [[ ${sample} == ${normal} ]]; then
                    echo "06: pileup finished for normal ${normal} interval ${index}." | tee -a main.log
                    # mv log and merge pileup logs
                    mv ${sample}.VarScan.pileup.${index}.log all_logfiles
                    cat $(ls all_logfiles/${sample}.VarScan.pileup.*.log | sort -V) > all_logfiles/${sample}.VarScan.pileup.log
                    # delete
                    rm all_logfiles/${sample}.VarScan.pileup.*.log
                    echo "06: pileup finished for normal sample ${normal}." | tee -a main.log
                else
                    echo "06: sample ${sample} is not in the tumors_and_normals.csv file."
                    exit 1
                fi
            fi
        else
            echo "06: pileup finished for normal ${normal} interval ${index}." | tee -a main.log
            # mv log and merge pileup logs
            mv ${sample}.VarScan.pileup.${index}.log all_logfiles
        fi
    else
        echo "06: pileup finished for normal ${normal} interval ${index}." | tee -a main.log
        # mv log and merge pileup logs
        mv ${sample}.VarScan.pileup.${index}.log all_logfiles
    fi
fi
