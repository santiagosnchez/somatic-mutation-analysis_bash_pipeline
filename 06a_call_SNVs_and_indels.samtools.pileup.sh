#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=30g,mem=30g,walltime=12:00:00
#PBS -e ${sample}.VarScan.pileup.${index}.log
#PBS -j eo
# scheduler settings

# set date to calculate running time
start=$(date)

# load modules
module load parallel/20210322
#module load sambamba/0.7.0
module load samtools/1.10

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# print start date
echo $start

# create output dirs
if [[ ! -e varscan ]]; then
    mkdir -p varscan/pileups
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

if [[ $(samtools quickcheck ${dir}/${sample}.bqsr.bam && echo 1) == 1 ]]; then
# run mpileup in parallel for 30 intervals
    samtools mpileup --reference ${reference} -l $bed30intervals/$bed ${dir}/${sample}.bqsr.bam > varscan/pileups/${sample}.${index}.pileup
    # check if finished
    check_finish=$?
else
    # log
    echo "06: ${sample} bam not found or truncated." | tee -a main.log
    # fail signal
    check_finish=1
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
        if [[ "${pileups}" == 29 || "${pileups}" == 30 ]]; then
            cat $(ls varscan/pileups/${sample}.*.pileup | sort -V) > varscan/pileups/${sample}.pileup
            if [[ "$?" == 0 ]]; then
                rm varscan/pileups/${sample}.*.pileup
                # check if both tumor and normal pileups are found
                cat tumors_and_normals.csv | grep "^${sample},"
                if [[ "$?" == 0 ]]; then
                    for line in `cat tumors_and_normals.csv | grep "^${sample},"`; do
                        # first element is tumor, second is normal
                        tumor=$(echo $line | sed 's/,.*//')
                        normal=$(echo $line | sed 's/^.*,//')
                        if [[ "${normal}" == "" || "${normal}" == "NA" || "${normal}" == "PON" || "${normal}" == "pon" ]]; then
                            echo "06: VarScan does not run on tumor-only mode for ${tumor}. Skipping..." | tee -a main.log
                        else
                            if [[ -e varscan/pileups/${tumor}.pileup && -e varscan/pileups/${normal}.pileup ]]; then
                                # log
                                echo "06: all pileup finished for tumor ${tumor} and normal ${normal}." | tee -a main.log
                                echo "06: submitting VarScan calls." | tee -a main.log
                                # estimate runtime
                                # first scatter
                                first_scatter_date=$(ls ${sample}.VarScan.pileup.${index}.log all_logfiles/${sample}.VarScan.pileup.[0-9]*.log | \
                                  parallel 'head -2 {} | tail -1' | parallel date --date={} +%s | sort -n | parallel date --date=@{} | head -1)
                                runtime=$( how_long "${first_scatter_date}" h )
                                # calc runtime
                                echo "06: VarScan pileup took ${runtime} hours for ${sample}" | tee -a main.log
                                # log to main
                                # submit calling step
                                qsub -v \
tumor=${tumor},\
normal=${normal},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/06b_call_SNVs_and_indels.varscan.sh
                            elif [[ -e varscan/pileups/${tumor}.pileup && ! -e varscan/pileups/${normal}.pileup ]]; then
                                # log
                                echo "06: waiting for normal ${normal} pileup to finish." | tee -a main.log
                                # wait for file
                                qsub -v \
file="varscan/pileups/${normal}.pileup",\
sample=${tumor},\
tumor=${tumor},\
normal=${normal},\
script=06b_call_SNVs_and_indels.varscan.sh,\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/wait_for_file.sh
                            fi
                        fi
                    done
                    # mv log and merge pileup logs
                    cp ${sample}.VarScan.pileup.${index}.log all_logfiles
                    cat $(ls all_logfiles/${sample}.VarScan.pileup.*.log | sort -V) > all_logfiles/${sample}.VarScan.pileup.log
                    # delete
                    rm all_logfiles/${sample}.VarScan.pileup.*.log
                    rm ${sample}.VarScan.pileup.${index}.log
                else
                    # if sample is normal, just tidyup
                    normal=$(grep ",${sample}$" tumors_and_normals.csv &> /dev/null && echo ${sample})
                    if [[ ${sample} == ${normal} ]]; then
                        echo "06: all pileup finished for normal ${normal}." | tee -a main.log
                        # mv log and merge pileup logs
                        cp ${sample}.VarScan.pileup.${index}.log all_logfiles
                        cat $(ls all_logfiles/${sample}.VarScan.pileup.*.log | sort -V) > all_logfiles/${sample}.VarScan.pileup.log
                        # delete
                        rm all_logfiles/${sample}.VarScan.pileup.*.log
                        echo "06: pileup finished for normal sample ${normal}." | tee -a main.log
                        rm ${sample}.VarScan.pileup.${index}.log
                    else
                        echo "06: sample ${sample} is not in the tumors_and_normals.csv file."
                        exit 1
                    fi
                fi
            else
                echo "06: error merging pileups ${sample}." | tee -a main.log
                exit 1
            fi
        else
            echo "06: pileup finished for ${sample} interval ${index}." | tee -a main.log
            # mv log and merge pileup logs
            mv ${sample}.VarScan.pileup.${index}.log all_logfiles
        fi
    else
        # finished scatter run
        echo "06: pileup finished for ${sample} interval ${index}." | tee -a main.log
        # mv log and merge pileup logs
        mv ${sample}.VarScan.pileup.${index}.log all_logfiles
    fi
fi
