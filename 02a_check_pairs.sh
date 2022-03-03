#PBS -l nodes=1:ppn=10,vmem=60g,mem=60g
#PBS -e ${sample}.${index}.checkpairs.log
#PBS -j eo
# scheduler settings


# load modules
module load sambamba/0.7.0
module load samtools/1.10
module load bwa/0.7.17

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# make directories
# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# create tmp dir
if [[ ! -e tmp ]]; then
    mkdir tmp
fi

# get walltime if not set
if [[ -z $wt ]]; then
    wt=$(qstat -f $PBS_JOBID | sed -rn 's/.*Resource_List.walltime = (.*)/\1/p' | sed 's/:.*//')
fi

# load reference path and other reference files
# for details check script
if [[ -z ${pipeline_dir} ]]; then
    source /hpf/largeprojects/tabori/shared/software/somatic-mutation-discovery/export_paths_to_reference_files.sh
else
    source ${pipeline_dir}/export_paths_to_reference_files.sh
fi

# check length of file names and run script
if [[ ${#forward} -gt 0 && ${#reverse} -gt 0 ]]; then
   # run script
   ${pipeline_dir}/fetch_fwd_rev_sing.pl ${forward} ${reverse} ${index} ${sample} | tee -a main.log
   if [[ "$?" == 0 ]]; then
       # get total singletons
       total_single=$(cat ${sample}.${index}.checkpairs.log | grep "Single: " | sed 's/.*: //')
       if [[ "${total_single}" -gt 0 ]]; then
          new_forward="tmp/${sample}.${index}.1.fastq.gz"
          new_reverse="tmp/${sample}.${index}.2.fastq.gz"
          singletons="tmp/${sample}.${index}.S.fastq.gz"
          # make tmp file_list file
          echo "${sample},${new_forward},${new_reverse}" >> tmp/${sample}_file_list.csv
          echo "${sample},${singletons}," >> tmp/${sample}_file_list.csv
          # make sure it's not duplicating
          cat tmp/${sample}_file_list.csv | sort -u > tmp/${sample}_file_list2.csv && mv tmp/${sample}_file_list2.csv tmp/${sample}_file_list.csv
          # calculate new walltime and read group
          wt=$(get_walltime $new_forward $new_reverse)
          # submit new jobs
          qsub -l walltime="${wt}":00:00 -v wt="${wt}",file_list="tmp/${sample}_file_list.csv",index=${index},sample=${sample},forward=${new_forward},reverse=${new_reverse},mode=${mode} ${pipeline_dir}/02b_align_and_sort_bam_to_ref.bwa.sh
          # for singletons
          wt=$(get_walltime $singletons)
          qsub -l walltime="${wt}":00:00 -v wt="${wt}",file_list="tmp/${sample}_file_list.csv",index="${index}s",sample=${sample},forward=${singletons},mode=${mode} ${pipeline_dir}/02b_align_and_sort_bam_to_ref.bwa.sh
      else
          # delete tmp files
          new_forward="tmp/${sample}.${index}.1.fastq.gz"
          new_reverse="tmp/${sample}.${index}.2.fastq.gz"
          singletons="tmp/${sample}.${index}.S.fastq.gz"
          rm $new_forward $new_reverse $singletons
          # proceed normally
          qsub -l walltime="${wt}":00:00 -v wt="${wt}",file_list="$file_list",index=${index},sample=${sample},forward=${forward},reverse=${reverse},mode=${mode} ${pipeline_dir}/02b_align_and_sort_bam_to_ref.bwa.sh
      fi
   else
      echo "fetch_fwd_rev_sing.pl failed with an error for ${sample}"
      exit 1
   fi
else
    # submit as single ended
     qsub -l walltime="${wt}":00:00 -v wt="${wt}",file_list="$file_list",index=${index},sample=${sample},forward=${forward},mode=${mode} ${pipeline_dir}/02b_align_and_sort_bam_to_ref.bwa.sh
fi

# final check
if [[ "$?" == 0 ]]; then
    # move log
    mv ${sample}.${index}.checkpairs.log all_logfiles
    echo "02: Done checking pairs for ${sample}.${index}" | tee -a main.log
fi
