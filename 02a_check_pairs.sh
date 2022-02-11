#PBS -l nodes=1:ppn=10,vmem=60g,mem=60g
#PBS -e ${sample}.checkpairs.log
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

#
if [[ ${#forward} -gt 0 && ${#reverse} -gt 0 ]]; then
   ${pipeline_dir}/fetch_fwd_rev_sing.pl ${forward} ${reverse} | tee -a main.log
   if [[ "$?" == 0 ]]; then
       total_single=$(cat ${sample}.checkpairs.log | grep "Single: " | sed 's/.*: //')
       if [[ "${total_single}" -gt 0 ]]; then
          new_forward=$(cat ${sample}.checkpairs.log | grep "tmp/.*.1.fastq.gz")
          new_reverse=$(cat ${sample}.checkpairs.log | grep "tmp/.*.2.fastq.gz")
          singletons=$(cat ${sample}.checkpairs.log | grep "tmp/.*.S.fastq.gz")
          # make tmp file_list file
          echo "${sample},${new_forward},${new_reverse}" >> tmp/${sample}_file_list.csv
          echo "${sample},${singletons}" >> tmp/${sample}_file_list.csv
          # calculate new walltime
          wt=$(get_walltime $new_forward $new_reverse)
          # submit new jobs
          qsub -l walltime="${wt}":00:00 -v file_list="tmp/${sample}_file_list.csv",index=${index},sample=${sample},rg="${rg}",forward=${new_forward},reverse=${new_reverse},mode=${mode} ${pipeline_dir}/02b_align_and_sort_bam_to_ref.bwa.sh
          # for singletons
          wt=$(get_walltime $singletons)
          qsub -l walltime="${wt}":00:00 -v file_list="tmp/${sample}_file_list.csv",index="${index}s",sample=${sample},rg="${rg}",forward=${singletons},reverse="",mode=${mode} ${pipeline_dir}/02b_align_and_sort_bam_to_ref.bwa.sh
      else
          # proceed normally
          qsub -l walltime="${wt}":00:00 -v file_list="$file_list",index=${index},sample=${sample},rg="${rg}",forward=${forward},reverse=${reverse},mode=${mode} ${pipeline_dir}/02b_align_and_sort_bam_to_ref.bwa.sh
      fi
   else
      echo "fetch_fwd_rev_sing.pl failed with an error for ${sample}"
   fi
else
    # submit as single ended
     qsub -l walltime="${wt}":00:00 -v file_list="$file_list",index=${index},sample=${sample},rg="${rg}",forward=${forward},reverse="",mode=${mode} ${pipeline_dir}/02b_align_and_sort_bam_to_ref.bwa.sh
fi
