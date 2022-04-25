#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=30g,mem=30g,walltime=10:00:00
#PBS -e ${tumor}__${normal}.mutect2.${index}.log
#PBS -j eo
# scheduler settings

# load modules
module load java/1.8
#module load gatk/4.2.2.0
module load samtools/1.10

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# create dir for contamination
if [[ ! -e contamination ]]; then
    mkdir contamination
fi

# create output dirs
if [[ ! -e mutect2 ]]; then
    mkdir -p mutect2/f1r2
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

if [[ ! -e mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf ]]; then
  if [[ "${normal}" == "PON" ]]; then
$gatk_path/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./.tmp" Mutect2 \
 -I ${dir}/${tumor}.bqsr.bam \
 -tumor ${tumor} \
 -R ${reference} \
 -O mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.${index}.vcf \
 -germline-resource $gnomad_resource \
 -pon ${gatk_pon} \
 --max-mnp-distance 0 \
 --f1r2-tar-gz mutect2/f1r2/${tumor}__${normal}.${index}.f1r2.tar.gz \
 -L ${bed30intervals}/${bed}
  elif [[ ${gnomad_resource} == 'none' ]]; then
$gatk_path/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./.tmp" Mutect2 \
 -I ${dir}/${tumor}.bqsr.bam \
 -tumor ${tumor} \
 -R ${reference} \
 -O mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.${index}.vcf \
 -pon ${gatk_pon} \
 --max-mnp-distance 0 \
 --f1r2-tar-gz mutect2/f1r2/${tumor}__${normal}.${index}.f1r2.tar.gz \
 -L ${bed30intervals}/${bed}
  else
# submit GetPileupSummaries for normal and tumor
    # tumor first
    if [[ ! -e contamination/${tumor}.getpileupsummaries.table ]]; then
$gatk_path/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./.tmp" GetPileupSummaries \
 -I ${dir}/${tumor}.bqsr.bam \
 -V ${gnomad_resource} \
 -L $bed30intervals/${bed} \
 -O contamination/${tumor}.getpileupsummaries.${index}.table
      echo "06: Done with GPS for ${tumor} interval ${index}" |  tee -a main.log
    fi
    # then normal
    if [[ ! -e contamination/${normal}.getpileupsummaries.table ]]; then
$gatk_path/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./.tmp" GetPileupSummaries \
 -I ${dir}/${tumor}.bqsr.bam \
 -V ${gnomad_resource} \
 -L $bed30intervals/${bed} \
 -O contamination/${normal}.getpileupsummaries.${index}.table
      echo "06: Done with GPS for ${normal} interval ${index}" |  tee -a main.log
    fi

# run gatk's mutect2
$gatk_path/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./.tmp" Mutect2 \
 -I ${dir}/${tumor}.bqsr.bam \
 -I ${dir}/${normal}.bqsr.bam \
 -tumor ${tumor} \
 -normal ${normal} \
 -R ${reference} \
 -O mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.${index}.vcf \
 -germline-resource $gnomad_resource \
 -pon ${gatk_pon} \
 --max-mnp-distance 0 \
 --f1r2-tar-gz mutect2/f1r2/${tumor}__${normal}.${index}.f1r2.tar.gz \
 -L ${bed30intervals}/${bed}
 #  -bamout mutect2/${tumor}__${normal}.${index}.bam \
 #  --create-output-bam-index \
  fi
else
    # force $? == 0
    ls &> /dev/null
fi

# check if finished
check_finish=$?

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    # check if all mutect2 operations finished
    # first check for files
    ls all_logfiles/${tumor}__${normal}.mutect2.[1-9]*.log &> /dev/null
    # if mutect2 still running
    if [[ "$?" == 0 ]]; then
        mutect_logfiles=$(ls all_logfiles/${tumor}__${normal}.mutect2.[1-9]*.log | wc -l)
        # try to wrap up in one go
        if [[ ! -e all_logfiles/${tumor}__${normal}.mutect2.log && "${mutect_logfiles}" == 29 ]]; then
            cat $(ls ${tumor}__${normal}.mutect2.${index}.log all_logfiles/${tumor}__${normal}.mutect2.[0-9]*.log | sort -V) > all_logfiles/${tumor}__${normal}.mutect2.log
            rm $(ls all_logfiles/${tumor}__${normal}.mutect2.[0-9]*.log )
            # gather vcffiles
            # generate list of files with their own -I flag
            vcffiles=$(ls mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.*.vcf | sort -V | sed 's/^/-I /')
            $gatk_path/gatk GatherVcfs $vcffiles -O mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf
            # gather stats files, needed for Filtering
            statsfiles=$(ls mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.*.vcf.stats | sort -V | sed 's/^/-stats /')
            $gatk_path/gatk MergeMutectStats $statsfiles -O mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf.stats
            # log to main
            if [[ "$?" == 0 ]]; then
                rm mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.[1-9]*.vcf*
            fi
            # log to main
            echo "06: ${tumor}__${normal} Mutect2 variant calling completed for interval ${index}." | tee -a main.log
            echo "06: Moving to read-orientation for ${tumor}__${normal}." | tee -a main.log
            # submit read orientation analysis
            qsub -v \
tumor=${tumor},\
normal=${normal},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/07_read_orientation.gatk.LearnReadOrientationModel.sh
            # do CalculateContamination for non tumor-only runs
            if [[ ${normal} != "PON" ]]; then
                # gather GPS tables merged
                gpsfiles_tumor=$(ls contamination/${tumor}.getpileupsummaries.*.table | sort -V | sed 's/^/-I /')
                $gatk_path/gatk GatherPileupSummaries ${gpsfiles_tumor} -O contamination/${tumor}.getpileupsummaries.table --sequence-dictionary ${reference_dict}
                gpsfiles_normal=$(ls contamination/${normal}.getpileupsummaries.*.table | sort -V | sed 's/^/-I /')
                $gatk_path/gatk GatherPileupSummaries ${gpsfiles_normal} -O contamination/${normal}.getpileupsummaries.table --sequence-dictionary ${reference_dict}
                if [[ "$?" == 0 ]]; then
                    rm contamination/${tumor}.getpileupsummaries.[1-9]*.table
                    rm contamination/${normal}.getpileupsummaries.[1-9]*.table
                fi
                # log to main
                echo "06: GetPileupSummaries completed ${tumor} and ${normal}." | tee -a main.log
                echo "06: Moving to CalculateContamination for ${tumor}__${normal}." | tee -a main.log
            else
                echo "06: Normal is PON for ${tumor}." | tee -a main.log
            fi
            qsub -v \
normal=${normal},\
tumor=${tumor},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/06c_check_crosscontamination.gatk.CalculateContamination.sh | tee -a main.log
            # move logfile
            rm ${tumor}__${normal}.mutect2.${index}.log
        else
            # log to main
            echo "06: ${tumor}__${normal} Mutect2 variant calling completed for interval ${index}." | tee -a main.log
            # move logfile
            mv ${tumor}__${normal}.mutect2.${index}.log all_logfiles
        fi
    # no scattered logfiles found
    else
        # check if mutect2 is running
        ls ${tumor}__${normal}.mutect2.[1-9]*.log &> /dev/null
        if [[ "$?" == 0 ]]; then
            # log to main
            echo "06: ${tumor}__${normal} Mutect2 variant calling completed for interval ${index}." | tee -a main.log
            # move logfile
            mv ${tumor}__${normal}.mutect2.${index}.log all_logfiles
        fi
    fi
fi
