#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=30g,mem=30g,walltime=15:00:00
#PBS -e ${tumor}__${normal}.mutect2.${index}.log
#PBS -j eo
# scheduler settings

# set date to calculate running time
start=$(date)

# load modules
module load java/1.8
#module load gatk/4.2.2.0
module load samtools/1.10
module load bcftools/1.11
module load tabix

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# log start
echo $start

# create dir for contamination
if [[ -z $make_pon ]]; then
    if [[ ! -e contamination ]]; then
        mkdir contamination
    fi
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

# check if make_pon exists
if [[ ${make_pon} == 1 ]]; then
# run MuTect2 for makeing PoN
$gatk_path/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./.tmp" Mutect2 \
 -I ${dir}/${normal}.bqsr.bam \
 -R ${reference} \
 --max-mnp-distance 0 \
 -O mutect2/${normal}.mutect2.pon.${index}.vcf \
 -L ${bed30intervals}/${bed}
else
  # run Mutect2 on tumor
  if [[ ! -e mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf ]]; then
    if [[ "${normal}" == "PON" ]]; then
      # if specified, use own PoN
      grep -a "Using specified PoN" main.log &> /dev/null
      if [[ $? == 0 ]]; then
          # define specific PoN
          this_pon=$(grep -a "Using specified PoN" main.log | sed 's/.*: //')
          gatk_pon=${gatk_pon_location}/${this_pon}
      fi
$gatk_path/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./.tmp" Mutect2 \
 -I ${dir}/${tumor}.bqsr.bam \
 -tumor ${tumor} \
 -R ${reference} \
 -O mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.${index}.vcf \
 -germline-resource $gnomad_resource \
 -pon ${gatk_pon} \
 --max-mnp-distance 0 \
 --f1r2-tar-gz mutect2/f1r2/${tumor}__${normal}.${index}.f1r2.tar.gz \
 --genotype-germline-sites true \
 --genotype-pon-sites true \
 -L ${bed30intervals}/${bed}
    elif [[ ${gnomad_resource} == 'null' && -e ${gatk_pon} ]]; then
$gatk_path/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./.tmp" Mutect2 \
 -I ${dir}/${tumor}.bqsr.bam \
 -I ${dir}/${normal}.bqsr.bam \
 -tumor ${tumor} \
 -normal ${normal} \
 -R ${reference} \
 -O mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.${index}.vcf \
 -pon ${gatk_pon} \
 --max-mnp-distance 0 \
 --f1r2-tar-gz mutect2/f1r2/${tumor}__${normal}.${index}.f1r2.tar.gz \
 --genotype-germline-sites true \
 --genotype-pon-sites true \
 -L ${bed30intervals}/${bed}
    elif [[ ${gnomad_resource} == 'null' && ! -e ${gatk_pon} ]]; then
$gatk_path/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./.tmp" Mutect2 \
 -I ${dir}/${tumor}.bqsr.bam \
 -I ${dir}/${normal}.bqsr.bam \
 -tumor ${tumor} \
 -normal ${normal} \
 -R ${reference} \
 -O mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.${index}.vcf \
 --max-mnp-distance 0 \
 --f1r2-tar-gz mutect2/f1r2/${tumor}__${normal}.${index}.f1r2.tar.gz \
 --genotype-germline-sites true \
 -L ${bed30intervals}/${bed}
    else
# submit GetPileupSummaries for normal and tumor
    # tumor first
      if [[ ! (-e contamination/${tumor}.getpileupsummaries.table || -e contamination/${tumor}.getpileupsummaries.${index}.table) ]]; then
$gatk_path/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./.tmp" GetPileupSummaries \
 -I ${dir}/${tumor}.bqsr.bam \
 -V ${gnomad_resource} \
 -L $bed30intervals/${bed} \
 -O contamination/${tumor}.getpileupsummaries.${index}.table
        echo "06: Done with GPS for ${tumor} interval ${index}" |  tee -a main.log
      fi
    # then normal
      if [[ ! (-e contamination/${normal}.getpileupsummaries.table || -e contamination/${normal}.getpileupsummaries.${index}.table) ]]; then
$gatk_path/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./.tmp" GetPileupSummaries \
 -I ${dir}/${normal}.bqsr.bam \
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
 --genotype-germline-sites true \
 --genotype-pon-sites true \
 -L ${bed30intervals}/${bed}
 #  -bamout mutect2/${tumor}__${normal}.${index}.bam \
 #  --create-output-bam-index \
    fi
  else
      # force $? == 0
      ls &> /dev/null
  fi
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
        if [[ "${mutect_logfiles}" == 29 ]]; then
            # gather vcffiles
            # generate list of files with their own -I flag
            if [[ ${make_pon} == 1 ]]; then
                vcffiles=$(ls mutect2/${normal}.mutect2.pon.*.vcf | sort -V | sed 's/^/-I /')
                $gatk_path/gatk GatherVcfs $vcffiles -O mutect2/${normal}.mutect2.pon.merged.vcf
                # delete if finished
                if [[ "$?" == 0 ]]; then
                    rm mutect2/${normal}.mutect2.pon.[1-9]*.vcf
                    rm mutect2/${normal}.mutect2.pon.[1-9]*.vcf.idx
                fi
                # gather stats files, needed for Filtering
                statsfiles=$(ls mutect2/${normal}.mutect2.pon.*.vcf.stats | sort -V | sed 's/^/-stats /')
                $gatk_path/gatk MergeMutectStats $statsfiles -O mutect2/${normal}.mutect2.pon.merged.vcf.stats
                # delete if finished
                if [[ "$?" == 0 ]]; then
                    rm mutect2/${normal}.mutect2.pon.[1-9]*.vcf.stats
                fi
                # compress and index
                index-vcf mutect2/${normal}.mutect2.pon.merged.vcf
            else
                vcffiles=$(ls mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.*.vcf | sort -V | sed 's/^/-I /')
                $gatk_path/gatk GatherVcfs $vcffiles -O mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf
                # delete if finished
                if [[ "$?" == 0 ]]; then
                    rm mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.[1-9]*.vcf
                    rm mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.[1-9]*.vcf.idx
                fi
                # gather stats files, needed for Filtering
                statsfiles=$(ls mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.*.vcf.stats | sort -V | sed 's/^/-stats /')
                $gatk_path/gatk MergeMutectStats $statsfiles -O mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf.stats
                # delete if finished
                if [[ "$?" == 0 ]]; then
                    rm mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.[1-9]*.vcf.stats
                fi
            fi
            # log to main
            echo "06: ${tumor}__${normal} Mutect2 variant calling completed for interval ${index}." | tee -a main.log
            # check if make_pon exists
            if [[ ${make_pon} == 1 ]]; then
                echo "06: VCF for PoN created successfuly for ${normal}." | tee -a main.log
                # get samples
                file_list=$(grep -a "^file list: " main.log | tail -1 | sed 's/^file list: //')
                samples=$(cat ${file_list} | cut -d, -f1 | sort -u)
                total_samples=$(echo "${samples}" | wc -l)
                # count VCF files
                total_vcfs=$(ls mutect2/*.mutect2.pon.merged.vcf.gz | grep "${samples/ /\|/}" | wc -l)
                if [[ ${total_samples} == ${total_vcfs} ]]; then
                    # submit last PoN job
                    qsub -v \
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/07_pon.gatk.CreateSomaticPanelOfNormals.sh
                    # last log
                    echo "06: All VCFs completed." | tee -a main.log
                    echo "06: Submitting final step for PoN." | tee -a main.log
                else
                    remaining=$(echo "${total_samples}-${total_vcfs}" | bc)
                    echo "06: $remaining VCF(s) remaining."
                fi
            else
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
                if [[ "${normal}" != "PON" && ${gnomad_resource} != "null" ]]; then
                    # gather GPS tables merged
                    gpsfiles_tumor=$(ls contamination/${tumor}.getpileupsummaries.*.table | sort -V | sed 's/^/-I /')
                    $gatk_path/gatk GatherPileupSummaries ${gpsfiles_tumor} -O contamination/${tumor}.getpileupsummaries.table --sequence-dictionary ${reference_dict}
                    if [[ "$?" == 0 ]]; then
                        rm contamination/${tumor}.getpileupsummaries.[1-9]*.table
                    fi
                    gpsfiles_normal=$(ls contamination/${normal}.getpileupsummaries.*.table | sort -V | sed 's/^/-I /')
                    $gatk_path/gatk GatherPileupSummaries ${gpsfiles_normal} -O contamination/${normal}.getpileupsummaries.table --sequence-dictionary ${reference_dict}
                    if [[ "$?" == 0 ]]; then
                        rm contamination/${normal}.getpileupsummaries.[1-9]*.table
                    fi
                    # log to main
                    echo "06: GetPileupSummaries completed ${tumor} and ${normal}." | tee -a main.log
                    echo "06: Moving to CalculateContamination for ${tumor}__${normal}." | tee -a main.log
                else
                    echo "06: Normal is PON or no gnomad resource (null) for ${tumor}." | tee -a main.log
                fi
                qsub -v \
normal=${normal},\
tumor=${tumor},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/06c_check_crosscontamination.gatk.CalculateContamination.sh | tee -a main.log
            fi
            # estimate runtime
            # first scatter
            first_scatter_date=$(ls ${tumor}__${normal}.mutect2.${index}.log all_logfiles/${tumor}__${normal}.mutect2.[1-9]*.log | \
              parallel 'head -2 {} | tail -1' | parallel date --date={} +%s | sort -n | parallel date --date=@{} | head -1)
            # calc runtime
            runtime=$( how_long "${first_scatter_date}" h )
            # log
            echo "06: ${tumor}__${normal} Mutect2 took ${runtime} hours" | tee -a main.log
            # concat logfiles
            cat $(ls ${tumor}__${normal}.mutect2.${index}.log all_logfiles/${tumor}__${normal}.mutect2.[0-9]*.log | sort -V) > all_logfiles/${tumor}__${normal}.mutect2.log
            rm $(ls all_logfiles/${tumor}__${normal}.mutect2.[0-9]*.log )
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
