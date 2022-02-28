#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=30g,mem=30g,walltime=12:00:00
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

# create output dirs
if [[ ! -e mutect2 ]]; then
    mkdir -p mutect2/f1r2
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
# change intervals to null if not WES
if [[ "${mode}" != "wes" ]]; then
    intervals=null
fi

if [[ ! -e mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf ]]; then
# run gatk's mutect2
$gatk_path/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./tmp" Mutect2 \
 -I ${dir}/${tumor}.bqsr.bam \
 -I ${dir}/${normal}.bqsr.bam \
 -tumor ${tumor} \
 -normal ${normal} \
 -R ${reference} \
 -O mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.${index}.vcf \
 -germline-resource $gnomad_resource \
 --f1r2-tar-gz mutect2/f1r2/${tumor}__${normal}.${index}.f1r2.tar.gz \
 -L ${bed30intervals}/${bed}
 #  -bamout mutect2/${tumor}__${normal}.${index}.bam \
 #  --create-output-bam-index \
else
 ls &> /dev/null
fi

# check if finished
check_finish=$?

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    # next round of jobs are submitted manually or not
    # check if all mutect2 operations finished
    # first check for files
    ls all_logfiles/${tumor}__${normal}.mutect2.[1-9]*.log &> /dev/null
    # if mutect2 still running
    if [[ "$?" == 0 ]]; then
        mutect_logfiles=$(ls all_logfiles/${tumor}__${normal}.mutect2.[1-9]*.log | wc -l)
        # try to wrap up in one go
        if [[ ! -e all_logfiles/${tumor}__${normal}.mutect2.log && "${mutect_logfiles}" == 29 ]]; then
            cat $(ls all_logfiles/${tumor}__${normal}.mutect2.[0-9]*.log | sort -V) > all_logfiles/${tumor}__${normal}.mutect2.log
            rm $(ls all_logfiles/${tumor}__${normal}.mutect2.[0-9]*.log )
            # gather vcffiles
            # generate list of files with their own -I flag
            vcffiles=$(ls mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.*.vcf | sort -V | sed 's/^/-I /')
            $gatk_path/gatk GatherVcfs $vcffiles -O mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf
            # gather stats files, needed for Filtering
            statsfiles=$(ls mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.*.vcf.stats | sort -V | sed 's/^/-stats /')
            $gatk_path/gatk MergeMutectStats $statsfiles -O mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf.stats
            if [[ "$?" == 0 ]]; then
                rm mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.[1-9]*.vcf*
            fi
            # log to main
            echo "06: ${tumor}__${normal} Mutect2 variant calling completed for interval ${index}." | tee -a main.log
            echo "06: Moving to read-orientation for ${tumor}__${normal}." | tee -a main.log
            # submit read orientation analysis
            qsub -v tumor=${tumor},normal=${normal},mode=${mode} ${pipeline_dir}/07_read_orientation.gatk.LearnReadOrientationModel.sh
            # move logfile
            mv ${tumor}__${normal}.mutect2.${index}.log all_logfiles
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
        else
            # check what's missing
            if [[ ! -e mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf ]]; then
                # gather vcffiles
                # generate list of files with their own -I flag
                vcffiles=$(ls mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.*.vcf | sort -V | sed 's/^/-I /')
                $gatk_path/gatk GatherVcfs $vcffiles -O mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf
                # gather stats files, needed for Filtering
                statsfiles=$(ls mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.*.vcf.stats | sort -V | sed 's/^/-stats /')
                $gatk_path/gatk MergeMutectStats $statsfiles -O mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf.stats
                if [[ "$?" == 0 ]]; then
                    rm mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.[1-9]*.vcf*
                fi
            elif [[ ! -e mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf.stats ]]; then
                # gather stats files, needed for Filtering
                statsfiles=$(ls mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.*.vcf.stats | sort -V | sed 's/^/-stats /')
                $gatk_path/gatk MergeMutectStats $statsfiles -O mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.merged.vcf.stats
                if [[ "$?" == 0 ]]; then
                    rm mutect2/${tumor}__${normal}.mutect2.unfiltered.${mode}.[1-9]*.vcf*
                fi
            fi
            # log to main
            echo "06: ${tumor}__${normal} Mutect2 variant calling completed." | tee -a main.log
            echo "06: Moving to read-orientation for ${tumor}__${normal}." | tee -a main.log
            # submit read orientation analysis
            qsub -v tumor=${tumor},normal=${normal},mode=${mode} ${pipeline_dir}/07_read_orientation.gatk.LearnReadOrientationModel.sh
            # move logfile
            mv ${tumor}__${normal}.mutect2.${index}.log all_logfiles
        fi
    fi
fi
