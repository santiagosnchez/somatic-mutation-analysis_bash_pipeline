#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=10g,mem=10g,walltime=5:00:00
#PBS -e ${tumor}__${normal}.annovar.${tissue}.log
#PBS -j eo
# scheduler settings

# set date to calculate running time
start=$(date)

# load modules
module load bcftools/1.11
module load tabix
module load parallel/20210322

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# create output dirs
if [[ ! -e annovar ]]; then
    mkdir annovar
fi

# create tmp dir
if [[ ! -e .tmp ]]; then
    mkdir .tmp
fi

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# load reference path and other reference files
# for details check script
source ${pipeline_dir}/00_export_pipeline_environment.sh ${organism} ${genome} ${mode}

# check that mutect2annovar.pl script is available
mutect2annovar.pl --help &> /dev/null
if [[ "$?" != 1 ]]; then # exits at 1 instead of 0 when calling --help
    echo "09: mutect2annovar.pl was not found in path or perl module not found." | tee -a main.log
    exit 1
fi

# check that mutect2annovar.pl script is available
table_annovar.pl --help &> /dev/null
if [[ "$?" != 1 ]]; then # exits at 1 instead of 0 when calling --help
    echo "09: table_annovar.pl was not found in path." | tee -a main.log
    exit 1
fi

# if [[ ! -e "varscan/${tumor}__${normal}.varscan.all.Germline.hc.${mode}.vcf.gz" ]]; then
#     echo "09: VarScan has not finished for ${tumor}__${normal}. Waiting..." | tee -a main.log
#     qsub -v file="varscan/${tumor}__${normal}.varscan.all.Germline.hc.${mode}.vcf.gz",tumor=${tumor},normal=${normal},mode=${mode},script=09_variant_annotation.snpEff-funcotator.sh ${pipeline_dir}/wait_for_file.sh
#     exit 0
# fi

if [[ "${tissue}" == "Somatic" ]]; then
    # what caller
    caller="mutect2"
    # make annovar table
    mutect2annovar.pl \
      --vcf ${caller}/${tumor}__${normal}.${caller}.filtered-norm.${mode}.vcf.gz \
      --output annovar/${tumor}__${normal}.${caller}.${tissue}.filtered-norm.${mode}.mutect2annovar_tbl.txt \
      --filter false \
      --header false \
      --tumour ${tumor} \
      --normal ${normal}
    # run annovar
    # fetch the file name of the bedfile (for both WGS and WES)
    bedfile=$(echo ${intervals_bed} | rev | cut -d/ -f1 | rev)
    # make a symlink not ready
    if [[ ! -e $annovar_db/${bedfile} ]]; then
        ln -s ${intervals_bed} $annovar_db/${bedfile}
    fi
    # run annovar script
    table_annovar.pl \
      annovar/${tumor}__${normal}.${caller}.${tissue}.filtered-norm.${mode}.mutect2annovar_tbl.txt \
      $annovar_db \
      --protocol refGene,ensGene,avsnp150,1000g2015aug_all,esp6500siv2_all,cosmic70,clinvar_20220320,exac03,bed \
      --operation g,g,f,f,f,f,f,f,r \
      --buildver ${genome} \
      --remove \
      --otherinfo \
      --bedfile ${bedfile} \
      --outfile annovar/${tumor}__${normal}.${caller}.${tissue}.filtered-norm.${mode}

elif [[ "${tissue}" == "Germline" ]]; then
    # exit, not for germline
    echo "09: annovar not set up for germline."
    exit 1
else
    echo "09: caller was not defined (mutect2/varscan) for ${tumor}__${normal}" | tee -a main.log
    exit 1
fi
# get maf from funcotator vcf

#ls vcf/${tumor}__${normal}.*annotated-funcotator* | parallel --plus '${pipeline_dir}/funcotator-vcf2maf.sh {} > {/.vcf/.maf}'

# check if finished
check_finish=$?

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    echo "09: Annovar finished for ${tumor}__${normal}." | tee -a main.log
    # calc runtime
    runtime=$( how_long "${start}" h )
    echo "09: Step ${tumor}__${normal}.annovar.${tissue}.log took ${runtime} hours" | tee -a main.log
    # move logfile
    mv ${tumor}__${normal}.annovar.${tissue}.log all_logfiles
fi
