#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=10g,mem=10g,walltime=5:00:00
#PBS -e ${tumor}__${normal}.annovar.${tissue}.log
#PBS -j eo
# scheduler settings

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
source ${pipeline_dir}/export_paths_to_reference_files.sh ${organism} ${genome} ${mode}

# check that mutect2annovar.pl script is available
mutect2annovar.pl --help &> /dev/null
if [[ "$?" != 1 ]]; then # exits at 1 instead of 0 when calling --help
    echo "09: mutect2annovar.pl was not found in path or perl module not found."
    exit 1
fi

# check that mutect2annovar.pl script is available
table_annovar.pl --help &> /dev/null
if [[ "$?" != 1 ]]; then # exits at 1 instead of 0 when calling --help
    echo "09: table_annovar.pl was not found in path."
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
    # extract SNPs
    # do on all variants
    # bcftools view -v snps -Oz \
    #   ${caller}/${tumor}__${normal}.${caller}.filtered.${mode}.vcf.gz > \
    #   annovar/${tumor}__${normal}.${caller}.snv.${tissue}.filtered.${mode}.vcf.gz
    # make annovar table
    mutect2annovar.pl \
      --vcf mutect2/${tumor}__${normal}.${caller}.filtered-norm.${mode}.vcf.gz \
      --output annovar/${tumor}__${normal}.${caller}.${tissue}.filtered-norm.${mode}.mutect2annovar_tbl.txt \
      --filter false \
      --header false \
      --tumour ${tumor} \
      --normal ${normal}
    # run annovar
    # fetch the file name of the bedfile (for both WGS and WES)
    bedfile=$(echo ${intervals_bed} | rev | cut -d/ -f1 | rev)
    # make a symlink
    ln -sf ${intervals_bed} $annovar_db/${bedfile}
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
    # bgzip and tabix vcf
    ls vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated*.vcf | parallel index-vcf {}
    # log to main
    echo "09: ${tissue} annotation with SnpEff and Funcotator completed for ${tumor}__${normal}." | tee -a main.log
    # run analyses
    if [[ "${tissue}" == "Somatic" && -e "all_logfiles/${tumor}__${normal}.annotation.Germline.log" ]]; then
        # submit last step
        qsub -v \
normal=${normal},\
tumor=${tumor},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/10_run_analyses.signatures_and_TBM.sh
    elif [[ "${tissue}" == "Germline" && -e "all_logfiles/${tumor}__${normal}.annotation.Somatic.log" ]]; then
        # submit last step
        qsub -v \
normal=${normal},\
tumor=${tumor},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/10_run_analyses.signatures_and_TBM.sh
    fi
    # move logfile
    mv ${tumor}__${normal}.annotation.${tissue}.log all_logfiles
fi