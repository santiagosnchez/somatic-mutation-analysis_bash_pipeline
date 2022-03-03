#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=10g,mem=10g,walltime=5:00:00
#PBS -e ${tumor}__${normal}.annotation.${tissue}.log
#PBS -j eo
# scheduler settings

# load modules
module load java/1.8
module load snpEff/4.11
module load bcftools/1.11
module load tabix
module load parallel/20210322

# set working dir
cd $PBS_O_WORKDIR

# create output dirs
if [[ ! -e vcf/snpEff ]]; then
    mkdir -p vcf/snpEff
fi

# create tmp dir
if [[ ! -e tmp ]]; then
    mkdir tmp
fi

# load reference path and other reference files
# for details check script
if [[ -z ${pipeline_dir} ]]; then
    source /hpf/largeprojects/tabori/shared/software/somatic-mutation-discovery/export_paths_to_reference_files.sh
else
    source ${pipeline_dir}/export_paths_to_reference_files.sh
fi

# change intervals to null if not WES
if [[ "${mode}" != "wes" ]]; then
    intervals=null
fi

# if [[ ! -e "varscan/${tumor}__${normal}.varscan.all.Germline.hc.${mode}.vcf.gz" ]]; then
#     echo "09: VarScan has not finished for ${tumor}__${normal}. Waiting..." | tee -a main.log
#     qsub -v file="varscan/${tumor}__${normal}.varscan.all.Germline.hc.${mode}.vcf.gz",tumor=${tumor},normal=${normal},mode=${mode},script=09_variant_annotation.snpEff-funcotator.sh ${pipeline_dir}/wait_for_file.sh
#     exit 0
# fi

if [[ "${tissue}" == "Somatic" ]]; then
    # what caller
    caller="mutect2"
    # make temporary header file with pedigree
    echo "##PEDIGREE=<Derived=${tumor},Original=${normal}>" > ${tumor}__${normal}.tmp.vcf.header.txt

    # normalize variants, reduce complex alleles and add header line
    bcftools annotate \
     -h ${tumor}__${normal}.tmp.vcf.header.txt \
     -Oz -o mutect2/${tumor}__${normal}.${caller}.normalized_head.${mode}.vcf.gz \
     ${caller}/${tumor}__${normal}.${caller}.normalized.${mode}.vcf.gz

    # delete tmp header file
    if [[ "$?" == 0 ]]; then
        rm ${tumor}__${normal}.tmp.vcf.header.txt
    fi

    # run snpEff on mutect2 vcf
    java -jar $snpeff_jar \
     -dataDir $snpeff_datadir \
     hg38 \
     -v \
     -canon \
     -cancer \
     -stats vcf/snpEff/${tumor}__${normal}.${tissue}.snpEff_summary.html \
     -csvStats vcf/snpEff/${tumor}__${normal}.${tissue}.snpEff_summary.csv \
     ${caller}/${tumor}__${normal}.${caller}.normalized_head.${mode}.vcf.gz > \
     vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-snpeff.${mode}.vcf

     # run gatk's funcotator on somatic mutations
     $gatk_path/gatk Funcotator \
      --variant ${caller}/${tumor}__${normal}.${caller}.selected.${mode}.vcf.gz \
      --reference $reference \
      --ref-version hg38 \
      --data-sources-path $funcotator_databases_s \
      --transcript-selection-mode CANONICAL \
      --output vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-funcotator.${mode}.vcf \
      --output-file-format VCF

elif [[ "${tissue}" == "Germline" ]]; then
    # define caller
    caller="varscan"
    # run snpEff on Varscan Germline calls
    java -jar $snpeff_jar \
     -dataDir $snpeff_datadir \
     hg38 \
     -v \
     -canon \
     -stats vcf/snpEff/${tumor}__${normal}.${tissue}.snpEff_summary.html \
     -csvStats vcf/snpEff/${tumor}__${normal}.${tissue}.snpEff_summary.csv \
     ${caller}/${tumor}__${normal}.${caller}.all.${tissue}.hc.${mode}.vcf.gz  > \
     vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-snpeff.${mode}.vcf

    # run gatk's funcotator on germline variants
    $gatk_path/gatk Funcotator \
     --variant ${caller}/${tumor}__${normal}.${caller}.all.${tissue}.hc.${mode}.vcf.gz \
     --reference $reference \
     --ref-version hg38 \
     --data-sources-path $funcotator_databases_g \
     --transcript-selection-mode CANONICAL \
     --output vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-funcotator.${mode}.vcf \
     --output-file-format VCF
else
    echo "09: caller was not define (mutect2/varscan) for ${tumor}__${normal}" | tee -a main.log
    exit 1
fi
# get maf from funcotator vcf

#ls vcf/${tumor}__${normal}.*annotated-funcotator* | parallel --plus '${pipeline_dir}/funcotator-vcf2maf.sh {} > {/.vcf/.maf}'

# check if finished
check_finish=$?

# create log dir
if [[ ! -e all_logfiles ]]; then
    mkdir all_logfiles
fi

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    # bgzip and tabix vcf
    ls vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated*.vcf | parallel index-vcf {}
    # log to main
    echo "09: ${tissue} annotation with SnpEff and Funcotator completed for ${tumor}__${normal}." | tee -a main.log
    # run analyses
    if [[ "${tissue}" == "Somatic" && -e "all_logfiles/${tumor}__${normal}.annotation.Germline.log" ]]; then
        # submit last step
        qsub -v normal=${normal},tumor=${tumor},mode=${mode} ${pipeline_dir}/10_run_analyses.signatures_and_TBM.sh
    elif [[ "${tissue}" == "Germline" && -e "all_logfiles/${tumor}__${normal}.annotation.Somatic.log" ]]; then
        # submit last step
        qsub -v normal=${normal},tumor=${tumor},mode=${mode} ${pipeline_dir}/10_run_analyses.signatures_and_TBM.sh
    fi
    # move logfile
    mv ${tumor}__${normal}.annotation.${tissue}.log all_logfiles
fi
