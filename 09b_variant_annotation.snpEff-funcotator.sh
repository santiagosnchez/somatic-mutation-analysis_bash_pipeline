#!/bin/bash
#PBS -l nodes=1:ppn=1,vmem=10g,mem=10g,walltime=8:00:00
#PBS -e ${tumor}__${normal}.annotation-snpeff-funcotator.${tissue}.log
#PBS -j eo
# scheduler settings

# set date to calculate running time
start=$(date)

# load modules
module load java/1.8
module load snpEff/4.11
module load vep/102
module load bcftools/1.11
module load tabix
module load parallel/20210322

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# create output dirs
if [[ ! -e vcf/snpEff ]]; then
    mkdir -p vcf/snpEff
fi
if [[ ! -e vcf/vep ]]; then
    mkdir -p vcf/vep
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

# if [[ ! -e "varscan/${tumor}__${normal}.varscan.all.Germline.hc.${mode}.vcf.gz" ]]; then
#     echo "09: VarScan has not finished for ${tumor}__${normal}. Waiting..." | tee -a main.log
#     qsub -v file="varscan/${tumor}__${normal}.varscan.all.Germline.hc.${mode}.vcf.gz",tumor=${tumor},normal=${normal},mode=${mode},script=09_variant_annotation.snpEff-funcotator.sh ${pipeline_dir}/wait_for_file.sh
#     exit 0
# fi


if [[ "${tissue}" == "Somatic" ]]; then
    # what caller
    caller="mutect2"
    # make temporary header file with pedigree
    echo "##PEDIGREE=<Derived=${tumor},Original=${normal}>" > ./.tmp/${tumor}__${normal}.tmp.vcf.header.txt

    # normalize variants, reduce complex alleles and add header line
    bcftools annotate \
     -h ./.tmp/${tumor}__${normal}.tmp.vcf.header.txt \
     -o ./.tmp/${tumor}__${normal}.${caller}.normalized_head.${mode}.vcf \
     ${caller}/${tumor}__${normal}.${caller}.filtered-norm.${mode}.vcf.gz

    # delete tmp header file
    if [[ "$?" == 0 ]]; then
        rm ./.tmp/${tumor}__${normal}.tmp.vcf.header.txt
        # index tmp
        index-vcf ./.tmp/${tumor}__${normal}.${caller}.normalized_head.${mode}.vcf
        # replace prev vcf
        mv ./.tmp/${tumor}__${normal}.${caller}.normalized_head.${mode}.vcf.gz \
           ${caller}/${tumor}__${normal}.${caller}.filtered-norm.${mode}.vcf.gz
        mv ./.tmp/${tumor}__${normal}.${caller}.normalized_head.${mode}.vcf.gz.tbi \
           ${caller}/${tumor}__${normal}.${caller}.filtered-norm.${mode}.vcf.gz.tbi
    fi

    # run annotators on selected only
    bcftools view -f PASS ${caller}/${tumor}__${normal}.${caller}.filtered-norm.${mode}.vcf.gz > ${caller}/${tumor}__${normal}.${caller}.selected.${mode}.vcf
    index-vcf ${caller}/${tumor}__${normal}.${caller}.selected.${mode}.vcf

    # run annotators on selected only, on no-obpriors
    bcftools view -f PASS ${caller}/${tumor}__${normal}.${caller}.filtered_no-obpriors-norm.${mode}.vcf.gz > ${caller}/${tumor}__${normal}.${caller}.selected_no-obpriors.${mode}.vcf
    index-vcf ${caller}/${tumor}__${normal}.${caller}.selected_no-obpriors.${mode}.vcf

    if [[ -e ${snpeff_datadir}/${genome} ]]; then
        # check if sample has matched normal
        if [[ "${normal}" != "PON" ]]; then
            # run snpEff on mutect2 vcf
            java -jar $snpeff_jar \
             -dataDir $snpeff_datadir \
             ${genome} \
             -v \
             -canon \
             -cancer \
             -stats vcf/snpEff/${tumor}__${normal}.${tissue}.snpEff_summary.html \
             -csvStats vcf/snpEff/${tumor}__${normal}.${tissue}.snpEff_summary.csv \
             ${caller}/${tumor}__${normal}.${caller}.selected.${mode}.vcf.gz > \
             vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-snpeff.${mode}.vcf
             # on no_ob
             java -jar $snpeff_jar \
              -dataDir $snpeff_datadir \
              ${genome} \
              -v \
              -canon \
              -cancer \
              -stats vcf/snpEff/${tumor}__${normal}.${tissue}.snpEff_summary.html \
              -csvStats vcf/snpEff/${tumor}__${normal}.${tissue}.snpEff_summary.csv \
              ${caller}/${tumor}__${normal}.${caller}.selected_no-obpriors.${mode}.vcf.gz > \
              vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-snpeff_no-obpriors.${mode}.vcf
        else
            # annotate unmatched on PON
            # run snpEff on mutect2 vcf
            java -jar $snpeff_jar \
             -dataDir $snpeff_datadir \
             ${genome} \
             -v \
             -canon \
             -stats vcf/snpEff/${tumor}__${normal}.${tissue}.snpEff_summary.html \
             -csvStats vcf/snpEff/${tumor}__${normal}.${tissue}.snpEff_summary.csv \
             ${caller}/${tumor}__${normal}.${caller}.selected.${mode}.vcf.gz > \
             vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-snpeff.${mode}.vcf
            # on no-ob
            java -jar $snpeff_jar \
             -dataDir $snpeff_datadir \
             ${genome} \
             -v \
             -canon \
             -stats vcf/snpEff/${tumor}__${normal}.${tissue}.snpEff_summary.html \
             -csvStats vcf/snpEff/${tumor}__${normal}.${tissue}.snpEff_summary.csv \
             ${caller}/${tumor}__${normal}.${caller}.selected_no-obpriors.${mode}.vcf.gz > \
             vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-snpeff_no-obpriors.${mode}.vcf

        fi
         # index vcf
         index-vcf vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-snpeff.${mode}.vcf
         index-vcf vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-snpeff_no-obpriors.${mode}.vcf
    else
         ln -s ${caller}/${tumor}__${normal}.${caller}.selected.${mode}.vcf.gz vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-snpeff.${mode}.vcf.gz
         ln -s ${caller}/${tumor}__${normal}.${caller}.selected.${mode}.vcf.gz.tbi vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-snpeff.${mode}.vcf.gz.tbi
         ln -s ${caller}/${tumor}__${normal}.${caller}.selected_no-obpriors.${mode}.vcf.gz vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-snpeff_no-obpriors.${mode}.vcf.gz
         ln -s ${caller}/${tumor}__${normal}.${caller}.selected_no-obpriors.${mode}.vcf.gz.tbi vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-snpeff_no-obpriors.${mode}.vcf.gz.tbi
    fi
     if [[ ! -z ${funcotator_databases_s} ]]; then
         # run gatk's funcotator on somatic mutations
         $gatk_path/gatk Funcotator \
          --variant ${caller}/${tumor}__${normal}.${caller}.selected.${mode}.vcf.gz \
          --reference $reference \
          --ref-version ${genome} \
          --data-sources-path $funcotator_databases_s \
          --transcript-selection-mode CANONICAL \
          --output vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-funcotator.${mode}.vcf \
          --output-file-format VCF
          # index vcf
          index-vcf vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-funcotator.${mode}.vcf
          # run gatk's funcotator on somatic mutations
          $gatk_path/gatk Funcotator \
           --variant ${caller}/${tumor}__${normal}.${caller}.selected_no-obpriors.${mode}.vcf.gz \
           --reference $reference \
           --ref-version ${genome} \
           --data-sources-path $funcotator_databases_s \
           --transcript-selection-mode CANONICAL \
           --output vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-funcotator_no-obpriors.${mode}.vcf \
           --output-file-format VCF
           # index vcf
           index-vcf vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-funcotator_no-obpriors.${mode}.vcf
    fi
    # if [[ ! -z ${vep_datadir} ]]; then
    ##   zcat ${caller}/${tumor}__${normal}.${caller}.selected.${mode}.vcf.gz | \
    #     vep -i ${caller}/${tumor}__${normal}.${caller}.selected.${mode}.vcf \
    #      -v \
    #      --cache \
    #      --dir_cache ${vep_datadir} \
    #      --refseq \
    #      --species ${vep_species} \
    #      -o vcf/vep/${tumor}__${normal}.vep.txt \
    #      --stats_file vcf/vep/${tumor}__${normal}.vep_summary.html \
    #      --canonical \
    #      --vcf \
    #      --stats_text
    # fi
fi

elif [[ "${tissue}" == "Germline" ]]; then
    # define caller
    #caller="varscan"
    caller="haplotypecaller"

    if [[ -e ${snpeff_datadir}/${genome} ]]; then
        if [[ -e ${caller}/haplotypecaller/${tumor}__${normal}.${caller}.selected.${mode}.vcf.gz ]]; then
          # run snpEff on Varscan Germline calls
          java -jar $snpeff_jar \
           -dataDir $snpeff_datadir \
           ${genome} \
           -v \
           -canon \
           -stats vcf/snpEff/${tumor}__${normal}.${tissue}.snpEff_summary.html \
           -csvStats vcf/snpEff/${tumor}__${normal}.${tissue}.snpEff_summary.csv \
           ${caller}/${tumor}__${normal}.${caller}.selected.${mode}.vcf.gz > \
           vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-snpeff.${mode}.vcf
           # ${caller}/${tumor}__${normal}.${caller}.all.${tissue}.hc.${mode}.vcf.gz
           # index vcf
           index-vcf vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-snpeff.${mode}.vcf
           # log
           echo "09: Done annotating with SnpEff ${caller}/${tumor}__${normal}.${caller}.selected.${mode}.vcf.gz" | tee -a main.log
       else
          # log error
          echo "09: HaplotypeCaller VCF (${caller}/${tumor}__${normal}.${caller}.selected.${mode}.vcf.gz) not found."
          check_finish=1
       fi
    else
         ln -s ${caller}/${tumor}__${normal}.${caller}.selected.${mode}.vcf.gz vcf/${tumor}__${normal}.${caller}.all.${tissue}.selected.${mode}.vcf.gz
         ln -s ${caller}/${tumor}__${normal}.${caller}.selected.${mode}.vcf.gz.tbi vcf/${tumor}__${normal}.${caller}.all.${tissue}.selected.${mode}.vcf.gz.tbi
    fi
    if [[ ! -z ${funcotator_databases_g} ]]; then
        # run gatk's funcotator on germline variants
        $gatk_path/gatk Funcotator \
         --variant ${caller}/${tumor}__${normal}.${caller}.selected.${mode}.vcf.gz \
         --reference $reference \
         --ref-version ${genome} \
         --data-sources-path $funcotator_databases_g \
         --transcript-selection-mode CANONICAL \
         --output vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-funcotator.${mode}.vcf \
         --output-file-format VCF
         # index vcf
         index-vcf vcf/${tumor}__${normal}.${caller}.all.${tissue}.annotated-funcotator.${mode}.vcf
    else
        echo "09: Funcotator database not found for sample ${tumor}__${normal}." | tee -a main.log
        echo "09: Organism is ${organism}." | tee -a main.log
    fi
else
    echo "09: caller was not defined (mutect2/haplotypecaller) for ${tumor}__${normal}" | tee -a main.log
    exit 1
fi
# get maf from funcotator vcf

#ls vcf/${tumor}__${normal}.*annotated-funcotator* | parallel --plus '${pipeline_dir}/funcotator-vcf2maf.sh {} > {/.vcf/.maf}'

# check if finished
check_finish=$?

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    # log to main
    echo "09: ${tissue} annotation with SnpEff and Funcotator completed for ${tumor}__${normal}." | tee -a main.log
    # run analyses
    # submit last step
    qsub -v \
normal=${normal},\
tumor=${tumor},\
tissue=${tissue},\
mode=${mode},\
pipeline_dir=${pipeline_dir},\
organism=${organism},\
genome=${genome} \
${pipeline_dir}/10_run_analyses.signatures_and_TBM.sh
    # calc runtime
    runtime=$( how_long "${start}" h )
    echo "09: ${tumor}__${normal}.annotation-snpeff-funcotator.${tissue}.log took ${runtime} hours" | tee -a main.log
    # move logfile
    mv ${tumor}__${normal}.annotation-snpeff-funcotator.${tissue}.log all_logfiles
fi
