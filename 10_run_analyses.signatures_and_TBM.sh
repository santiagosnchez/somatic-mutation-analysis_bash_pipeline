#!/bin/bash
#PBS -l nodes=1:ppn=8,vmem=10g,mem=10g,walltime=5:00:00
#PBS -e ${tumor}__${normal}.analyses.${tissue}.log
#PBS -j eo
# scheduler settings

# set date to calculate running time
start=$(date)

# load modules
module load java/1.8
#module load gatk/4.0.1.2
module load samtools/1.10
module load bcftools/1.11
#module load R/4.1.0
#module load R/4.0.2
module load R/4.1.2
module load parallel/20210322

# set working dir
cd $PBS_O_WORKDIR

# print jobid to 1st line
echo $PBS_JOBID

# create output dirs
if [[ ! -e analyses ]]; then
    mkdir analyses
fi
if [[ ! -e analyses/signatures ]]; then
    mkdir analyses/signatures
fi
if [[ ! -e analyses/annotations ]]; then
    mkdir analyses/annotations
fi

if [[ ! -e analyses/no-obpriors ]]; then
    mkdir analyses/no-obpriors
fi
if [[ ! -e analyses/no-obpriors/signatures ]]; then
    mkdir analyses/no-obpriors/signatures
fi
if [[ ! -e analyses/no-obpriors/annotations ]]; then
    mkdir analyses/no-obpriors/annotations
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

# list of MMR genes:
# MLH1 MSH2 MSH6 PMS2 POLD1 POLE IDH1 TP53 NF1
# grep -w "MLH1\|MSH2\|MSH6\|PMS2\|POLD1\|POLE\|IDH1\|TP53\|NF1"

# pull germline and somatic missense (nonsynonymous) mutations
# look for MMR genes
# germline on VarScan calls

# ${pipeline_dir}/get_gene_annotations_from_vcf-funcotator.sh \
#  vcf/${tumor}__${normal}.varscan.all.Germline.annotated-funcotator.${mode}.vcf.gz \
#  MLH1 \
#  MSH2 \
#  MSH6 \
#  PMS2 \
#  POLD1 \
#  POLE \
#  IDH1 \
#  TP53 \
#  NF1 \
#  > analyses/${tumor}__${normal}.germline_MMR_mutations.genes.csv
# add IDH4

# ${pipeline_dir}/get_gene_annotations_from_vcf.sh \
# vcf/${tumor}__${normal}.varscan.all.Germline.annotated-snpeff.${mode}.vcf.gz \
# POLD1 \
# POLD2 \
# POLD3 \
# POLD4 \
# POLE \
# POLE2 > analyses/${tumor}__${normal}.germline_POL_mutations.genes.csv
#
# # somatic on Mutect2
# ${pipeline_dir}/get_gene_annotations_from_vcf-funcotator.sh \
#   vcf/${tumor}__${normal}.mutect2.all.Somatic.annotated-funcotator.${mode}.vcf.gz \
#   MLH1 \
#   MSH2 \
#   MSH6 \
#   PMS2 \
#   POLD1 \
#   POLE \
#   IDH1 \
#   TP53 \
#   NF1 \
#   > analyses/${tumor}__${normal}.somatic_MMR_mutations.genes.csv
#
# ${pipeline_dir}/get_gene_annotations_from_vcf.sh \
# vcf/${tumor}__${normal}.mutect2.annotated-snpeff.${mode}.vcf.gz \
# POLD1 \
# POLD2 \
# POLD3 \
# POLD4 \
# POLE \
# POLE2 > analyses/${tumor}__${normal}.somatic_POL_mutations.genes.csv

echo "10: Fetching all variant annotations." | tee -a main.log

# get all annotations into csv
# funcotator Somatic
if [[ ${tissue} == "Somatic" ]]; then
    if [[ -e vcf/${tumor}__${normal}.mutect2.all.Somatic.annotated-funcotator.${mode}.vcf.gz ]]; then
        bash ${pipeline_dir}/scripts/funcotator-vcf2maf2.sh \
        vcf/${tumor}__${normal}.mutect2.all.Somatic.annotated-funcotator.${mode}.vcf.gz \
        ${tumor} ${normal} Somatic \
        > analyses/annotations/${tumor}__${normal}.annotations_funcotator_somatic.csv
    fi
    if [[ -e vcf/${tumor}__${normal}.mutect2.all.Somatic.annotated-snpeff.${mode}.vcf.gz ]]; then
        # get all annotations into csv
        # snpeff Somatic
        bash ${pipeline_dir}/scripts/snpeff-vcf2tbl.sh \
        vcf/${tumor}__${normal}.mutect2.all.Somatic.annotated-snpeff.${mode}.vcf.gz \
        ${tumor} ${normal} Somatic \
        > analyses/annotations/${tumor}__${normal}.annotations_snpeff_somatic.csv
    fi
    # get stats for no-ob file annotated-snpeff_no-obpriors
    if [[ ${mode} != "wgs" ]]; then
        if [[ -e vcf/${tumor}__${normal}.mutect2.all.Somatic.annotated-funcotator_no-obpriors.${mode}.vcf.gz ]]; then
            bash ${pipeline_dir}/scripts/funcotator-vcf2maf2.sh \
            vcf/${tumor}__${normal}.mutect2.all.Somatic.annotated-funcotator_no-obpriors.${mode}.vcf.gz \
            ${tumor} ${normal} Somatic \
            > analyses/no-obpriors/annotations/${tumor}__${normal}.annotations_funcotator_somatic_no-obpriors.csv
        fi
        if [[ -e vcf/${tumor}__${normal}.mutect2.all.Somatic.annotated-snpeff_no-obpriors.${mode}.vcf.gz ]]; then
            # get all annotations into csv
            # snpeff Somatic
            bash ${pipeline_dir}/scripts/snpeff-vcf2tbl.sh \
            vcf/${tumor}__${normal}.mutect2.all.Somatic.annotated-snpeff_no-obpriors.${mode}.vcf.gz \
            ${tumor} ${normal} Somatic \
            > analyses/no-obpriors/annotations/${tumor}__${normal}.annotations_snpeff_somatic_no-obpriors.csv
        fi
    fi
    echo "10: done extracting somatic variants from vcf file for ${tumor}__${normal}" | tee -a main.log
fi

# check if tumor-only mode or germline
if [[ "${normal}" != "PON" && ${tissue} == "Germline" ]]; then
    if [[ -e vcf/${tumor}__${normal}.haplotypecaller.all.Germline.annotated-funcotator.${mode}.vcf.gz ]]; then
        bash ${pipeline_dir}/scripts/funcotator-vcf2maf2.sh \
        vcf/${tumor}__${normal}.haplotypecaller.all.Germline.annotated-funcotator.${mode}.vcf.gz \
        ${tumor} ${normal} Germline \
        > analyses/annotations/${tumor}__${normal}.annotations_funcotator_germline.csv
    fi
    if [[ -e vcf/${tumor}__${normal}.haplotypecaller.all.Germline.annotated-snpeff.${mode}.vcf.gz ]]; then
        bash ${pipeline_dir}/scripts/snpeff-vcf2tbl.sh \
        vcf/${tumor}__${normal}.haplotypecaller.all.Germline.annotated-snpeff.${mode}.vcf.gz \
        ${tumor} ${normal} Germline \
        > analyses/annotations/${tumor}__${normal}.annotations_snpeff_germline.csv
    fi
    echo "10: done extracting germline variants from vcf file for ${tumor}__${normal}" | tee -a main.log
fi

if [[ ${tissue} == "Somatic" ]]; then
    # log
    echo "10: Running signature analysis (${tumor}__${normal})" | tee -a main.log

    # run variant analysis
    Rscript ${pipeline_dir}/scripts/variant_analysis.nofigs2.R ${mode} ${tumor}__${normal} ${organism} "signatures"
    # for no-ob
    if [[ ${mode} != "wgs" ]]; then
        Rscript ${pipeline_dir}/scripts/variant_analysis.nofigs2.R ${mode} ${tumor}__${normal} ${organism} "no-obpriors/signatures"

        sig_prop_nob=$(echo $(cat analyses/no-obpriors/signatures/${tumor}__${normal}.COSMIC_v3.2.signatures.csv | \
        awk -v FS="," '$3 ~ /^(SBS1|SBS6|SBS10a|SBS10b|SBS10c|SBS10d|SBS14|SBS15|SBS20|SBS26|SBS44)$/' | \
        cut -d, -f1) | tr ' ' ',')
    fi

    # get sigs into a single CSV line
    # header first
    # sig_heads=$(echo $(cat analyses/signatures/${tumor}__${normal}.COSMIC_v3.2.signatures.csv | \
    # awk -v FS="," '$3 ~ /^(SBS1|SBS6|SBS10a|SBS10b|SBS10c|SBS10d|SBS14|SBS15|SBS20|SBS26|SBS44)$/' | \
    # cut -d, -f3) | tr ' ' ',')

    # proportions
    sig_prop=$(echo $(cat analyses/signatures/${tumor}__${normal}.COSMIC_v3.2.signatures.csv | \
    awk -v FS="," '$3 ~ /^(SBS1|SBS6|SBS10a|SBS10b|SBS10c|SBS10d|SBS14|SBS15|SBS20|SBS26|SBS44)$/' | \
    cut -d, -f1) | tr ' ' ',')

    echo "10: done with COSMIC signatures (${tumor}__${normal})" | tee -a main.log

    # log
    echo "10: calculating observed coverage, SNVs and indels (${tumor}__${normal})" | tee -a main.log

    # add header to analyses/coverage_tmb_and_mmr_sigs
    if [[ ! -e analyses/coverage_tmb_and_mmr_sigs.csv ]]; then
        echo "Tumor,Normal,Observed_coverage,Expected_coverage,SNV,Indels,TMB_SNV,TMB_indels,SBS1,SBS6,SBS10a,SBS10b,SBS10c,SBS10d,SBS14,SBS15,SBS20,SBS26,SBS44" > analyses/coverage_tmb_and_mmr_sigs.csv
    fi

    # add header to analyses/old_output.tmbs.tsv
    if [[ ! -e analyses/old_output.tmbs.tsv ]]; then
        echo -e "sample\ttotal_SNV\tTMB" > analyses/old_output.tmbs.tsv
    fi

    # expected coverage
    expected=$(cat $intervals_bed | awk '{ count = count + ($3 - ($2 + 1)) } END { print count }')
    # find observed coverage
    if [[ ! -e ${dir}/${tumor}.bqsr.bam ]]; then
        # if no bam, just use expected
        coverage=${expected}
    else
        if [[ "${normal}" != "PON" ]]; then
            coverage=$(samtools depth -b $intervals_bed -q20 -Q20 -d1000 ${dir}/${tumor}.bqsr.bam ${dir}/${normal}.bqsr.bam | awk '$3 >= 4 && $4 >= 4' | wc -l)
        else
            coverage=$(samtools depth -b $intervals_bed -q20 -Q20 -d1000 ${dir}/${tumor}.bqsr.bam | awk '$3 >= 4' | wc -l)
        fi
    fi

    # estimate tumor mutation burden (TMB)
    # use prev coverage estimate
    # total snvs
    total_snvs=$(bcftools view -H --types snps mutect2/${tumor}__${normal}.mutect2.selected.${mode}.vcf.gz | wc -l)

    # total indels
    total_indels=$(bcftools view -H --types indels mutect2/${tumor}__${normal}.mutect2.selected.${mode}.vcf.gz | wc -l)

    # calc TMB
    TMB_snvs=$( echo "scale=2; ${total_snvs}/(${coverage}/1000000)" | bc | sed 's/^\./0\./')
    TMB_indels=$( echo "scale=2; ${total_indels}/(${coverage}/1000000)" | bc | sed 's/^\./0\./' )

    # output
    echo "${tumor},${normal},${coverage},${expected},${total_snvs},${total_indels},${TMB_snvs},${TMB_indels},${sig_prop}" >> analyses/coverage_tmb_and_mmr_sigs.csv
    echo -e "${tumor}\t${normal}\t${total_snvs}\t${TMB_snvs}" >> analyses/old_output.tmbs.tsv

    # get stats for no-ob file
    if [[ ${mode} != "wgs" ]]; then

        total_snvs_nob=$(bcftools view -H --types snps mutect2/${tumor}__${normal}.mutect2.selected_no-obpriors.${mode}.vcf.gz | wc -l)
        total_indels_nob=$(bcftools view -H --types indels mutect2/${tumor}__${normal}.mutect2.selected_no-obpriors.${mode}.vcf.gz | wc -l)
        TMB_snvs_nob=$( echo "scale=2; ${total_snvs_nob}/(${coverage}/1000000)" | bc | sed 's/^\./0\./')
        TMB_indels_nob=$( echo "scale=2; ${total_indels_nob}/(${coverage}/1000000)" | bc | sed 's/^\./0\./' )

        # output
        if [[ ! -e analyses/no-obpriors/coverage_and_tmb_no-obpriors.csv ]]; then
            echo "tumor,normal,obs_coverage,exp_coverage,snvs,indels,tmb_snvs,tmb_indels,"${sig_heads} > analyses/no-obpriors/coverage_and_tmb_no-obpriors.csv
        else
            echo "${tumor},${normal},${coverage},${expected},${total_snvs_nob},${total_indels_nob},${TMB_snvs_nob},${TMB_indels_nob},${sig_prop_nob}" >> analyses/no-obpriors/coverage_tmb_and_mmr_sigs_no-obpriors.csv
        fi
    fi

    echo "10: Done calculating tumor mutation burden (TMB)." | tee -a main.log

    # look at differences in calls between varscan, mutect2 with ob-priors and without.
    # varscan_snvs=$(bcftools view -H -v snps varscan/${tumor}__${normal}.all.Somatic.hc.vcf.gz | wc -l)
    # mutect2_all_filters_snvs=$(bcftools view -H -v snps -f PASS mutect2/${tumor}__${normal}.mutect2.selected.${mode}.vcf | wc -l)
    # mutect2_all_filters_snvs=$(bcftools view -H -v snps -f PASS mutect2/${tumor}__${normal}.mutect2.selected_no-obpriors.${mode}.vcf | wc -l)

fi

# add to archive
# zip -ru ${tumor}__${normal}.analyses.zip analyses/${tumor}__${normal}.*

# check if finished
check_finish=$?

# check if command finished
if [[ "$check_finish" == 0 ]]; then
    # log finished sample
    echo "10: Done for ${tumor}__${normal} in ${tissue}" | tee -a main.log
    # check if finished for sample
    if [[ ${mode} == "wgs" || ${normal} == "PON" || (-e all_logfiles/${tumor}__${normal}.analyses.Germline.log && ${tissue} == "Somatic") || (-e all_logfiles/${tumor}__${normal}.analyses.Somatic.log && ${tissue} == "Germline") ]]; then
        # final logs and tidy up dir, ie. gargabe collection
        # run final COSMIC signature analysis
        echo "${tumor},${normal}" >> finished.csv
        # get time it took
        total_time_in_hours=$( how_long main.log h )
        total_time_in_days=$( how_long main.log d )
        echo "pipeline finished for ${tumor}__${normal} in ${total_time_in_hours} hrs or ${total_time_in_days} days." | tee -a main.log
        # check if all samples finished
        finished=$( cat finished.csv | sort -u | wc -l )
        started=$( cat tumors_and_normals.csv | grep -v "^#" | sort -u | wc -l )
        if [[ "$finished" -eq "$started" ]]; then
            # log and fetch MMR genes in annotations
            if [[ ${normal} == "PON" ]]; then
                echo "10: Fetching somatic variants of interest (${tumor}__${normal})" | tee -a main.log
            else
                echo "10: Fetching germline and somatic variants of interest (${tumor}__${normal})" | tee -a main.log
            fi

            # function to extrac rows that match MMR/tumor genes
            fetch_mmr_ann(){
              skip_rows=2
              col=10
              genes_screened="MLH1|MSH2|MSH6|PMS2|POLD1|POLE|IDH1|TP53|NF1|MCM8"
              if [[ $(echo $1 | grep "funcotator" &> /dev/null && echo 1) == 1 ]]; then
                skip_rows=3
                col=4
              fi
              awk -v FS="," -v SR=$skip_rows -v COL=${col} -v GENES=${genes_screened} \
              'BEGIN { print "# genes screened: "GENES}
              {
                if (NR >= SR){
                  if ($COL ~ /^(MLH1|MSH2|MSH6|PMS2|POLD1|POLE|IDH1|TP53|NF1|MCM8)$/){
                    print $0
                  }
                } else {
                 print $0
               }
              }' $1
              }
            export -f fetch_mmr_ann

            # get all MMR annotations into a CSV table
            ls analyses/annotations/*.annotations_funcotator_somatic.csv | parallel 'if [[ {#} == 1 ]]; then fetch_mmr_ann {}; else fetch_mmr_ann {} | tail -n +4; fi' > analyses/all_mmr_annotations_funcotator_somatic.csv
            ls analyses/annotations/*.annotations_snpeff_somatic.csv | parallel 'if [[ {#} == 1 ]]; then fetch_mmr_ann {}; else fetch_mmr_ann {} | tail -n +4; fi' > analyses/all_mmr_annotations_snpeff_somatic.csv
            ls analyses/annotations/*.annotations_funcotator_germline.csv | parallel 'if [[ {#} == 1 ]]; then fetch_mmr_ann {}; else fetch_mmr_ann {} | tail -n +4; fi' > analyses/all_mmr_annotations_funcotator_germline.csv
            ls analyses/annotations/*.annotations_snpeff_germline.csv | parallel 'if [[ {#} == 1 ]]; then fetch_mmr_ann {}; else fetch_mmr_ann {} | tail -n +3; fi' > analyses/all_mmr_annotations_snpeff_germline.csv

            # export to zip file
            today=$(date -I)
            zip -r export_results.${today}.zip annovar/*_multianno.maf \
                                               analyses/all_mmr_annotations_funcotator_somatic.csv \
                                               analyses/all_mmr_annotations_funcotator_germline.csv \
                                               analyses/all_mmr_annotations_snpeff_somatic.csv \
                                               analyses/all_mmr_annotations_snpeff_germline.csv \
                                               analyses/old_output*

            # add tmb_and_coverage to archive
            #zip -ru analyses.zip analyses/coverage_and_tmb.csv
            # tidyup and clean working dir
            if [[ ! -e bam ]]; then
                mv BQSR bam
            else
                mv BQSR/* bam
                rm -rf BQSR
            fi
            #
            # if [[ -e varscan/pileups ]]; then
            #   rm -rf varscan/pileups
            # fi
            # move stats to all_logfiles
            mv mutect2/*.filteringStats.tsv all_logfiles
            #mv mutect2/*.mutect2.unfiltered.wes.merged.vcf.stats all_logfiles
            # delete all other vcf files except filtered and unfiltered VCFs
            #rm $(ls mutect2/*vcf* | grep -v "filtered\|selected")
            # bgzip and tabix all vcf files
            ls mutect2/*.vcf | parallel --tmpdir ./.tmp "index-vcf {}"
            # delete directories with bam data
            rm -rf preprocessed_bam aligned_bam .tmp
            if [[ -e unmapped_bam ]]; then
                rm -rf unmapped_bam
            fi
            # write final logs and change permissions
            # final log
            echo "10: pipeline finished for batch." | tee -a main.log
            # log final
            date | tee -a main.log
            # get date pipeline started
            total_time_in_hours=$( how_long main.log h )
            total_time_in_days=$( how_long main.log d )
            # final log
            echo -e "\n10: pipeline took ${total_time_in_days} days or ${total_time_in_hours} hrs for batch to complete" | tee -a main.log
            # change/adjust permisions
            # this configuration allows the main user and the users in the tabori group to
            # read/write/excecute
            # dirs first
            echo "10: changing permissions" | tee -a main.log
            find . -type d -user `whoami` -exec chmod 775 {} \;
            find . -type f -user `whoami` -exec chmod 664 {} \;
            echo "@@@@ All done. See you next time! @@@@" | tee -a main.log
        fi
    fi
    # run mutational signature analysis
    # Rscript ${pipeline_dir}/cosmic_signature_analysis.R ${mode} ${tumor}__${normal}.

    # last move logfile to dir
    if [[ -e ${tumor}__${normal}.analyses.${tissue}.log ]]; then
        mv ${tumor}__${normal}.analyses.${tissue}.log all_logfiles
    fi
fi
