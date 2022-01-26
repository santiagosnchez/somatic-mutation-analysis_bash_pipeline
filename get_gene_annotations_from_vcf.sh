#!/bin/bash

args=($@)
vcf_file=${args[0]}
sample=$( echo $vcf_file | sed 's/\..*//')
sample=$( echo $sample | rev | cut -d'/' -f1 | rev)

for i in `seq 2 ${#args[@]}`; do
    i=$(( i - 1 ))
    gene=${args[i]}
    if [[ "${vcf_file}" = *".gz" ]]; then
        bcftools view -H $vcf_file | \
        grep "|${gene}|" | grep -v "synonymous" | \
        awk -v OFS="," -v g=$gene -v s=$sample '{
          gt = $NF;
          gsub(":.*","",gt);
          gsub("__",",",s);
          split($8, x, ";");
          split(x[length(x)], y, "=");
          split(y[2], z, ",");

          for (i=1; i <= length(z); ++i){
              split(z[i], a, "|");
              if (a[4] == g){
                  print s,$1,$2,$4,$5,gt,a[2],a[3],a[4],a[6],a[7],a[10],a[11];
              }
          }
        }'
    else
        bcftools view -H $vcf_file | \
        grep "|${gene}|" | \
        awk -v OFS="," -v g=$gene -v s=$sample '{
            gt = $NF;
            gsub(":.*","",gt);
            gsub("__",",",s);
            split($8, x, ";");
            split(x[length(x)], y, "=");
            split(y[2], z, ",");

            for (i=1; i <= length(z); ++i){
                split(z[i], a, "|");
                if (a[2] ~ /missense_variant/ && a[4] == g){
                    print s,$1,$2,$4,$5,gt,a[2],a[3],a[4],a[6],a[7],a[10],a[11];
                }
            }
        }'
    fi
done
