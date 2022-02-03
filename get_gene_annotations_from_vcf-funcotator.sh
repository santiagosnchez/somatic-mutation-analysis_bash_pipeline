#!/bin/bash

args=($@)
vcf_file=${args[0]}
sample=$( echo $vcf_file | sed 's/\..*//')
sample=$( echo $sample | rev | cut -d'/' -f1 | rev)
export sample

for i in `seq 2 ${#args[@]}`; do
    i=$(( i - 1 ))
    gene=${args[i]}
    bcftools view -H $vcf_file | \
    perl -ne '
    @fields = split /\t/, $_;
    $fields[7] =~ m/FUNCOTATION=\[(.+?)\]/;
    @sample = split "__", $ENV{sample};
    @maf = split /\|/, $1;
    @data = @sample;
    push @data, @fields[0,1,3,4,9];
    $data[6] =~ s/:.*//;
    push @data, @maf[0,5,7,12,13,15,16,18];
    print join(",", @data) . "\n";' | \
    grep ",${gene},"
done
