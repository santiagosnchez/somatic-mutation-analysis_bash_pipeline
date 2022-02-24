#!/bin/bash

vcf_file=$1
bcftools view $vcf_file | \
perl -ne '
  if ($_ =~ /^#/){
    if ($_ =~ m/##Funcotator /){
      print $_;
    }
    if ($_ =~ m/##INFO=<ID=ANN/){
      $_ =~ m/Functional annotations: \'(.+?)\' \">/;
      @heads = split /\ | /, $1;
      print join("\t", @heads) . "\n";
    }
  } else {
    @fields = split /\t/, $_;
    $fields[7] =~ m/;*ANN=(.+?);*/;
    @maf = split /\|/, $1;
    print join("\t", @maf) . "\n";
  }
'
