#!/bin/bash

vcf_file=$1
bcftools view $vcf_file | \
perl -ne '
  if ($_ =~ /^#/){
    if ($_ =~ m/##Funcotator /){
      print $_;
    }
    if ($_ =~ m/##INFO=<ID=FUNCOTATION/){
      $_ =~ m/Funcotation fields are: (.+?)\"/;
      @heads = split /\|/, $1;
      map { $_ =~ s/Gencode_\d+_//; substr($_, 0, 1) = uc substr($_, 0, 1) } @heads;
      $heads[0] = "Hugo_Symbol";
      print join(",", @heads) . "\n";
    }
  } else {
    @fields = split /\t/, $_;
    $fields[7] =~ m/FUNCOTATION=\[(.+?)\]/;
    @maf = split /\|/, $1;
    print join(",", @maf) . "\n";
  }
' | sed "3,$ s/^/${2},${3},${4},/; 2,2 s/^/Tumor,Normal,Source,/1p' | sed '2d'
