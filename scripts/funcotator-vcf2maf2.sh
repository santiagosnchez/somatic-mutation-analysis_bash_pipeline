#!/bin/bash

vcf_file=$1
export tumor=$2
export normal=$3
export tissue=$4

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
      $new_head = "Tumor|Normal|Source|" . join("|", @heads) . "|TumorAlleleFrequency\n";
      print $new_head;
    }
  } else {
    @fields = split /\t/, $_;
    $fields[7] =~ m/FUNCOTATION=\[(.+?)\]/;
    $maf = $1;
    @fmt = split /:/, $fields[8];
    @tgt = split /:/, $fields[10];
    %fmt_gt_map = {};
    foreach (0 .. (scalar(@fmt)-1)){ $fmt_gt_map{$fmt[$_]} = $tgt[$_] };
    $af = $fmt_gt_map{"AF"};
    $outline = $ENV{tumor} . "|" . $ENV{normal} . "|" . $ENV{tissue} . "|" . $maf . "|" . $af . "\n";
    print $outline;
  }
' | tr '|' ','
