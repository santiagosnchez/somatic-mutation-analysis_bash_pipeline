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
      $num_fields = scalar(@heads);
      push @heads, "TumorAlleleFrequency";
      print join(",", @heads) . "\n";
    }
  } else {
    @fields = split /\t/, $_;
    $fields[7] =~ m/FUNCOTATION=\[(.+?)\]/;
    @maf = split /\|/, $1;
    @fmt = split /:/, $fields[8];
    @tgt = split /:/, $fields[10];
    %fmt_gt_map = {};
    foreach (0 .. (scalar(@fmt)-1)){ $fmt_gt_map{$fmt[$_]} = $tgt[$_] };
    $num_maf_fields = scalar(@maf);
    if ($num_maf_fields == $num_fields){
      print join(",", @maf) . "," . $fmt_gt_map{"AF"} . "\n";
    } else {
      $add_extra_commas = $num_fields - $num_maf_fields;
      print join(",", @maf) . ("," x $add_extra_commas) . "," . $fmt_gt_map{"AF"} . "\n";
    }
  }
' | sed "3,$ s/^/${2},${3},${4},/; 2,2 s/^/Tumor,Normal,Source,/1p" | sed '2d'
