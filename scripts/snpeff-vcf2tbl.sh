#!/bin/bash

vcf_file=$1
bcftools view $vcf_file | \
perl -ne '
  if ($_ =~ /^#/){
    if ($_ =~ m/##INFO=<ID=ANN/){
      $_ =~ m/Functional annotations: (.+?)\">/;
      @heads = split / \| /, $1;
      $str_head = join(",", @heads);
      print "Chromosome,Position,Reference," . substr($str_head, 1, length($str_head)-2) . ",TumorAlleleFrequency" . "\n";
    }
  } else {
    @fields = split /\t/, $_;
    @annfields = split /;/, $fields[7];
    %annfieldsmap = {};
    foreach (@annfields){
      @lab_ann = split /=/, $_;
      $annfieldsmap{$lab_ann[0]} = $lab_ann[1];
    }
    @fmt = split /:/, $fields[8];
    @tgt = split /:/, $fields[10];
    %fmt_gt_map = {};
    foreach (0 .. (scalar(@fmt)-1)){ $fmt_gt_map{$fmt[$_]} = $tgt[$_] };
    @annsubfields = split(/,/, $annfieldsmap{"ANN"});
    foreach (@annsubfields){
       $_ =~ s/\|/,/g;
       print join(",", @fields[0,1,3]) . "," . $_ . "," . $fmt_gt_map{"AF"} . "\n";
    }

  }
' | sed "2,$ s/^/${2},${3},${4},/; 1 s/^/Tumor,Normal,Source,/1p" | sed '1d'
