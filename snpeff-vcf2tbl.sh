#!/bin/bash

vcf_file=$1
bcftools view $vcf_file | \
perl -ne '
  if ($_ =~ /^#/){
    if ($_ =~ m/##INFO=<ID=ANN/){
      $_ =~ m/Functional annotations: (.+?)\">/;
      @heads = split / \| /, $1;
      $str_head = join(",", @heads);
      print "Chromosome,Position,Reference," . substr($str_head, 1, length($str_head)-2) . "\n";
    }
  } else {
    @fields = split /\t/, $_;
    @annfields = split /;/, $fields[7];
    %annfieldsmap = {};
    foreach (@annfields){
      @lab_ann = split /=/, $_;
      $annfieldsmap{$lab_ann[0]} = $lab_ann[1];
    }
    @annsubfields = split(/,/, $annfieldsmap{"ANN"});
    foreach (@annsubfields){
       $_ =~ s/\|/,/g;
       print join(",", @fields[0,1,3]) . "," . $_ . "\n";
    }
  }
'
