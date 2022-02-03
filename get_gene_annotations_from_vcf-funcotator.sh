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
    chomp($_);
    @fields = split /\t/, $_;
    $fields[7] =~ m/FUNCOTATION=\[(.+?)\]/;
    @sample = split "__", $ENV{sample};
    @maf = split /\|/, $1;
    @vcf = @sample;
    push @vcf, @fields[0,1,3,4,10];
    if ($fields[8] =~ m/:FREQ:*/){
      @format = split /:/, $fields[8];
      @gt = split /:/, $fields[10];
      @idx = grep { $format[$_] eq "FREQ" } 0 .. $#format;
      $vcf[6] = $gt[0];
      $gt[$idx] =~ s/%//;
      push @vcf, ($gt[$idx[0]]/100);
    }
    elsif ($fields[8] =~ m/:AF:*/){
      @format = split /:/, $fields[8];
      @gt = split /:/, $fields[10];
      @idx = grep { $format[$_] eq "AF" } 0 .. ($#format-1);
      $vcf[6] = $gt[0];
      push @vcf, ($gt[$idx[0]]);
    }
    else {
      $vcf[6] =~ s/:.+//;
      push @vcf, ("NA");
    }
    print join(",", @vcf) . "," . join(",", @maf[0,5,7,12,13,15,16,18]) . "\n";
    ' | \
    grep ",${gene},"
done
