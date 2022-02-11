#!/bin/perl

# Description:
# Parses two gzipped fastq files and checks if
# fwd and rev reads are properly [paired
# and saves singletons to separate file

# record start time
my $start = time();

# check args
if (scalar(@ARGV) != 2){
   die "Needs forward and reverse read files (only two gzipped files)\n";
}

# location of fastq files
my $fastq_fwd = $ARGV[0];
my $fastq_rev = $ARGV[1];
# split path
my @path1 = split /\//, $fastq_fwd;
my @path2 = split /\//, $fastq_rev;
# shorter base names for file
my $fastq_SNF;
my $fastq_SNR;

# find the last element in path and replace suffix
if (substr($fastq_fwd, 0, 1) eq "/"){
    $path1[ $#path1 ] =~ s/_R1.*//;
    $fastq_SNF = $path1[ $#path1 ];
} else {
    $path1[ $#path1 - 1 ] =~ s/_R1.*//;
    $fastq_SNF = $path1[ $#path1 - 1 ];
}
if (substr($fastq_rev, 0, 1) eq "/"){
    $path2[ $#path2 ] =~ s/_R2.*//;
    $fastq_SNR = $path2[ $#path2 ];
} else {
    $path2[ $#path2 - 1 ] =~ s/_R2.*//;
    $fastq_SNR = $path2[ $#path2 -1 ];
}

# initialize filehandles
my $IN_FWD;
my $IN_REV;
my $OUT_FWD;
my $OUT_REV;
my $OUT_SIN;

# open gzippped files (input and output)
open(IN_FWD, "-|", "/usr/bin/gunzip", "-c", $fastq_fwd);
open(IN_REV, "-|", "/usr/bin/gunzip", "-c", $fastq_rev);
open(OUT_FWD, "|-", "/usr/bin/gzip > tmp/$fastq_SNF.1.fastq.gz");
open(OUT_REV, "|-", "/usr/bin/gzip > tmp/$fastq_SNR.2.fastq.gz");
open(OUT_SIN, "|-", "/usr/bin/gzip > tmp/$fastq_SNF.S.fastq.gz");

# start fwd and rev dictionaries/hashes
# and other variables
my %fwd={};
my %rev={};
my $block1 = "";
my $block2 = "";
my $head1;
my $head2;
my $total_paried=0;
my $total_sing=0;

# parse both files simultaneously
while ($line_fwd = <IN_FWD>, $line_rev = <IN_REV>){
    # if first line is header
    if (($. % 4) == 1){
        $block1 .= $line_fwd;
        $block2 .= $line_rev;
        chomp($line_fwd);
        chomp($line_rev);
        $head1 = $line_fwd;
        $head2 = $line_rev;
        $head1 =~ s/ .*|\/[12]$//g;
        $head2 =~ s/ .*|\/[12]$//g;
    }
    # if line is last for block (i.e., quality scores)
    elsif (($. % 4) == 0){
        # complete block
        $block1 .= $line_fwd;
        $block2 .= $line_rev;
        # if headers match write to file
        if ($head1 eq $head2){
           print OUT_FWD $block1;
           print OUT_REV $block2;
           # count records
           $total_paired += 1;
        }
        else {
           # look for fwd header in rev dictionary
           # save to file if present and delete entry
	   if (exists $rev{$head1}){
               print OUT_FWD $block1;
               print OUT_REV $rev{$head1};
               # saves memory
               delete $rev{$head1};
               # count records
               $total_paired += 1;
           }
           else {
               # otherwise add to dictionary
               $fwd{$head1} = $block1;
           }
           # same as previous block, but for reverse
           if (exists $fwd{$head2}){
               print OUT_FWD $fwd{$head2};
               print OUT_REV $block2;
               delete $fwd{$head2};
               $total_paired += 1;
           }
           else {
               $rev{$head2} = $block2;
           }
        }
        # reset read blocks for fwd and rev
        $block1 = "";
        $block2 = "";
    }
    # if line is not header or last
    # append to block
    else {
        $block1 .= $line_fwd;
        $block2 .= $line_rev;
    }
}
# close open filehandles
close(OUT_FWD);
close(OUT_REV);
close(IN_FWD);
close(IN_REV);

# write unmatched reads to file
# starting with fwd
if (scalar(keys %fwd) > 1){
    foreach $head (keys %fwd){
        print OUT_SIN $fwd{$head};
        delete $fwd{$head};
        # count singletons
        $total_sing += 1;
    }
}

# then for rev
if (scalar(keys %rev) > 1){
    foreach $head (keys %rev){
        print OUT_SIN $rev{$head};
        delete $rev{$head};
        # count singletons
        $total_sing += 1;
    }
}

# final time
my $end = time();
my $total_hours = ($end - $start)/60/60;

# log to screen
print "Paried: $total_paired\n";
print "Single: $total_sing\n";
print "Files written to:
tmp/$fastq_SNF.1.fastq.gz
tmp/$fastq_SNR.2.fastq.gz
tmp/$fastq_SNF.S.fastq.gz\n";
printf "%.5f hours to finish\n", $total_hours;
