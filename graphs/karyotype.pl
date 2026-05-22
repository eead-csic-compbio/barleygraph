#!/usr/bin/env perl
use strict;
use warnings;

# see https://circos.ca/documentation/tutorials/ideograms/karyotypes

if(!$ARGV[0]){ die "# usage: $0 hapIDranges.tsv\n" }

my ($chr, $end, $color, %coords, %ranges, @order);

open(TSV,"<",$ARGV[0]) ||
  die "# cannot read $ARGV[0]\n";
while(<TSV>) {
  #chr1H 1 76743 ...
  if(/^(\S+)\t\d+\t(\d+)\t/) {
    ($chr,$end) = ($1, $2);

    if(!defined($coords{$chr})) {
      push(@order, $chr);
    }

    $coords{$chr} = $end;
    $ranges{$chr}++;
  }
}
close(TSV);

foreach $chr (@order) {

  $color = $chr;
  $color =~ s/H//;

  printf("chr\t-\t%s\t%s\_ranges=%d\t0\t%d\t%s\n",
    $chr,$chr,$ranges{$chr},$coords{$chr},$color);	  
}
