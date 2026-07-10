#!/usr/bin/env perl
use strict;
use warnings;

my $MINDEPTH = 5;

if(!$ARGV[3]) {
  die "# usage: $0 <input VCF> <input sample file> <out dist file> <out count file>\n"
}

my ($VCFfile, $sampleFile, $distFile, $countFile) = @ARGV;

if(!-e $VCFfile) {
  die "# ERROR: need a valid .vcf.gz file\n";
} elsif(!-e $sampleFile) {
  die "# ERROR: need a valid text file with samples in .vcf.gz file\n";
}

open(DIST,">$distFile") || 
  die "# ERROR: need a valid dist matrix file path\n";
open(COUNT,">$countFile") || 
  die "# ERROR: need a valid count matrix file path\n";

####################################################

my ($b1,$b2,$cmd);
my ($al11,$al12,$d1,$al21,$al22,$d2);

my (%VCFcolumns, @barleys);
my $coln = 10;
open(SAMPLES,"<",$sampleFile);
while(<SAMPLES>) {
  if(/(\S+)/) {
    $VCFcolumns{$1} = $coln;
    $coln++;
    push(@barleys,$1);
  }
}  
close(SAMPLES);

foreach $b1 (@barleys) {

  print DIST "$b1";
  print COUNT "$b1";

    foreach $b2 (@barleys) {
      if($b1 eq $b2) { 
	print DIST "\t0.0000";
        print COUNT "\t0.0000";
      }
      else {

        my ($eq,$neq) = (0,0); #print "$VCFcolumns{$b1},$VCFcolumns{$b2}\n";
        open(VCF,"zcat $VCFfile | cut -f $VCFcolumns{$b1},$VCFcolumns{$b2} | ") || 
          die "# ERROR: cannot parse zcat $VCFfile | cut -f $VCFcolumns{$b1},$VCFcolumns{$b2}\n";
        while(<VCF>) {
          if(/([^\/])\/([^\/]):(\S+)\t([^\/])\/([^\/]):(\S+)/){
            ($al11,$al12,$d1,$al21,$al22,$d2) = ($1,$2,$3,$4,$5,$6);
	        
            next if($al11 eq '.' && $al12 eq '.'); # missing
            next if($al21 eq '.' && $al22 eq '.');

            next if($al11 ne $al12); # heterozygous
            next if($al21 ne $al22);

            next if($d1 < $MINDEPTH || $d2 < $MINDEPTH); # low depth

	    if($al11 eq $al21){ $eq++ } 
   	    else { $neq++}
          }
        }
        close(VCF);
				
        printf(DIST "\t%1.4f",1-($eq/($eq+$neq))); 
	printf(COUNT "\t%d",($eq+$neq));
      }
    }
    print DIST "\n";
    print COUNT "\n";
}

