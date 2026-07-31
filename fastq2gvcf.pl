#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Std;
use File::Basename qw/ fileparse /; 

# Maps reads in FASTQ file to indexed PHG pangenome graph in order to imputate and call
# haplotypes across graph ranges, producing a PHG-imputed FASTA sequence which is ultimately
# mapped back to the graph reference (MorexV3 in barley) to call variants a produce a gVCF file
#
# J Sarria, B Contreras-Moreira
# Copyright [2026] Estacion Experimental de Aula Dei-CSIC

my ($cmd,$root,$key,$val,%opts,%config);
my ($fqfile,$redo) = ('', 0);
my ($threads,$chunksize,$minmatch,$outdir) = (5, 500, 101, '/tmp');

getopts('hRf:c:o:m:t:k:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)) {
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-f input FASTQ file                                (example: -f sample.fastq.gz)\n";
  print "-c graph config file                               (example: -c /path/to/graph/config.impute.yaml)\n";
  print "-o output folder                                   (optional, example: -o mysample, default -o $outdir)\n";
  print "-m min match length                                (optional, example: -m 150, default -m $minmatch)\n";
  print "-t threads                                         (optional, example: -t 12, default -d $threads)\n";
  print "-k chunk size                                      (optional, example: -k 1000, default -k $chunksize\n";
  print "-R redo all steps even if results are in place     (optional, by default previous results are re-used)\n";
  #print "\nPrimary citation: https://www.biorxiv.org/content/10.1101/2025.07.17.665301v1\n";
  exit(0);
}

if(!defined($opts{'f'})) {
  die "# ERROR: need input FASTQ file (-f)\n";
} else {
  $fqfile = $opts{'f'};
  $root = fileparse($fqfile, qr/\..*/);
}

if(!defined($opts{'c'})) {
  die "# ERROR: need graph config file (-c)\n";
} else {
  open(CONFIG,"<",$opts{'c'}) ||
    die "# ERROR: cannot read graph config file (-c)\n";
  while(<CONFIG>) {
    next if(/^#/);
    if(/([^:]+):\s+(\S+)/) {
      ($key,$val) = ($1, $2); 
      if(!-e $val) {
        die "# ERROR: bad value for '$key' in $opts{'c'}\n";
      } elsif($key eq 'kmer-index-file') {
        if(!-e $val.'.ssa') {
        die "# ERROR: need also kmer index file $val\.ssa, check $opts{'c'} ($key)\n";
        } elsif(!-e $val.'.len.gz') {
        die "# ERROR: need also kmer index file $val\.len.gz, check $opts{'c'} ($key)\n";
        } 
      } elsif($key eq 'hvcf-dir') {
        opendir(HVCF,$val) ||
          die "# ERROR: cannot list $val\n";
        my @files = grep{/\.h.vcf.gz/} readdir(HVCF);
        closedir(HVCF);

        if(scalar(@files) == 0) {
          die "# ERROR: cannot find h.vcf.gz files in $val\n";
        }
      } 

      $config{$key} = $val	
    }
  } 
  close(CONFIG);
} 

if(defined($opts{'o'})) {
  $outdir = $opts{'o'};

  if(-d $opts{'o'}) {
    warn "# WARNING : folder '$outdir' exists, files might be overwritten\n";
  } else {
    if ( !mkdir($outdir) ) {
      die "# ERROR: cannot create $outdir\n";
    }	  
  }
}

if(defined($opts{'m'}) && $opts{'m'} >= 0) {
  $minmatch = $opts{'m'}
}

if(defined($opts{'t'}) && $opts{'t'} >= 0) {
  $threads = $opts{'t'}
}

if(defined($opts{'k'}) && $opts{'k'} >= 0) {
  $chunksize = $opts{'k'}
}

if(defined($opts{'R'})) {
  $redo = 1
}

warn "## $0 -f $fqfile -c $opts{'c'} -o $outdir -m $minmatch -t $threads -k $chunksize -R $redo\n\n";

# 0) check phg is loaded (conda)
if(!`which phg`) {
  die "# ERROR: cannot find phg in \$PATH, perhaps you need to conda activate it?\n";  
}

# 1) map reads
my $mapfile = "$outdir/$root" . '_1_readMapping.txt';
if($redo == 1 || !-e $mapfile) {
  $cmd = "phg map-reads --index $config{'kmer-index-file'} --read-files $fqfile -o $outdir --hvcf-dir $config{'hvcf-dir'} " .
    "--threads $threads --min-mem-length $minmatch "; #--conda-env-prefix $config{'conda-env-prefix'}";
  run_cmd("# 1 Running phg map-reads ...", $cmd);

} else {
  print "# 1 re-using $mapfile\n";
}	

# 2) find paths
my $hvcf_file = "$outdir/$root" . '_1.h.vcf';
if($redo == 1 || !-e $hvcf_file) {
  $cmd = "phg find-paths --read-files $mapfile --output-dir $outdir --hvcf-dir $config{'hvcf-dir'} --path-type haploid " .
    "--threads $threads --reference-genome $config{'reference-fasta'}";
  run_cmd("# 2 Running phg find-paths ...", $cmd);

} else {
  print "# 2 re-using $hvcf_file\n";
}

# 3) create fasta from hvcf
my $composite_file = "$outdir/$root" . '_1_composite.fa';
if($redo == 1 || !-e $composite_file) {
  $cmd = "phg create-fasta-from-hvcf --hvcf-file $hvcf_file -o $outdir --db-path $config{'db-path'} " .
    "--range-bedfile $config{'range-bed'} --fasta-type composite"; #--conda-env-prefix $config{'conda-env-prefix'}";
  run_cmd("# 3 Running phg create-fasta-from-hvcf ...", $cmd);
  
} else {
  print "# 3 re-using $composite_file\n";
}

# 4) index fasta
my $composite_index_file = $composite_file . '.fai';
if($redo == 1 || !-e $composite_index_file) {
  $cmd = "samtools faidx $composite_file";
  run_cmd("# 4 Indexing composite fasta ...", $cmd);

} else {
  print "# 4 re-using $composite_index_file\n";
}



sub run_cmd {
  my ($message, $cmd) = @_;

  warn $message;
  system($cmd);
  if($? != 0) {
      die "# EXIT: command failed ($cmd)\n";
  }
}


