#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Std;
use File::Basename qw/ fileparse /; 

# Map reads in FASTQ file(s) on indexed PHG pangenome graph in order to imputate and call
# haplotypes across graph ranges. By default produces hVCF output.
# Uses the following Linux system binaries: find, sort
#
# J Sarria, B Contreras-Moreira
# Copyright [2026] Estacion Experimental de Aula Dei-CSIC
#
#
# (left out due to space constraints, see code commented below)
# Optionally can produce a PHG-composite FASTA sequence which is mapped back to the graph reference
# to call variants finally output in a gVCF file.
# NOTE: heterozygous sites are explicitely excluded from output VCF as might be artificially 
# introduced by chunks aligning to same reference regions.

my ($cmd,$root,$key,$val,$cfile,%opts,%config,@temp);
my ($fqfiles,$output_file,$dogVCF,$redo,$allsites) = ('', '', 0, 0, 0);
my ($threads,$chunksize,$minmatch,$outdir) = (5, 500, 101, '/tmp');
my $graph_db = "/graph_db";
my $graph_list_file = $graph_db . '/graph_list.tsv';

my $agcEXE      = 'agc';
#my $samtoolsEXE = 'samtools';
#my $bwaEXE      = 'minibwa';
#my $bedtoolsEXE = 'bedtools';
#my $bcftoolsEXE = 'bcftools';

getopts('hlARgB:1:2:G:c:o:m:t:k:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)) {
  print "Maps reads in FASTQ file to imputate and call haplotypes from graph.\n\n";
  print "Usage: $0 [options]\n\n";
  print "-h this message\n";
  print "-1 input FASTQ file                                (example: -1 sample1.fq.gz)\n";
  print "-2 paired FASTQ file                               (optional, example: -2 sample2.fq.gz)\n";
  print "-l list supported graphs                           (optional, checks $graph_list_file)\n";
  print "-G graph name                                      (-G or -c, example: -G Pan20-mmap-pro)\n";
  print "-c graph config file                               (-G or -c, example: -c /path/graph/config.impute.yaml)\n";
  #print "-g produce gVCF file                               (optional, requires vcf_dbs in config)\n";
  print "-o output folder                                   (optional, example: -o mysample, default -o $outdir)\n";
  print "-m min match length                                (optional, example: -m 150, default -m $minmatch)\n";
  print "-t threads                                         (optional, example: -t 12, default -d $threads)\n";
  #print "-k chunk size in bases                             (optional, example: -k 1000, default -k $chunksize\n";
  #print "-B path to [mini]bwa binary                        (optional, example: -B /path/to/[mini]bwa)\n";
  #print "-A all sites in output VCF, not just variants      (optional, by default only variants are considered)\n";
  print "-R redo all steps even if results are in place     (optional, by default previous results are re-used)\n";
  #print "\nPrimary citation:\n";
  exit(0);
}

if(defined($opts{'l'})) {
  open(LIST,"<",$graph_list_file) ||
      die "# ERROR: no supported graphs ($graph_list_file), run setup_graph first\n";
  while(<LIST>) {
    @temp = split;	  
    print "$temp[1]\n";
  }   
  close(LIST);
  exit(0);
}

if(!defined($opts{'1'})) {
  die "# ERROR: need input FASTQ file (-1)\n";
} else {
  $fqfiles = $opts{'1'};
  $root = fileparse($fqfiles, qr/\..*/);
}

if(defined($opts{'2'})) {
  $fqfiles .= ",$opts{'2'}";
}

if(defined($opts{'G'})) {
  
  open(LIST,"<",$graph_list_file) ||
      die "# ERROR: no supported graphs ($graph_list_file), run setup_graph first\n";
  while(<LIST>) {
    @temp = split;
    if($temp[1] eq $opts{'G'}) {
      $cfile = "$graph_db/$temp[0]/$temp[1]/$temp[1].yaml";
      last;
    }
  }
  close(LIST);

  if($cfile eq '') {
    die "# ERROR: cannot find $opts{'G'} in ($graph_list_file)\n";
  } 
} elsif(defined($opts{'c'})) {
  $cfile = $opts{'c'};
} else {
  die "# ERROR: need either graph name (-G) or config file (-c)\n";
}

# actually read config file
open(CONFIG,"<",$cfile) ||
  die "# ERROR: cannot read config file $cfile, try -l to list supported graphs\n";
while(<CONFIG>) {
  next if(/^#/);
  if(/([^:]+):\s+(\S+)/) {
    ($key,$val) = ($1, $2); 
    if($key ne "reference_name" && $key ne "reference_fasta" &&	!-e $val) {
      die "# ERROR: bad value for '$key' in $opts{'c'}\n";
    } elsif($key eq 'kmer_index') {
      if(!-e $val.'.ssa') {
        die "# ERROR: need also kmer index file $val\.ssa, check $opts{'c'} ($key)\n";
      } elsif(!-e $val.'.len.gz') {
        die "# ERROR: need also kmer index file $val\.len.gz, check $opts{'c'} ($key)\n";
      } 
    } elsif($key eq 'hvcf_bed') {
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

if(defined($opts{'g'})) {
  $dogVCF = 1
}

if(defined($opts{'m'}) && $opts{'m'} >= 0) {
  $minmatch = $opts{'m'}
}

if(defined($opts{'t'}) && $opts{'t'} >= 0) {
  $threads = $opts{'t'}
}

#if(defined($opts{'k'}) && $opts{'k'} >= 0) {
#  $chunksize = $opts{'k'}
#}
#if(defined($opts{'B'})) {
#  $bwaEXE = $opts{'B'}
#}
#if(defined($opts{'A'})) {
#  $allsites = 1
#}

if(defined($opts{'R'})) {
  $redo = 1
}

$output_file = "$outdir/$root" . ".$chunksize.vcf.gz";
if($dogVCF == 0) {
  $output_file = "$outdir/$root" . '_1.h.vcf';
}

print "## read file(s): $fqfiles\n";
print "## config file: $cfile\n";
print "## params: -o $outdir -m $minmatch -t $threads -R $redo\n\n";

# 0.1) check output
if($redo == 0 && -e $output_file) {
  print "## output: $output_file\n";
  exit(0);
}

# 0.2) check phg is loaded (conda)
if(!`which phg`) {
  die "# ERROR: cannot find phg in \$PATH, perhaps you need to conda activate it?\n";  
}

# 1) map reads
my $mapfile = "$outdir/$root" . '_1_readMapping.txt';
my $hvcf_file = "$outdir/$root" . '_1.h.vcf';
push(@temp, $mapfile);
if($redo == 1 || !-e $hvcf_file) {
  $cmd = "phg map-reads --index $config{'kmer_index'} --read-files $fqfiles -o $outdir " .
    "--hvcf_bed $config{'hvcf_bed'} --threads $threads --min-mem-length $minmatch "; 
  run_cmd($cmd, "# 1 Running phg map-reads ...");

} else {
  print "# 1 re-using $mapfile\n";
}	

# 2) find paths
if($redo == 1 || !-e $hvcf_file) {

  # 2.1) make sure reference FASTA exists
  if(!-e $config{'reference_fasta'}) {
    my $rname = `agc listref $config{'agc_assemblies'}`;
    chomp $rname; 
    $cmd = "agc getset $config{'agc_assemblies'} $rname > $config{'reference_fasta'}";
    run_cmd($cmd, "# 2.1 Extracting reference sequence (1st time only)...");

  } elsif(!-e $config{'agc_assemblies'}) {
    die "# ERROR: cannot find $config{'agc_assemblies'}, please fix $cfile\n";
  }	  

  $cmd = "phg find-paths --read-files $mapfile --output-dir $outdir --hvcf_bed $config{'hvcf_bed'} " .
    "--path-type haploid --threads $threads --reference-genome $config{'reference_fasta'}";
  run_cmd($cmd, "# 2 Running phg find-paths ...");

} else {
  print "# 2 re-using $hvcf_file\n";
}

if($dogVCF == 0) {
  print "## output: $output_file\n";
  unlink(@temp);
  exit(0);
}

# 3) create fasta from hvcf
#my $composite_file = "$outdir/$root" . '_1_composite.fa';
#push(@temp, $composite_file);
#if($redo == 1 || !-e $composite_file) {
#  $cmd = "phg create-fasta-from-hvcf --hvcf-file $hvcf_file -o $outdir --db-path $config{'db_path'} " .
#    "--range_bedfile $config{'range_bed'} --fasta-type composite";
#  run_cmd($cmd, "# 3 Running phg create-fasta-from-hvcf ...");  
#} else {
#  print "# 3 re-using $composite_file\n";
#}

# 4) index fasta
#my $composite_index_file = $composite_file . '.fai';
#push(@temp, $composite_index_file);
#if($redo == 1 || !-e $composite_index_file) {
#  $cmd = "$samtoolsEXE faidx $composite_file";
#  run_cmd($cmd, "# 4 Indexing composite fasta ...");
#} else {
#  print "# 4 re-using $composite_index_file\n";
#}

# 5) chunk composite sequence
#my $composite_chunks_bam = "$outdir/$root" . "_1_composite.$chunksize.bam";
#push(@temp, $composite_chunks_bam, "$composite_chunks_bam.csi");
#if($redo == 1 || !-e $composite_chunks_bam) {
#  my $comp_genome_bed  = "$outdir/$root" . '_1_comp_genome.bed';
#  my $comp_chunks_bed  = "$outdir/$root" . '_1_comp_chunks.bed';
#  my $comp_chunks_file = "$outdir/$root" . '_1_comp_chunks.fa';
#  my $ref_index_file   = $config{'reference_fasta'} . '.mbw';
#  my @temp2 = ($comp_genome_bed, $comp_chunks_bed, $comp_chunks_file);
#
#  $cmd = "perl -lane 'print \"\$F[0]\t0\t\$F[1]\"' $composite_index_file > $comp_genome_bed";
#  run_cmd($cmd) if(!-e $comp_genome_bed || $redo == 1);
#  
#  $cmd = "$bedtoolsEXE makewindows -b $comp_genome_bed -w $chunksize > $comp_chunks_bed";
#  run_cmd($cmd) if(!-e $comp_chunks_bed || $redo == 1);
#  
#  $cmd = "$bedtoolsEXE getfasta -fi $composite_file -bed $comp_chunks_bed -fo $comp_chunks_file";
#  run_cmd($cmd, "# 5.1 Chunking composite fasta ...") if(!-e $comp_chunks_file || $redo == 1);
#
#  $cmd = "$bwaEXE index -t $threads $config{'reference_fasta'}";
#  if($bwaEXE !~ /minibwa/) {
#    $ref_index_file   = $config{'reference_fasta'} . '.sa';
#    $cmd = "$bwaEXE index $config{'reference_fasta'}";
#  }
#  run_cmd($cmd, "# 5.2 Indexing reference fasta ...") if(!-e $ref_index_file || $redo == 1);
#
#  $cmd = "$bwaEXE map -t $threads $config{'reference_fasta'} $comp_chunks_file | " .
#    "$samtoolsEXE sort -@ $threads -o $composite_chunks_bam && $samtoolsEXE index -c -@ $threads $composite_chunks_bam";
#  if($bwaEXE !~ /minibwa/) {
#    $cmd = "$bwaEXE mem -t $threads $config{'reference_fasta'} $comp_chunks_file | " .
#    "$samtoolsEXE sort -@ $threads -o $composite_chunks_bam && $samtoolsEXE index -c -@ $threads $composite_chunks_bam";
#  }
#  run_cmd($cmd, "# 5.3 Mapping composite fasta ...");
#  unlink(@temp2);
#  
#} else {
#  print "# 5 re-using $composite_chunks_bam\n";
#}

# 6) variant call, produces final gVCF output file with het sites filtered out
#if($redo == 1 || !-e $output_file) {
#  $cmd = "$bcftoolsEXE mpileup -Ou -f $config{'reference_fasta'} $composite_chunks_bam | ";
#  if($allsites == 1) {
#    $cmd .= "$bcftoolsEXE call -m --ploidy 2 -Oz | bcftools view -g ^het -o $output_file";
#  } else {
#    $cmd .= "$bcftoolsEXE call -m --ploidy 2 -v -Oz | bcftools view -g ^het -o $output_file";
#  }
#  run_cmd($cmd, "# 6 Variant calling ...");
#  print "## output: $output_file\n";
#} else {
#  print "## output $output_file (re-used)\n";
#}

# clean only if successful
#unlink(@temp);
#exit(0);


sub run_cmd {
  my ($cmd, $message) = @_;

  print "$message\n" if($message);
  system($cmd);
  if($? != 0) {
      die "# EXIT: command failed ($cmd)\n";
  }
}
