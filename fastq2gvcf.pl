#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Std;
use File::Basename qw/ fileparse /; 

# Maps reads in FASTQ file to indexed PHG pangenome graph in order to imputate and call
# haplotypes across graph ranges, producing a PHG-imputed FASTA sequence which is ultimately
# mapped back to the graph reference (MorexV3 in barley) to call variants a produce a gVCF file.
#
# NOTE: heterozygous sites are explicitely excluded from output VCF as might be artificially 
# introduced by chunks aligning to same reference regions

#
# J Sarria, B Contreras-Moreira
# Copyright [2026] Estacion Experimental de Aula Dei-CSIC

my ($cmd,$root,$key,$val,%opts,%config,@temp);
my ($fqfile,$output_file,$redo,$allsites) = ('', '', 0, 0);
my ($threads,$chunksize,$minmatch,$outdir) = (5, 500, 101, '/tmp');
my $samtoolsEXE = 'samtools';
my $bwaEXE      = 'minibwa';
my $bedtoolsEXE = 'bedtools';
my $bcftoolsEXE = 'bcftools';

getopts('hARB:f:c:o:m:t:k:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)) {
  print "Maps reads in FASTQ file to imputate and call variants in output gVCF file\n\n";
  print "Usage: $0 [options]\n\n";
  print "-h this message\n";
  print "-f input FASTQ file                                (example: -f sample.fastq.gz)\n";
  print "-c graph config file                               (example: -c /path/to/graph/config.impute.yaml)\n";
  print "-o output folder                                   (optional, example: -o mysample, default -o $outdir)\n";
  print "-m min match length                                (optional, example: -m 150, default -m $minmatch)\n";
  print "-t threads                                         (optional, example: -t 12, default -d $threads)\n";
  print "-k chunk size in bases                             (optional, example: -k 1000, default -k $chunksize\n";
  print "-B path to [mini]bwa binary                        (optional, example: -B /path/to/[mini]bwa)\n";
  print "-A all sites in output VCF, not just variants      (optional, by default only variants are considered)\n";
  print "-R redo all steps even if results are in place     (optional, by default previous results are re-used)\n";
  #print "\nPrimary citation:\n";
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

if(defined($opts{'B'})) {
  $bwaEXE = $opts{'B'}
}

if(defined($opts{'A'})) {
  $allsites = 1
}

if(defined($opts{'R'})) {
  $redo = 1
}

$output_file = "$outdir/$root" . ".$chunksize.vcf.gz";

warn "## $0 -f $fqfile -c $opts{'c'} -o $outdir -m $minmatch -t $threads -k $chunksize -B $bwaEXE -A $allsites -R $redo\n\n";

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
push(@temp, $mapfile);
if($redo == 1 || !-e $mapfile) {
  $cmd = "phg map-reads --index $config{'kmer-index-file'} --read-files $fqfile -o $outdir --hvcf-dir $config{'hvcf-dir'} " .
    "--threads $threads --min-mem-length $minmatch "; #--conda-env-prefix $config{'conda-env-prefix'}";
  run_cmd($cmd, "# 1 Running phg map-reads ...");

} else {
  print "# 1 re-using $mapfile\n";
}	

# 2) find paths
my $hvcf_file = "$outdir/$root" . '_1.h.vcf';
push(@temp, $hvcf_file);
if($redo == 1 || !-e $hvcf_file) {
  $cmd = "phg find-paths --read-files $mapfile --output-dir $outdir --hvcf-dir $config{'hvcf-dir'} --path-type haploid " .
    "--threads $threads --reference-genome $config{'reference-fasta'}";
  run_cmd($cmd, "# 2 Running phg find-paths ...");

} else {
  print "# 2 re-using $hvcf_file\n";
}

# 3) create fasta from hvcf
my $composite_file = "$outdir/$root" . '_1_composite.fa';
push(@temp, $composite_file);
if($redo == 1 || !-e $composite_file) {
  $cmd = "phg create-fasta-from-hvcf --hvcf-file $hvcf_file -o $outdir --db-path $config{'db-path'} " .
    "--range-bedfile $config{'range-bed'} --fasta-type composite"; #--conda-env-prefix $config{'conda-env-prefix'}";
  run_cmd($cmd, "# 3 Running phg create-fasta-from-hvcf ...");
  
} else {
  print "# 3 re-using $composite_file\n";
}

# 4) index fasta
my $composite_index_file = $composite_file . '.fai';
push(@temp, $composite_index_file);
if($redo == 1 || !-e $composite_index_file) {
  $cmd = "$samtoolsEXE faidx $composite_file";
  run_cmd($cmd, "# 4 Indexing composite fasta ...");

} else {
  print "# 4 re-using $composite_index_file\n";
}

# 5) chunk composite sequence
my $composite_chunks_bam = "$outdir/$root" . "_1_composite.$chunksize.bam";
push(@temp, $composite_chunks_bam);
if($redo == 1 || !-e $composite_chunks_bam) {

  my $comp_genome_bed  = "$outdir/$root" . '_1_comp_genome.bed';
  my $comp_chunks_bed  = "$outdir/$root" . '_1_comp_chunks.bed';
  my $comp_chunks_file = "$outdir/$root" . '_1_comp_chunks.fa';
  my $ref_index_file   = $config{'reference-fasta'} . '.mbw';
  my @temp2 = ($comp_genome_bed, $comp_chunks_bed, $comp_chunks_file, $ref_index_file);

  $cmd = "perl -lane 'print \"\$F[0]\t0\t\$F[1]\"' $composite_index_file > $comp_genome_bed";
  run_cmd($cmd) if(!-e $comp_genome_bed || $redo == 1);
  
  $cmd = "$bedtoolsEXE makewindows -b $comp_genome_bed -w $chunksize > $comp_chunks_bed";
  run_cmd($cmd) if(!-e $comp_chunks_bed || $redo == 1);
  
  $cmd = "$bedtoolsEXE getfasta -fi $composite_file -bed $comp_chunks_bed -fo $comp_chunks_file";
  run_cmd($cmd, "# 5.1 Chunking composite fasta ...") if(!-e $comp_chunks_file || $redo == 1);

  $cmd = "$bwaEXE index -t $threads $config{'reference-fasta'}";
  if($bwaEXE !~ /minibwa/) {
    $ref_index_file   = $config{'reference-fasta'} . '.sa';
    $cmd = "$bwaEXE index $config{'reference-fasta'}";
  }
  run_cmd($cmd, "# 5.2 Indexing reference fasta ...") if(!-e $ref_index_file || $redo == 1);

  $cmd = "$bwaEXE map -t $threads $config{'reference-fasta'} $comp_chunks_file | " .
    "$samtoolsEXE sort -@ $threads -o $composite_chunks_bam && $samtoolsEXE index -c -@ $threads $composite_chunks_bam";
  if($bwaEXE !~ /minibwa/) {
    $cmd = "$bwaEXE mem -t $threads $config{'reference-fasta'} $comp_chunks_file | " .
    "$samtoolsEXE sort -@ $threads -o $composite_chunks_bam && $samtoolsEXE index -c -@ $threads $composite_chunks_bam";
  }
  run_cmd($cmd, "# 5.3 Mapping composite fasta ...");

  unlink(@temp2);
  
} else {
  print "# 5 re-using $composite_chunks_bam\n";
}

# 6) variant call, produces final gVCF output file with het sites filtered out
if($redo == 1 || !-e $output_file) {
  $cmd = "$bcftoolsEXE mpileup -Ou -f $config{'reference-fasta'} $composite_chunks_bam | ";
  if($allsites == 1) {
    $cmd .= "$bcftoolsEXE call -m --ploidy 2 -Oz | bcftools view -g ^het -o $output_file";
  } else {
    $cmd .= "$bcftoolsEXE call -m --ploidy 2 -v -Oz | bcftools view -g ^het -o $output_file";
  }

  run_cmd($cmd, "# 6 Variant calling ...");
  print "## output: $output_file\n";

} else {
  print "## output $output_file (re-used)\n";
}

# clean only if successful
unlink(@temp);

exit(0);




sub run_cmd {
  my ($cmd, $message) = @_;

  print $message if($message);
  system($cmd);
  if($? != 0) {
      die "# EXIT: command failed ($cmd)\n";
  }
}
