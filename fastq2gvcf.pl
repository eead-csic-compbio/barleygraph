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

my ($cmd,$root,%opts);
my ($fqfile,$hvcfdir,$rfile,$kindexfile,$redo) = ('', '', '', '', 0);
my ($threads,$chunksize,$minmatch,$outdir,$condaenvpath) = (5, 500, 101, '/tmp', '');

getopts('hRf:r:i:v:o:m:t:k:c:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)) {
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-f input FASTQ file                                (example: -f sample.fastq.gz)\n";
  print "-r reference FASTA file                            (example: -r MorexV3.fna)\n";
  print "-i kmer index file                                 (example: -i Pan20_proali.fmd)\n";
  print "-v graph hvcf-dir                                  (must contain h.vcf.gz files)\n";
  print "-o output folder                                   (optional, example: -o mysample, default -o /tmp)\n";
  print "-c conda phg path                                  (optional, example: -c /home/user/miniconda3/envs/phgv2.4)\n";
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

if(!defined($opts{'r'})) {
  die "# ERROR: need reference FASTA file (-r)\n";
} else {
  $rfile = $opts{'r'};
}

if(!defined($opts{'i'})) {
  die "# ERROR: need kmer index file (-i)\n";
} else {
  $kindexfile = $opts{'i'};

  if(!-e $kindexfile.'.ssa') {
    die "# ERROR: need also kmer index file $kindexfile\.ssa\n";
  } elsif(!-e $kindexfile.'.len.gz') {
    die "# ERROR: need also kmer index file $kindexfile\.len.gz\n";
  }  
}

if(!defined($opts{'v'})) {
  die "# ERROR: need valid hvcf-dir with h.vcf.gz files (-v)\n";
} else {
  $hvcfdir = $opts{'v'};

  opendir(HVCF,$hvcfdir) || 
    die "# ERROR: cannot list $hvcfdir\n";
  my @files = grep{/\.h.vcf.gz/} readdir(HVCF);
  closedir(HVCF);

  if(scalar(@files) == 0) {
    die "# ERROR: cannot find h.vcf.gz files in $hvcfdir\n";
  }
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

if(defined($opts{'c'})) {
  $condaenvpath = $opts{'c'};
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

warn "## $0 -f $fqfile -r $rfile -i $kindexfile -o $outdir -c $condaenvpath -m $minmatch -t $threads -k $chunksize -R $redo\n\n";

#INDEX="$PAN20_ROOT/output/vcf_files/proali/Pan20_proali.fmd"
#HVCF_DIR="$PAN20_ROOT/output/vcf_files/proali"
#REFERENCE_GENOME="$PAN20_ROOT/data/MorexV3.fa"
#DB_PATH="$PAN20_ROOT/vcf_dbs"
#RANGE_BEDFILE="$PAN20_ROOT/output/ref_ranges.bed"
#MOREX_REF="/agave/compbio/references/barley/Barley_Morex_V3/GCA_904849725.1_MorexV3_pseudomolecules.chrnames.fna"
#ASSEMBLY_DIR="/agave/compbio/references/barley/pangenome56/Assemblies"
#MIN_MEM_LENGTH=101

#DATA_DIR="$SCRIPT_DIR/data"
#READMAPPINGS_DIR="$SCRIPT_DIR/readmappings"
#IMPUTED_DIR="$SCRIPT_DIR/imputed_hvcf"
#COMPOSITE_DIR="$SCRIPT_DIR/composite_fastas"
#ALIGNMENTS_DIR="$SCRIPT_DIR/alignments"
#FLAGSTAT_DIR="$SCRIPT_DIR/flagstat"
#LOG_DIR="$SCRIPT_DIR/logs"
#VCF_DIR="$SCRIPT_DIR/vcf"
#mkdir -p "$READMAPPINGS_DIR" "$IMPUTED_DIR" "$COMPOSITE_DIR" "$ALIGNMENTS_DIR" "$FLAGSTAT_DIR" "$LOG_DIR" "$VCF_DIR"

# 0) check phg is loaded (conda)
if(!`which phg`) {
  die "# ERROR: cannot find phg in \$PATH, perhaps you need to conda activate it?\n";  
}

# 1) map reads
my $mapfile = "$outdir/$root" . '_1_readMapping.txt';
if($redo == 1 || !-e $mapfile) {
  $cmd = "phg map-reads --index $kindexfile --read-files $fqfile -o $outdir --hvcf-dir $hvcfdir " .
    "--threads $threads --min-mem-length $minmatch";
  if($condaenvpath ne '') { 
    $cmd .= " --conda-env-prefix $condaenvpath";
  }  

  warn "# 1 Running phg map-reads ...";
  system($cmd);
  if($? != 0) {
      die "# EXIT: failed ($cmd)\n";
  }
} else {
  print "# re-using $mapfile\n";
}	


