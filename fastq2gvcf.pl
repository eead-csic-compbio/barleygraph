#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Std;
use File::Basename qw/ basename /; 

# Maps reads in FASTQ file to indexed PHG pangenome graph in order to imputate and call
# haplotypes across graph ranges, producing a PHG-imputed FASTA sequence which is ultimately
# mapped back to the graph reference (MorexV3 in barley) to call variants a produce a gVCF file
#
# J Sarria, B Contreras-Moreira
# Copyright [2026] Estacion Experimental de Aula Dei-CSIC

my ($cmd,$root,$ref,%opts);
my ($fqfile,$rfile,$kindexfile,$redo) = ('', '', '', 0);
my ($threads,$chunksize,$minmatch,$outdir) = (4, 500, 101, '/tmp');

getopts('hRf:r:i:o:m:t:k:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)) {
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-f input FASTQ file                                (example: -f sample.fastq.gz)\n";
  print "-r reference FASTA file                            (example: -r MorexV3.fna)\n";
  print "-i kmer index file                                 (example: -i Pan20_proali.fmd)\n";
  print "-o output folder                                   (optional, example: -o mysample, default -o /tmp)\n";
  print "-m min match length                                (optional, example -m 150, default -m $minmatch)\n";
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
  $root = basename($fqfile);
}

if(!defined($opts{'r'})) {
  die "# ERROR: need reference FASTA file (-r)\n";
} else {
  $rfile = $opts{'r'};
  $ref = basename($rfile);
}

if(!defined($opts{'k'})) {
  die "# ERROR: need kmer index file (-k)\n";
} else {
  $kindexfile = $opts{'k'}
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

warn "## $0 -f $fqfile -r $rfile -i $kindexfile -o $outdir -m $minmatch -t $threads -k $chunksize -R $redo\n\n";

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

# check phg is loaded (conda)
if(!`which phg`) {
  die "# ERROR: cannot find phg in \$PATH, perhaps you need to conda activate it?\n";  
}


# [1/6] map-reads
my $mapfile = $root .'.'. $ref . '.readMapping.txt';
if($redo == 1 || !-e $mapfile) {
	#$cmd = run_cmd phg map-reads --index "$INDEX" --read-files "$fq" -o "$READMAPPINGS_DIR" --hvcf-dir "$HVCF_DIR" --threads "$THREADS" --min-mem-length "$MIN_MEM_LENGTH"	
  #warn "# [1/6] Running phg map-reads..."
} else {

}	

