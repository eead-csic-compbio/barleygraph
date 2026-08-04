#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Std;

# Download, unpack and setup graph to be used in Docker container
#
# J Sarria, B Contreras-Moreira
# Copyright [2026] Estacion Experimental de Aula Dei-CSIC

# To add a new graph in the future, simply add a line with this format:
# 'GRAPHNAME' => 'WGET PATTERN'
my %graphs = (
  'Pan20-mmap-pro' =>  
    'https://github.com/eead-csic-compbio/barleygraph/releases/download/Pan20-mmap-pro-1.0/Pan20-mmap-pro{00..12}.part',
);

my $target_path = '/barleygraph_databases/';

my ($graph,$tgzfile,$partial_file_pat,$cmd,%opts,@temp);

getopts('hlG:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)) {
  print "Usage: $0 [options]\n\n";
  print "-h this message\n";
  print "-l list available pangenome graphs          (optional)\n";
  print "-G graph name, should be in supported list  (required, example: -G Pan20-mmap-pro)\n";
  #print "\nPrimary citation:\n";
  exit(0);
}

if(defined($opts{'l'})) {
  foreach $graph (sort keys(%graphs)) {
    print "$graph => $graphs{$graph}\n";
  }
  exit(0);
}

if(defined($opts{'G'})) {

  $graph = $opts{'G'};
  $tgzfile = "$graph.tgz";

  if(!defined($graphs{$graph})) {
    die "# ERROR: unsupported graph ($graph), please run $0 -l to see available graphs\n"; 
  }	  
  
  if(-e $tgzfile) { 
    $cmd = "wget $graphs{$graph}";
    run_cmd($cmd, "# 1 Downloading graph parts ...");
    @temp  = glob("$graph.part");

    $cmd = "cat $graph* > $tgzfile && md5sum $tgzfile";
    run_cmd($cmd, "# 2 Joining graph parts ...");
  }

  $cmd = "tar xvfz $tgzfile -C $target_path";
  run_cmd($cmd, "# 3 Unpacking graph ..."); 
  #push(@temp, $tgzfile);

  print "# Setup complete for $graph\n";
  # clean only if successful
  unlink(@temp);
}

exit(0);




sub run_cmd {
  my ($cmd, $message) = @_;

  print "$message\n" if($message);
  system($cmd);
  if($? != 0) {
      die "# EXIT: command failed ($cmd)\n";
  }
}

#PARENT_DIR="$DATABASE_ROOT/$PANGENOME_NAME"
#TAR_FILE="$DATABASE_ROOT/${DATASET_NAME}.tgz"
#EXTRACT_DIR="$PARENT_DIR/$DATASET_NAME"

#echo "Using Database Root: $DATABASE_ROOT"

# Create the parent directory if it doesn't exist
#if ! mkdir -p "$PARENT_DIR" 2>/dev/null; then
#    echo "ERROR: Failed to create directory $PARENT_DIR."
#    echo "Do you have write permissions? You may need to run this script with 'sudo'."
#    exit 1
#fi

#  mkdir -p "$EXTRACT_DIR"
  
#  # Extract directly into the temporary folder
#  if ! tar -xzf "$TAR_FILE" -C "$EXTRACT_DIR"; then
#    echo "ERROR: Failed to extract ${DATASET_NAME}.tgz."
#    exit 1
#  fi

#  # Merge contents into the parent folder (skipping existing/duplicate files)
#  echo "Merging contents into shared $PANGENOME_NAME folder (skipping duplicate files)..."
  
#  # cp -rn recursively merges folders and copies files without overwriting.
#  cp -rn "$EXTRACT_DIR"/* "$PARENT_DIR"/ 2>/dev/null || true
  
  # Remove the temporary folder
  #  rm -rf "$EXTRACT_DIR"
