#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Std;
use Cwd 'abs_path';

# Download, unpack and setup graphs to be used in Docker container.
# Uses the following system's binaries: cat, wget, md5sum
#
# J Sarria, B Contreras-Moreira
# Copyright [2026] Estacion Experimental de Aula Dei-CSIC

# Due to their large size, graphs are split in parts; to add new graphs you need:
# 1) URL of the first part (00), 2) natural of the last part, 3) md5sum of joined parts
my %graphs;
$graphs{'Pan20-mmap-pro'}{'subfolder'} = 'Pan20';
$graphs{'Pan20-mmap-pro'}{'URL'} = 
  'https://github.com/eead-csic-compbio/barleygraph/releases/download/Pan20-mmap-pro-1.0/Pan20-mmap-pro00.part';
$graphs{'Pan20-mmap-pro'}{'last'} = 12;
$graphs{'Pan20-mmap-pro'}{'md5sum'} = '10b88c32b8cd01dfc214487ed1fe7cae';

# should match Docker
my $target_path = abs_path('graph_db/');

my ($graph,$tgzfile,$part,$partfile,$url,$cmd,$sum,$path);
my (%opts,@temp);

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
    print "$graph => $graphs{$graph}{'URL'}\n";
  }
  exit(0);
}

if(defined($opts{'G'})) {

  $graph = $opts{'G'};
  $tgzfile = "$graph.tgz";

  if(!defined($graphs{$graph})) {
    die "# ERROR: unsupported graph ($graph), please run $0 -l to see available graphs\n"; 
  }	  

  if(!-e $tgzfile) {

    # expand list of URLS for all parts
    foreach $part (0 .. $graphs{$graph}{'last'}) {
      $url = $graphs{$graph}{'URL'};
      if(length($part) == 1) {
        $url =~ s/00.part/0$part.part/
      } else {
        $url =~ s/00.part/$part.part/
      } 
      $partfile = (split(/\//,$url))[-1];
      push(@temp, $partfile);
  
      $cmd = "wget -q -c $url";
      run_cmd($cmd, "# 1 Downloading $url");
    }
  
    $cmd = "cat $graph*.part > $tgzfile";
    run_cmd($cmd, "# 2 Joining graph parts ...");

    print "# 3 Checking download sum ...\n";
    $sum = `md5sum $tgzfile`;
    if($sum =~ /$graphs{$graph}{'md5sum'}/) {
       unlink(@temp);
       print "# Download correct\n";	    
    } else {
      push(@temp, $tgzfile);
      unlink(@temp);
      die "# Download failed, incorrect checksum, please re-try\n";
    }
  }

  $path = "$target_path/$graphs{'Pan20-mmap-pro'}{'subfolder'}";
  if(!-d $path) {
    if(!mkdir($path)) {
      die "# ERROR: cannot create folder $path , check $target_path exists or permissions\n";
    }
  }

  $cmd = "tar xvfz $tgzfile -C $path";
  run_cmd($cmd, "# 3 Unpacking graph ..."); 

  # last operations on unpacked data
  
  # required by create-fasta-from-hvcf (imputation -g)
  mkdir("$path/$graph/vcf_dbs");
  symlink("$path/assemblies.agc", "$path/$graph/vcf_dbs/assemblies.agc");

  print "# Setup complete ($graph)\n";
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
