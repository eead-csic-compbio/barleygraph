#!/usr/bin/env python3

# Aligns arbitrary DNA sequences to a PHG pangenome graph:
# i) uses gmap to hierarchically match sequences to contributing genomes, 
# ii) queries PHG haplotypes to obtain precomputed matches in other genomes
#
#TODO: cambiar si sumamos paso pangenes
#
# J Sarria, B Contreras-Moreira
# Copyright [2024-25] Estacion Experimental de Aula Dei-CSIC

# %%
def parseFASTASeqs(fasta):
    """Takes a FASTA genome filename and parses individual sequences.
    Returns i) list of sequence names in input order, ii) dictionary with sequences."""

    names = []
    sequences = {}
    name = ''

    try:
        file = open(fasta)
    except OSError as error:
        print("# ERROR: cannot open/read file:", fasta, error)
        return [],{}

    for line in file:
        header = re.search(r"^>", line)
        if header:
            seq_name_match = re.search(r"^>(\S+)", line)
            if seq_name_match:
                name = seq_name_match.group(1)
                names.append(name)
            else:
                print("# ERROR: cannot parse FASTA header:", header)
        else:
            if name in sequences:
                sequences[name] = sequences[name] + line
            else:
                sequences[name] = line

    return names, sequences

# %%
def checkGmapDBVersion(gmap_db, ref_name):
    """Returns version of reference Gmap db."""

    version_db = '?'

    gmap_version_file=f'{gmap_db}/{ref_name}/{ref_name}.version'
    
    if not os.path.isfile(gmap_version_file):
        print(f"# ERROR(checkGmapVersion: file {gmap_version_file} does not exist")
        return version_exe, version_db

    else:
        with open (gmap_version_file) as f:
            for line in f:
                if not line.startswith(ref_name):
                    version_db = line.strip()

    return version_db

# %%
def checkGmapVersion(gmap_path):
    """Returns version of Gmap binary/executable."""

    version_exe = '?'
    command = f"{gmap_path} --version"
    try:
        result = subprocess.run(command,
                                shell=True,text=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.DEVNULL)

        #Part of GMAP package, version 2024-11-20
        data = result.stdout.splitlines()
        version_exe = data[2]

    except subprocess.CalledProcessError as e:
        print(f'# ERROR(checkGmapVersion): {e.cmd} failed: {e.stderr}')

    return version_exe

# %%
def getTempFASTAFile(query_name, fasta_sequence, prefix_string="temp",
                     temp_path="/tmp/", suffix_string=".fna",
                     verbose=False):
    """Takes up to 5 strings + 1 Boolean: 
    i) query name, ii) sequence, iii) temp prefix, iv) temp path, v) temp suffix/extension, 
    vi) verbose.
    Creates a temporary FASTA file with passed query sequence.
    Returns path to created temp file."""
   
    with tempfile.NamedTemporaryFile(prefix=prefix_string,suffix=suffix_string,
                                     dir=temp_path,delete=False) as temp:
        temp.write(b">" + query_name.encode() + b"\n" + fasta_sequence.encode()) 
        temp.close()

    file_path = Path(temp.name)
    if file_path.exists() and file_path.stat().st_size > 0:
        if verbose == True:
            print(f"# INFO(getTempFASTAFile) created {temp.name}")
        return temp.name
    else:
        print("# ERROR(getTempFASTAFile): cannot write temp file")
        return ''

# %%
def validMatch(gff_file, min_identity, min_coverage, verbose=False):
    """Checks if GFF3 contains gene & mRNA features with required identity and coverage.
    Returns: i) Boolean, ii) float ident, iii) float cover"""

    gene_match = False
    gmap_match = False 
    identity = 0.0
    coverage = 0.0

    if not os.path.isfile(gff_file):
        print(f"# ERROR(validMatch): file {gff_file} does not exist")
        return gmap_match
    
    else:
        with open (gff_file) as f:
            for line in f:
                if not line.startswith('#'):
                    fields = line.split("\t")

                    if fields[2] == "gene":
                        gene_match = True
                        continue
                    
                    #ID=query.mrna1;Name=name;Parent=genome.path1;Dir=na;coverage=100.0;identity=98.9;
                    elif fields[2] == "mRNA" and gene_match == True:
                        attributes = fields[-1].split(";")  
                        coverage = float(attributes[4].split("=")[1])   
                        identity = float(attributes[5].split("=")[1])   
                        if identity >= min_identity and coverage >= min_coverage:
                            gmap_match = True
                            break

        if verbose == True:                    
            print(f"# Match = {gmap_match} (ident={identity}, cover={coverage})")                                                                            
    return gmap_match, identity, coverage


# %%
def checkRangesKeyInRef(reference_gff3,hapIDranges,bed_folder_path,coverage=75.0,
                        bedtools_path='bedtools',grep_path='grep',verbose=False):
    """Retrieves PHG keys for range overlapping gmap match in reference genome.
    Passed coverage is used to intersect ranges and match. Overlap does not consider strandness. 
    Returns: i) string with matched coords in TSV format, ii) dictionary of keys/genomic ranges
    of all matched genomes (useful to compute alignments down the road).
    Column order in TSV: graph_chr, graph_start, graph_end, graph_strand (. if absent in ref), 
    multiple_mappings (Yes/No), match_genome, match_chr, match_start, match_end, match_strand"""

    keys = {}
    match_tsv = ''
    mult_mappings = 'No'
    chrom, genome, start, end, strand = '','','','',''

    with open(reference_gff3) as f:
        for line in f:
            if not line.startswith('#'):
                fields = line.split("\t")
                if fields[2] == "gene":
                    if genome == '':
                        chrom = fields[0]
                        genome = fields[1]
                        start = fields[3]
                        end = fields[4]
                        strand = fields[6]
                    else:
                        mult_mappings = 'Yes'
                        break
    f.close()                

    if verbose == True:
        print(f"# Checking match for {chrom}:{start}-{end} at reference")

    # read genome list from hapIDranges
    genomes = []
    with open(hapIDranges) as f:
        for line in f:
            if line.startswith('#'):
                genomes = line.split("\t")
                genomes = genomes[3:]
                genomes[-1] = genomes[-1].strip()
                break
    f.close()

    # prepare bedtools intersect command to find overlapping range
    command = f"{bedtools_path} intersect -a {hapIDranges} -b stdin -nonamecheck -F 0.5"             

    # BED-format interval of gmap match
    match_interval = f'{chrom}\t{start}\t{end}'

    # actually run the bedtools command   
    try:
        result = subprocess.run(command,
                                shell=True,check=True,text=True,input=match_interval,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                        
        intersections = result.stdout.splitlines()
        if verbose == True:
            print(f"# INFO(checkRangesKeyInRef): {result.stdout}")

        if len(intersections) > 1:
            if(verbose == True):
                print(f"# WARN(checkRangesKeyInRef): more than 1 range overlaps: {intersections}")
            intersections = intersections[0] 
                    
        for feature in intersections:
            feature = str(feature).split("\t")
            feature[-1] = feature[-1].strip() 
            match_tsv = (
                f'{feature[0]}\t{feature[1]}\t{feature[2]}\t{strand}\t{mult_mappings}'
                f'\t{genome}\t{feature[0]}\t{feature[1]}\t{feature[2]}\t{strand}')
            
            for c in range(0,len(genomes)):
                k = feature[c+3]
                if k == ".":
                    continue
                else:
                    clean_k = k[1:-1] #remove <>
                    if clean_k not in keys:
                        keys[clean_k] = genomes[c]

    except subprocess.CalledProcessError as e:
        print(f'# ERROR(checkRangesKeyInRef): {e.cmd} failed: {e.stderr}')

    # retrieve genomic ranges matching these keys
    for k in keys:

        range_bedfile = f"{bed_folder_path}/{keys[k]}.h.bed"
        command = f"{grep_path} {k} {range_bedfile}"
        try:
            result = subprocess.run(command,
                                shell=True,check=True,text=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)

            bed_data = result.stdout.splitlines()
            if len(bed_data) > 1:
                if(verbose == True):
                    print(f"# WARN(checkRangesKeyInRef): more than 1 key matches: {results.stdout}")
                bed_data = bed_data[0]
            
            bed_data = bed_data[0].split("\t")
            keys[k] = f'{keys[k]}:{bed_data[0]}:{bed_data[1]}-{bed_data[2]}({bed_data[3]})'

        except subprocess.CalledProcessError as e:
            print(f'# ERROR(checkRangesKeyInRef): {e.cmd} failed: {e.stderr}')

    return match_tsv, keys


# %%
def sortGenomesByRangeNumber(folder_path, hap_table_file, verbose=False):
    """Sorts the genomes in pangenome folder by number of ranges
    by parsing hap_table_file (max to min). 
    Returns list with genome names, without .fa extension."""

    # count how many ranges are supported by each genome contributing haplotypes
    hapnum = {}
    with open(hap_table_file) as hap:
        for line in hap:
            cols = line.split("\t")
            for g in cols[1].strip().split(","):
                if not g in hapnum:
                    hapnum[g] = 1
                else:
                    hapnum[g] += 1
    hap.close()
    
    # sort genomes and check they are in folder
    pangenome_genomes = []
    for g in sorted(hapnum, key=hapnum.get, reverse=True):
        if(os.path.isfile(f'{folder_path}/{g}.fa')):
            pangenome_genomes.append(g)
            if verbose == True:
                print(f'# {g} {hapnum[g]}')

    #pangenome_genomes = ['HOR_13942'] #debug

    return pangenome_genomes

# %%
def runGmapAgainstGenomes(pangenome_genomes, gmap_path, gmap_db, fasta_filename, 
                         min_identity, min_coverage,
                         cores=4, prefix='temp', path='/tmp/', verbose=False):
    """Iteratively gmaps FASTA file against genomes included in pangenome.
    Returns i) Boolean match, ii) string with matched genome name, iii) path to GFF3 file."""
    
    for genome in pangenome_genomes:

        gff_filename = f"{path}/{genome}.{prefix}.gff3"
        if verbose == True:
            print(f"# Running gmap against {genome}")

        # try dafault gmap 
        pangenome_gmap_command = (
            f"{gmap_path} -D {gmap_db} -d {genome} -t {cores} {fasta_filename} -f gff3_gene > {gff_filename}")

        if verbose == True:
            print(pangenome_gmap_command)

        try:
            result = subprocess.run(pangenome_gmap_command, shell=True, check=True, 
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if verbose == True:
                print(result.stdout.decode())

        except subprocess.CalledProcessError as e:

            # regex to decide whether gmapl is needed (genomes > 4Gbp)
            pattern = r'For big genomes of more than'
            match = re.search(pattern, e.stderr.decode())

            # try gmapl instead
            if match:
                if verbose == True:
                    print(f"# Running gmapl for {genome}")

                pangenome_gmap_command = (
                    f"{gmap_path}l -D {gmap_db} -d {genome} -t {cores} {fasta_filename} -f gff3_gene > {gff_filename}")

                try:
                    result = subprocess.run(pangenome_gmap_command, shell=True, check=True, 
                                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    if verbose == True:
                        print(result.stdout.decode())

                except subprocess.CalledProcessError as e:
                    print(f"\nERROR(runGmapAgainsGenomes): '{e.cmd}' returned non-zero exit status {e.returncode}.")

            else:    
                print(f"\nERROR(runGmapAgainsGenomes): '{e.cmd}' returned non-zero exit status {e.returncode}.")
            
        gmap_match, ident, cover = validMatch(gff_filename,min_identity,min_coverage,verbose=verbose)
        if gmap_match == True:            
            break
        else:
            os.remove(gff_filename)    

    return gmap_match, genome, gff_filename

# %%
def checkRangesKeysInPangenome(gff3,hapIDranges,bedfile,bed_folder_path,coverage=75.0,
                               bedtools_path='bedtools',grep_path='grep',verbose=False):
    """Retrieves PHG keys for range overlapping gmap match in 1st matched pangenome assembly.
    BED file is usually a .h.bed file with sorted ranges extracted from PHG .h.vcf.gz files.
    Passed coverage is used to intersect ranges and match. Overlap does not consider strandness.
    Returns: i) string with matched coords in TSV format, ii) dictionary of keys/genomic ranges
    of all matched genomes (useful to compute alignments down the road).
    Column order in TSV: graph_chr, graph_start, graph_end, graph_strand (. if absent in ref),
    multiple_mappings (Yes/No), match_genome, match_chr, match_start, match_end, match_strand"""

    keys = {}
    match_tsv = ''
    graph_key = ''
    mult_mappings = 'No'
    chrom, genome, start, end, strand = '','','','',''

    with open(gff3) as f:
        for line in f:
            if not line.startswith('#'):
                fields = line.split("\t")
                if fields[2] == "gene":
                    if genome == '':
                        chrom = fields[0]
                        genome = fields[1]
                        start = fields[3]
                        end = fields[4]
                        strand = fields[6]
                    else:
                        mult_mappings = 'Yes'
                        break
    f.close()

    if verbose == True:
        print(f"Checking match for {chrom}:{start}-{end} at {bedfile}")

    # read genome list from hapIDranges
    genomes = []
    with open(hapIDranges) as f:
        for line in f:
            if line.startswith('#'):
                genomes = line.split("\t")
                genomes = genomes[3:]
                genomes[-1] = genomes[-1].strip()
                break
    f.close()

    # prepare bedtools intersect command to find overlapping range,
    # bedfile should contain lines like this:
    # chr1H_OX460222.1 1 69 + 9c51... HOR_12184 chr1H_LR890096.1 9 66 21c7...
    command = f"{bedtools_path} intersect -a {bedfile} -b stdin -nonamecheck -F 0.5"             

    # BED-format interval of gmap match
    match_interval = f'{chrom}\t{start}\t{end}'

    # actually run the bedtools command
    try:
        result = subprocess.run(command,
                                shell=True,check=True,text=True,input=match_interval,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                        
        intersections = result.stdout.splitlines()
        if verbose == True:
            print(f"# INFO(checkRangesKeyInPangenome): {result.stdout}")

        if len(intersections) > 1:
            if(verbose == True):
                print(f"# WARN(checkRangesKeyInPangenome): more than 1 range overlaps: {intersections}")
            intersections = intersections[0] 
                    
        for feature in intersections:
            feature = str(feature).split("\t")
            feature[-1] = feature[-1].replace("\n", "") 
            match_tsv = ( # strand unknown as match is missing from reference genome
                f'{feature[6]}\t{feature[7]}\t{feature[8]}\t.\t{mult_mappings}'
                f'\t{genome}\t{feature[0]}\t{feature[1]}\t{feature[2]}\t{feature[3]}')
            graph_key = feature[4]
            
        # look for this key within graph ranges (grep)
        command = f"{grep_path} {graph_key} {hapIDranges}"
        try:
            result = subprocess.run(command,
                                shell=True,check=True,text=True,input=match_interval,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE)

            graph_data = result.stdout.splitlines()
            if len(graph_data) > 1:
                if(verbose == True):
                    print(f"# WARN(checkRangesKeyInPangenome): more than 1 graph matches: {results.stdout}")
                graph_data = graph_data[0]
 
            feature = graph_data[0].split("\t")
            feature[-1] = feature[-1].strip()
            for c in range(0,len(genomes)):
                k = feature[c+3]
                if k == ".":
                    continue
                else:
                    clean_k = k[1:-1] #remove <>
                    if clean_k not in keys:
                        keys[clean_k] = genomes[c]

        except subprocess.CalledProcessError as e:
            print(f'# ERROR(checkRangesKeyInPangenome): {e.cmd} failed: {e.stderr}')
                            
    except subprocess.CalledProcessError as e:
        print(f'ERROR(checkRangesKeyInPangenome): {e.cmd} failed: {e.stderr}')

    # retrieve genomic ranges matching these keys
    for k in keys:

        range_bedfile = f"{bed_folder_path}/{keys[k]}.h.bed"
        command = f"{grep_path} {k} {range_bedfile}"
        try:
            result = subprocess.run(command,
                                shell=True,check=True,text=True,input=match_interval,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE)

            bed_data = result.stdout.splitlines()
            if len(bed_data) > 1:
                if(verbose == True):
                    print(f"# WARN(checkRangesKeyInPangenome): more than 1 key matches: {results.stdout}")
                bed_data = bed_data[0]

            bed_data = bed_data[0].split("\t")
            keys[k] = f'{keys[k]}:{bed_data[0]}:{bed_data[1]}-{bed_data[2]}({bed_data[3]})'

        except subprocess.CalledProcessError as e:
            print(f'# ERROR(checkRangesKeyInPangenome): {e.cmd} failed: {e.stderr}')

    return match_tsv, keys



# %%
def main():

    # hard-coded defaults
    grep_exe = 'grep' # assumed to be available
    keepTempFiles  = False

    parser = argparse.ArgumentParser(
        description="Tool to map sequences in barley pangenome.\n",
        epilog="Citation:\n",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "config_file",
        help="path to YAML file with pangenome config (required)",
    )

    parser.add_argument(
        "fasta_file",
        help="path to FASTA file with sequences to map (required)"
    )

    parser.add_argument(
        "--tmp_path",
        default='/tmp/',
        help="path to writable folder for temporary files, default: /tmp"
    )

    parser.add_argument(
        "--bedtools_exe",
        default='bedtools',
        help="path to bedtools executable, default: bedtools"
    )

    parser.add_argument(
        "--cor", default=4, help="number of cores for gmap, default: 4"
    )

    parser.add_argument(
        "--minident",
        default=98.0,
        help="min %identity of gmap matches, default: 98.0",
    )

    parser.add_argument(
        "--mincover",
        default=95.0,
        help="min %coverage of gmap matches and pangenome ranges, default: 95.0",
    )

    parser.add_argument('--verb', help='Increase verbosity in output')

    parser.add_argument('--add_ranges', help='Add all pangenome ranges matching input sequences')

    args = parser.parse_args()

    fasta_file = args.fasta_file

    # parse pangenome config file
    with open(args.config_file, "r") as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

        vcf_dbs = config['vcf_dbs']
        gmap_db = config['gmap_db']
        reference_name = config['reference_name']
        pangenome_fastas_folder = config['pangenome_fastas']
        hapIDtable  = config['hapIDtable']
        hapIDranges = config['hapIDranges']
        gmap_exe = config['gmap_exe']

    # get optional params
    bedtools_exe = args.bedtools_exe
    ncores       = args.cor
    min_identity = args.minident
    min_coverage = args.mincover
    temp_path    = args.tmp_path
    verbose_out  = args.verb
    add_ranges   = args.add_ranges


    ######################################################

    gmap_db_version = checkGmapDBVersion(gmap_db, reference_name)
    gmap_version = checkGmapVersion(gmap_exe)

    print(f"# Gmap version: {gmap_version})")
    print(f"# Gmap database version: {gmap_db_version} ({reference_name})\n")

    print(f"# minimum identity %: {min_identity}")
    print(f"# minimum coverage %: {min_coverage}\n")

    # parse input FASTA file and process sequences one by one,
    # note it would be more efficient to process them in batch to avoid
    # repeated calls to gmap
    names, seqs = parseFASTASeqs(fasta_file)
    for seqname in names:
    
        keys = []
        matched_coords = ""

        fasta_sequence = seqs[seqname]

        # define prefix for temp & output files out of this query
        tmp_filename_prefix = uuid.uuid4().hex

        temp_fasta_query = getTempFASTAFile(seqname, fasta_sequence, 
                                            tmp_filename_prefix, temp_path, 
                                            verbose=verbose_out)
    
        gmap_match, genome, reference_gff = runGmapAgainstGenomes(
                                               [reference_name], gmap_exe, gmap_db, 
                                                temp_fasta_query, min_identity, min_coverage,
                                                ncores, tmp_filename_prefix, temp_path, 
                                                verbose=verbose_out)
                                
        # MATCH IN REFERENCE, shortcut to avoid pangenome search
        if gmap_match == True:
            matched_coords, keys = checkRangesKeyInRef(
                                            reference_gff, hapIDranges, 
                                            f'{vcf_dbs}hvcf_files/', coverage=min_coverage,
                                            bedtools_path=bedtools_exe, grep_path=grep_exe, 
                                            verbose=verbose_out)

            # TODO: match gmap output format
            print(f"{seqname}\t{matched_coords}")
            if add_ranges:
                for key in keys:
                    print(keys[key])

            if keepTempFiles == False:
                os.remove(temp_fasta_query)
                os.remove(reference_gff)

        # NO MATCH IN REFERENCE, must scan pangenome hierarchically, slower
        elif gmap_match == False:
            if verbose_out == True:
                print(f"# No match found in reference genome ({reference_name}), look up the pangenome")

            pangenome_genomes = sortGenomesByRangeNumber(pangenome_fastas_folder, hapIDtable, verbose=verbose_out)    
        
            gmap_match, genome, pangenome_gff = runGmapAgainstGenomes(
                                                pangenome_genomes, gmap_exe, gmap_db, 
                                                temp_fasta_query, 
                                                min_identity, min_coverage,
                                                ncores, tmp_filename_prefix, temp_path,
                                                verbose=verbose_out)

            if gmap_match == False and verbose_out == True:
                print("# No matches in pangenome either\n")
            
            else:
                if verbose_out == True:
                    print(f"# Match found in {genome}")

            range_bedfile = f"{vcf_dbs}hvcf_files/{genome}.h.bed"    

            matched_coords, keys = checkRangesKeysInPangenome(
                                                pangenome_gff,hapIDranges,range_bedfile,
                                                f'{vcf_dbs}hvcf_files/',coverage=min_coverage,
                                                bedtools_path=bedtools_exe,grep_path=grep_exe,
                                                verbose=verbose_out)

            print(f"{seqname}\t{matched_coords}")
            if add_ranges:
                for key in keys:
                    print(keys[key])

            if keepTempFiles == False:
                os.remove(temp_fasta_query)
                os.remove(pangenome_gff)


# %%
if __name__ == "__main__":

    import argparse
    import subprocess
    import os
    from pathlib import Path
    import sys
    import re
    import tempfile
    import uuid
    import yaml

    main()


