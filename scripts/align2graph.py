#!/usr/bin/env python3

# Aligns arbitrary DNA sequences to a gmap+PHG pangenome in up to 2 steps:
# i) hierarchically gmap input sequences to included genomes, starting from reference,
# ii) queries PHG haplotypes to obtain precomputed matches in other included genomes
# Returns reference-based coordinates of matched sequence
#
# J Sarria, B Contreras-Moreira
# Copyright [2024-25] Estacion Experimental de Aula Dei-CSIC

# %%
def parse_fasta_file(fasta, verbose=False):
    """Takes a FASTA filename and parses sequence names before 1st space.
    Returns: i) list of sequence names in input order, 
    ii) dictionary with sequence names as keys"""

    if verbose == True:
        print("# INFO(parse_fasta_names) parsing FASTA file:", fasta)

    names = []
    sequences = {}
    name = ''

    try:
        file = open(fasta)
    except OSError as error:
        print("# ERROR(parse_fasta_names): cannot open/read file:", fasta, error)
        return [],{}

    for line in file:
        header = re.search(r"^>", line)
        if header:
            seq_name_match = re.search(r"^>(\S+)", line)
            if seq_name_match:
                name = seq_name_match.group(1)
                names.append(name)
            else:
                print("# ERROR(parse_fasta_names): cannot parse header:", header)
        else:
            if name in sequences:
                sequences[name] = sequences[name] + line
            else:
                sequences[name] = line
                
    return names, sequences

# %%
def check_gmapdb_version(gmap_db, ref_name):
    """Returns version of reference Gmap db, '?' by default.
    Note that old gmap versions do not include version string in file.version."""

    version_db = '?'

    gmap_version_file=f'{gmap_db}/{ref_name}/{ref_name}.version'

    if not os.path.isfile(gmap_version_file):
        print(f"# ERROR(check_gmapdb_version: file {gmap_version_file} does not exist")
        return version_db

    else:
        with open (gmap_version_file) as f:
            for line in f:
                if not line.startswith(ref_name):
                    version_db = line.strip()

    return version_db

# %%
def check_gmap_version(gmap_path):
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
        version_exe = data[2].split()[5]

    except subprocess.CalledProcessError as e:
        print(f'# ERROR(check_gmap_version): {e.cmd} failed: {e.stderr}')

    return version_exe

# %%
def write_temp_fasta_file(names, seqs, prefix_string="temp",
                     temp_path="/tmp/", suffix_string=".fna",
                     verbose=False):
    """Takes up to 5 strings + 1 Boolean: 
    i) list if sequence names, ii) dictionary of sequences, 
    iii) temp prefix, iv) temp path, v) temp suffix/extension, 
    vi) verbose
    Creates a temporary FASTA file with passed sequences.
    Returns path to created temp file."""
   
    with tempfile.NamedTemporaryFile(prefix=prefix_string,suffix=suffix_string,
                                     dir=temp_path,delete=False) as temp:
        for name in names:
            temp.write(b">" + name.encode() + b"\n" + seqs[name].encode()) 
        temp.close()

    file_path = Path(temp.name)
    if file_path.exists() and file_path.stat().st_size > 0:
        if verbose == True:
            print(f"# INFO(get_temp_fasta_file) created {temp.name}")
        return temp.name
    else:
        print("# ERROR(get_temp_fasta_file): cannot write temp file")
        return ''

# %%
def valid_matches(gff_file, min_identity, min_coverage, verbose=False):
    """Checks if GFF3 contains gene & mRNA features with required identity and coverage.
    Assumes 1st match is best, but counts all matches satisfying identity and coverage.
    Returns dictionary with sequence names as 1ary key and the following 2ary keys:
    i) matches (int), 
    ii) ident (float), 
    iii) cover (float), 
    iv) chrom (string),
    v) start (int),
    vi) end (int),
    vii) strand (string)"""

    matches = {}

    if not os.path.isfile(gff_file):
        print(f"# ERROR(valid_matches): file {gff_file} does not exist")
        return matches
    
    else:
        with open (gff_file) as f:
            for line in f:
                if not line.startswith('#'):
                    fields = line.split("\t")
                    attributes = fields[-1].split(";")
                    seqname = attributes[1].split("=")[1]
                    
                    #gmap-2024-11-20
                    #ID=query.mrna1;Name=name;Parent=genome.path1;Dir=na;coverage=100.0;identity=98.9;
                    #gmap-2013-08-31
                    #ID=chr1start.mrna4;Name=chr1start;Parent=chr1start.path4;coverage=99.5;identity=82.4
                    if fields[2] == "mRNA":
                        
                        if(attributes[3].startswith("Dir=")):    
                            coverage = float(attributes[4].split("=")[1])
                            identity = float(attributes[5].split("=")[1])   
                        else:
                            identity = float(attributes[3].split("=")[1])   
                            coverage = float(attributes[4].split("=")[1])

                        if identity >= min_identity and coverage >= min_coverage:
                            if seqname not in matches:
                                matches[seqname] = {}
                                matches[seqname]['matches'] = 0

                            matches[seqname]['matches'] = matches[seqname]['matches'] + 1 

                            if(matches[seqname]['matches'] == 1):
                                matches[seqname]['ident'] = identity
                                matches[seqname]['cover'] = coverage
                                matches[seqname]['chrom'] = fields[0]
                                matches[seqname]['start'] = int(fields[3])
                                matches[seqname]['end'] = int(fields[4])
                                matches[seqname]['strand'] = fields[6]
                                if verbose == True:                    
                                    print("#",line)

    return matches


# %%
def get_overlap_ranges_reference(gmap_match,hapIDranges,bed_folder_path,
                                coverage=0.75,all_graph_matches=False,
                                bedtools_path='bedtools',grep_path='grep',
                                verbose=False):
    """Retrieves PHG keys for ranges overlapping gmap match in reference genome.
    Passed coverage is used to intersect ranges and match. Overlap does not consider strandness. 
    Returns: i) string with matched coords in TSV format, ii) dictionary of keys/genomic ranges
    of all matched genomes (useful to compute alignments down the road).
    Column order in TSV: graph_chr, graph_start, graph_end, graph_strand (. if absent in ref), 
    multiple_mappings (Yes/No), match_genome, match_chr, match_start, match_end, match_strand, other_matches"""

    keys = {}
    match_tsv = ''
    mult_mappings = 'No'
    all_ranges = '.'

    chrom = gmap_match['chrom']
    genome = gmap_match['genome']
    start = gmap_match['start']
    end = gmap_match['end']
    strand = gmap_match['strand']
    if(gmap_match['matches'] > 1):
        mult_mappings = 'Yes'

    if verbose == True:
        print(f"# Checking match for {chrom}:{start}-{end} within reference")

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
    command = f"{bedtools_path} intersect -a {hapIDranges} -b stdin -nonamecheck -e -F {coverage} -f {coverage}"             

    # BED-format interval of gmap match
    match_interval = f'{chrom}\t{start}\t{end}'

    # actually run the bedtools command   
    try:
        result = subprocess.run(command,
                                shell=True,check=True,text=True,input=match_interval,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                        
        intersections = result.stdout.splitlines()
        if verbose == True:
            print(f"# INFO(get_overlap_ranges_reference): {result.stdout}")

        if(len(intersections) == 0):
            match_tsv = f'{chrom}\t{start}\t{end}\t.\t{mult_mappings}\t{genome}\t.\t.\t.\t.\t'
            return match_tsv + all_ranges

        elif len(intersections) > 1:
            if(verbose == True):
                print(f"# WARN(get_overlap_ranges_reference): several overlaps, take 1st: {intersections}")
            intersections = intersections[:1]

        for feature in intersections:
            feature = str(feature).split("\t")
            feature[-1] = feature[-1].strip() 
            match_tsv = (
                f'{feature[0]}\t{feature[1]}\t{feature[2]}\t{strand}\t{mult_mappings}'
                f'\t{genome}\t{feature[0]}\t{feature[1]}\t{feature[2]}\t{strand}\t')
            
            if all_graph_matches == True:
                for c in range(0,len(genomes)):
                    k = feature[c+3]
                    if k == ".":
                        continue
                    else:
                        clean_k = k[1:-1] #remove <>
                        if clean_k not in keys:
                            keys[clean_k] = genomes[c]

    except subprocess.CalledProcessError as e:
        print(f'# ERROR(get_overlap_ranges_reference): {e.cmd} failed: {e.stderr}')

    # retrieve genomic ranges matching these keys
    if all_graph_matches == True:
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
                    bed_data = bed_data[:1]
                    if(verbose == True):
                        print(f"# WARN(get_overlap_range_references): several key matches, take 1st: {result.stdout}")
                    
                bed_data = bed_data[0].split("\t")
                keys[k] = f'{keys[k]}:{bed_data[0]}:{bed_data[1]}-{bed_data[2]}({bed_data[3]})'

            except subprocess.CalledProcessError as e:
                print(f'# ERROR(get_overlap_ranges_reference): {e.cmd} failed: {e.stderr}')
        all_ranges = ";".join(sorted(keys.values()))

    return match_tsv + all_ranges


# %%
def sort_genomes_by_range_number(folder_path, hap_table_file, verbose=False):
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
    if verbose == True:
        print("\n# genomes sorted by number of ranges:")

    pangenome_genomes = []
    for g in sorted(hapnum, key=hapnum.get, reverse=True):
        if(os.path.isfile(f'{folder_path}/{g}.fa')):
            pangenome_genomes.append(g)
            if verbose == True:
                print(f'# {g} {hapnum[g]}')

    #pangenome_genomes = ['HOR_13942'] #debug

    return pangenome_genomes

# %%
def run_gmap_genomes(pangenome_genomes, gmap_path, gmap_db, fasta_filename, 
                         min_identity, min_coverage,
                         cores=4, prefix='temp', path='/tmp/', verbose=False):
    """Iteratively gmaps input FASTA file against list of genomes.
    Returns dictionary with GMAP matches with sequence names as 1ary keys.
    For each input sequence the following 2ary keys are created: 
    i) integer with number of matches,    
    ii) string with matched genome name, 
    iii) string with GFF3 of match"""
    
    gmap_matches = {}

    # parse sequences and init dictionary of matches
    seqnames, sequences = parse_fasta_file(fasta_filename)
    for seqname in seqnames:
        gmap_matches[seqname] = {}        
        gmap_matches[seqname]['matches'] = 0
        gmap_matches[seqname]['genome'] = ''
        gmap_matches[seqname]['gff3'] = ''
     
    # loop over genomes hierarchically
    for genome in pangenome_genomes:
    
        # check which sequences still need to be gmapped
        g_seqnames = []
        g_sequences = {}
        for seqname in seqnames:
            if gmap_matches[seqname]['matches'] == 0:
                g_seqnames.append(seqname)
                g_sequences[seqname] = sequences[seqname]    

        if len(g_seqnames) == 0:
            if verbose == True:
                print(f"# INFO(run_gmap_genomes): all sequences already mapped\n")
            break
        else:
            if verbose == True:
                print(f"\n# INFO(run_gmap_genomes): gmap {len(g_seqnames)} sequences against {genome}")

        # create temp FASTA file with sequences to gmap        
        g_prefix = uuid.uuid4().hex
        g_fasta_filename = write_temp_fasta_file(g_seqnames, g_sequences, g_prefix, path) 

        g_gff_filename = f"{path}/{genome}.{g_prefix}.gff3"
        
        # try default gmap 
        gmap_command = (
            f"{gmap_path} -D {gmap_db} -d {genome} -t {cores} {g_fasta_filename} -f gff3_gene > {g_gff_filename}")
        if verbose == True:
            print(gmap_command)

        try:
            result = subprocess.run(gmap_command, shell=True, check=True, 
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

                gmap_command = (
                    f"{gmap_path}l -D {gmap_db} -d {genome} -t {cores} {g_fasta_filename} -f gff3_gene > {g_gff_filename}")

                try:
                    result = subprocess.run(gmap_command, shell=True, check=True, 
                                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    if verbose == True:
                        print(result.stdout.decode())

                except subprocess.CalledProcessError as e:
                    print(f"\nERROR(run_gmap_genomes): '{e.cmd}' returned non-zero exit status {e.returncode}.")

            else:    
                print(f"\nERROR(run_gmap_genomes): '{e.cmd}' returned non-zero exit status {e.returncode}.")

        genome_matches = valid_matches(g_gff_filename,min_identity,min_coverage,verbose=verbose)
        
        for seqname in genome_matches:
            gmap_matches[seqname]['genome'] = genome
            gmap_matches[seqname]['matches'] = genome_matches[seqname]['matches']
            gmap_matches[seqname]['chrom'] = genome_matches[seqname]['chrom']
            gmap_matches[seqname]['start'] = genome_matches[seqname]['start']
            gmap_matches[seqname]['end'] = genome_matches[seqname]['end']
            gmap_matches[seqname]['strand'] = genome_matches[seqname]['strand']

        # clean up temp files
        os.remove(g_gff_filename)
        os.remove(g_fasta_filename)
        
    return gmap_matches

# %%
def get_overlap_ranges_pangenome(gmap_match,hapIDranges,bedfile,bed_folder_path,
                                coverage=0.75,all_graph_matches=False,
                                bedtools_path='bedtools',grep_path='grep',
                                verbose=False):
    """Retrieves PHG keys for ranges overlapping gmap match in 1st matched pangenome assembly.
    BED file is usually a .h.bed file with sorted ranges extracted from PHG .h.vcf.gz files.
    Passed coverage is used to intersect ranges and match. Overlap does not consider strandness.
    Returns: i) string with matched coords in TSV format, ii) dictionary of keys/genomic ranges
    of all matched genomes (useful to compute alignments down the road).
    Column order in TSV: graph_chr, graph_start, graph_end, graph_strand (. if absent in ref),
    multiple_mappings (Yes/No), match_genome, match_chr, match_start, match_end, match_strand, other_matches"""

    keys = {}
    match_tsv = ''
    graph_key = ''
    mult_mappings = 'No'
    all_ranges = '.'
    
    chrom = gmap_match['chrom']
    genome = gmap_match['genome']
    start = gmap_match['start']
    end = gmap_match['end']
    strand = gmap_match['strand']
    if(gmap_match['matches'] > 1):
        mult_mappings = 'Yes'

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
    command = f"{bedtools_path} intersect -a {bedfile} -b stdin -nonamecheck -e -F {coverage} -f {coverage}"             

    # BED-format interval of gmap match
    match_interval = f'{chrom}\t{start}\t{end}'

    # actually call bedtools
    try:
        result = subprocess.run(command,
                                shell=True,check=True,text=True,input=match_interval,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                        
        intersections = result.stdout.splitlines()
        if verbose == True:
            print(f"# INFO(get_overlap_ranges_pangenome): {result.stdout}")


        if(len(intersections) == 0):
            match_tsv = f'{chrom}\t{start}\t{end}\t.\t{mult_mappings}\t{genome}\t.\t.\t.\t.\t'
            return match_tsv + all_ranges

        elif len(intersections) > 1:
            if(verbose == True):
                print(f"# WARN(get_overlap_ranges_pangenome): > several overlaps, take 1st: {intersections}")
            intersections = intersections[:1]
                    
        for feature in intersections:
            feature = str(feature).split("\t")
            feature[-1] = feature[-1].replace("\n", "") 
            match_tsv = ( # strand unknown as match is missing from reference genome
                f'{feature[6]}\t{feature[7]}\t{feature[8]}\t.\t{mult_mappings}'
                f'\t{genome}\t{feature[0]}\t{feature[1]}\t{feature[2]}\t{feature[3]}\t')
            graph_key = feature[4]

        # look for this key within graph ranges (grep)
        if all_graph_matches == True:
            command = f"{grep_path} {graph_key} {hapIDranges}"
            try:
                result = subprocess.run(command,
                                shell=True,check=True,text=True,input=match_interval,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE)

                graph_data = result.stdout.splitlines()
                if len(graph_data) > 1:
                    if(verbose == True):
                        print(f"# WARN(get_overlap_ranges_pangenome): several graph matches: {result.stdout}")
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
                print(f'# ERROR(get_overlap_ranges_pangenome): {e.cmd} failed: {e.stderr}')
                            
    except subprocess.CalledProcessError as e:
        print(f'ERROR(get_overlap_ranges_pangenome): {e.cmd} failed: {e.stderr}')

    # retrieve genomic ranges matching these keys
    if all_graph_matches == True:
        for k in keys:
            range_bedfile = f"{bed_folder_path}/{keys[k]}.h.bed"
            command = f"{grep_path} {k} {range_bedfile}"
            try:
                result = subprocess.run(command,
                                shell=True,check=True,text=True,input=match_interval,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE)

                bed_data = result.stdout.splitlines()
                if len(bed_data) > 1:
                    bed_data = bed_data[:1]
                    if(verbose == True):
                        print(f"# WARN(get_overlap_ranges_pangenome): several key matches, take 1st: {result.stdout}")
                    
                bed_data = bed_data[0].split("\t")
                keys[k] = f'{keys[k]}:{bed_data[0]}:{bed_data[1]}-{bed_data[2]}({bed_data[3]})'

            except subprocess.CalledProcessError as e:
                print(f'# ERROR(get_overlap_ranges_pangenome): {e.cmd} failed: {e.stderr}')
        
        all_ranges = ";".join(sorted(keys.values()))

    return match_tsv + all_ranges


# %%
def main():

    # hard-coded defaults
    grep_exe = 'grep' # assumed to be available

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
        help="min%coverage of gmap matches, default: 95.0",
    )

    parser.add_argument(
        "--mincover_range",
        default=75.0,
        help="min%coverage of gmap matches and pangenome ranges, default: 75.0",
    )

    #guess this is done from bmap_align_to_graph?
    #parser.add_argument('--show-unmapped', action='store_true', help='Show unmapped sequences in output')

    parser.add_argument(
        '--verb', 
        action='store_true', 
        help='Increase verbosity in output')
    
    parser.add_argument('--add_ranges', 
        action='store_true', 
        help='Add all pangenome ranges matching input sequences')

    args = parser.parse_args()

    fasta_file = args.fasta_file

    # parse YAML pangenome config file
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
    ncores       = int(args.cor)
    min_identity = float(args.minident)
    min_coverage = float(args.mincover)
    min_coverage_range = float(args.mincover_range)
    temp_path    = args.tmp_path
    verbose_out  = args.verb
    add_ranges   = args.add_ranges


    ######################################################

    gmap_db_version = check_gmapdb_version(gmap_db, reference_name)
    gmap_version = check_gmap_version(gmap_exe)

    print(f"# Gmap version: {gmap_version}")
    if(gmap_db_version != '?'):
        print(f"# Gmap database version: {gmap_db_version} ({reference_name})\n")
    else:
        print(f"# Gmap database version: unknown\n")

    print(f"# minimum identity %: {min_identity}")
    print(f"# minimum coverage %: {min_coverage}")
    print(f"# minimum coverage range %: {min_coverage_range}\n")

    if verbose_out == True:
        print(f"# verbose: {verbose_out}")
    
    pangenome_genomes = sort_genomes_by_range_number(
        pangenome_fastas_folder, 
        hapIDtable, 
        verbose=verbose_out) 
    
    # define prefix for temp & output files
    temp_prefix = uuid.uuid4().hex
        
    # match all sequences in one batch    
    gmap_matches = run_gmap_genomes(
        pangenome_genomes, 
        gmap_exe, gmap_db, 
        fasta_file, 
        min_identity, min_coverage,
        ncores, 
        temp_prefix, temp_path,
        verbose=verbose_out)
    

    # print header
    print(f'sequence\tchr\t\tstart\tend\t\tstrand\tmult_mappings\tgenome\t'
          'graph_chr\tgraph_start\tgraph_end\tgraph_strand\tgraph_ranges')
    
    # compute graph coordinates for matched sequences
    for seqname in gmap_matches:
        if gmap_matches[seqname]['matches'] > 0:

            if(gmap_matches[seqname]['genome'] == reference_name):

                matched_coords = get_overlap_ranges_reference(
                    gmap_matches[seqname], 
                    hapIDranges, 
                    f'{vcf_dbs}hvcf_files/', 
                    coverage=min_coverage_range/100,
                    all_graph_matches=add_ranges,
                    bedtools_path=bedtools_exe, 
                    grep_path=grep_exe, 
                    verbose=verbose_out)
                
            else:  
                matched_coords = get_overlap_ranges_pangenome(
                    gmap_matches[seqname],
                    hapIDranges,
                    f"{vcf_dbs}hvcf_files/{gmap_matches[seqname]['genome']}.h.bed",
                    f'{vcf_dbs}hvcf_files/',
                    coverage=min_coverage_range/100,
                    all_graph_matches=add_ranges,
                    bedtools_path=bedtools_exe,
                    grep_path=grep_exe,
                    verbose=verbose_out)

            # print output coordinates    
            print(f"{seqname}\t{matched_coords}")

    # else: 
    #     if args.show_unmapped == True:
    #         for seqname in gmap_matches:
    #             print(f"# {seqname} not mapped")

         


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



