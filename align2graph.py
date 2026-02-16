#!/usr/bin/env python3

# Aligns arbitrary DNA sequences to a gmap+PHG pangenome in up to 2 steps:
# i) hierarchically gmap input sequences to included genomes, starting from reference,
# ii) queries PHG haplotypes to obtain precomputed matches in other included genomes
# Returns reference-based coordinates of matched sequence
#
# J Sarria, B Contreras-Moreira
# Copyright [2024-26] Estacion Experimental de Aula Dei-CSIC

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
def check_gmap_version(gmap_path):
    """Returns version of GMAP binary/executable."""

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
            temp.write(b">" + name.encode() + b"\n" + seqs[name].encode() + b"\n") 
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
    pattern1 = r'coverage=([^;]+);identity=([^;]+)'
    pattern2 = r'identity=([^;]+);coverage=([^;]+)'
    pattern3 = r'Name=([^;]+)'

    if not os.path.isfile(gff_file):
        print(f"# ERROR(valid_matches): file {gff_file} does not exist")
        return matches
    
    else:
        with open (gff_file) as f:
            for line in f:
                if not line.startswith('#'):
                             
                    fields = line.split("\t")
                    
                    if fields[2] == "mRNA":
                    #gmap-2024-11-20
                    #ID=query.mrna1;Name=name;Parent=genome.path1;Dir=na;coverage=100.0;identity=98.9;
                    #gmap-2013-08-31
                    #ID=m84096...mrna1;Name=m84096..;Parent=m84096...path1;coverage=38.1;identity=99.3

                        match1 = re.search(pattern1, fields[-1])
                        if match1:
                            coverage = float(match1.group(1))
                            identity = float(match1.group(2))

                        else: # try other order
                            match2 = re.search(pattern2, fields[-1])
                            if match2: 
                                identity = float(match2.group(1))
                                coverage = float(match2.group(2)) 
                            else:
                                coverage = -1.0
                                identity = -1.0
                                if verbose == True: 
                                    print(f"# ERROR(valid_matches): cannot parse coverage/identity: {fields[-1]}")                               

                        if identity >=0 and identity >= min_identity and coverage >= min_coverage:
                            
                            match3 = re.search(pattern3, fields[-1])
                            if match3:
                                seqname = match3.group(1)
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
def get_rank(range_str, ranked_genomes):
    """
    Helper function to determine the rank of a range string based on a sorted genome list.
    range_str format is typically: chr2H@GDB_136:452737005-452738007(+)
    """
    if not ranked_genomes:
        return 0

    # extract genome: GDB_136
    match = re.search(r"@([^:]+):", range_str)
    if match:
        genome = match.group(1)
        try:
            return ranked_genomes.index(genome)
        except ValueError:
            return 99999    # Send it to the end if genome not found in the ranked list
    return 99999


# %%
def align_sequence_to_ranges(agc_path, agc_db_path, gmap_path, 
                             sequence, agc_ranges, ranked_genomes=None, 
                             aligner_tool='gmap', minimap2_path='minimap2', 
                             verbose=False):
    """
    Calls agc to cut genome fragments corresponding to list of ranges in appropriate agc format.
    Then aligns input sequence to these ranges using either GMAP (pairwise) or Minimap2 (splice).
    Returns alignment-corrected ranges in CSV string.
    """

    aligned_ranges  = []
    range_seqnames  = []
    range_sequences = {}
    name = ''
    seq  = ''

    if verbose:
        print(f"# INFO(align_sequence_to_ranges): agc_ranges: {agc_ranges} using {aligner_tool}")

    # --- Step 1: Cut sequences from AGC database and save in temp file ---
    command = f"{agc_path} getctg -p {agc_db_path} {' '.join(agc_ranges)}"
    try:
        result = subprocess.run(command,
                                shell=True, check=True, text=True, 
                                stdout=subprocess.PIPE,
                                stderr=subprocess.DEVNULL)

        for line in result.stdout.splitlines():
            header = re.search(r"^>", line)
            if header:
                name = line.replace('>','')
                name = name.replace(' ','_')
                range_seqnames.append(name)
                range_sequences[name] = ''
            else:
                range_sequences[name] = range_sequences[name] + line    

    except subprocess.CalledProcessError as e:
        print(f'# ERROR(align_sequence_to_ranges): agc failed, return unchecked ranges')
        return ";".join(agc_ranges)

    # --- Step 2: Align using selected tool ---
    if aligner_tool == 'minimap2':  # OPTION A: MINIMAP2 LOGIC
        
        query_temp_file_name = write_temp_fasta_file(['query'], {'query':sequence})

        for seqname in range_names:
            
            ref_temp_file_name = write_temp_fasta_file([seqname], {seqname:range_sequences[seqname]})

            # Run Minimap2
            command = [minimap2_path, "-ax", "splice", "--secondary=no", ref_temp_file_name, query_temp_file_name]
            
            try: 
                result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, check=True)
                
                if verbose:
                    print(f"# INFO(align_sequence_to_ranges): minimap2 finished for {seqname}")

                # Parse SAM Output
                for line in result.stdout.splitlines():
                    if line.startswith("@"): continue
                    
                    cols = line.split("\t")
                    if len(cols) < 6: continue
                    
                    flag = int(cols[1])
                    if flag & 4: continue # Unmapped
                    
                    ref_name = cols[2]
                    pos = int(cols[3])
                    cigar = cols[5]
                    
                    # Calculate alignment length on reference from CIGAR
                    ref_len = 0
                    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
                        if op in "MDN=X":
                            ref_len += int(length)
                    
                    # Parse Header and Calculate Global Coords
                    rangematch = re.search(r"(chr\d[a-zA-Z])_sampleName=([^:]+):(\d+)-(\d+)", seqname)
                    if rangematch:
                        contig = rangematch.group(1)
                        genome = rangematch.group(2)
                        chunk_start = int(rangematch.group(3))
                        
                        final_start = chunk_start + pos - 1
                        final_end = final_start + ref_len - 1
                        strand = '-' if (flag & 16) else '+'
                        
                        range_str = f'{contig}@{genome}:{final_start}-{final_end}({strand})'
                        aligned_ranges.append(range_str)

            except subprocess.CalledProcessError:
                print(f'# ERROR(align_sequence_to_ranges): minimap2 failed')
                os.remove(query_temp_file_name)
                os.remove(ref_temp_file_name)

            os.remove(ref_temp_file_name)
        
        os.remove(query_temp_file_name)

    else: # OPTION B: GMAP LOGIC

        for seqname in range_seqnames:
            range_str = '' 
            range_start = 0
            range_end = 0     
            range_strand = '+'

            temp_file_name = write_temp_fasta_file([seqname,'query'], 
                                                {seqname:range_sequences[seqname], 'query':sequence})

            command = f"cat {temp_file_name} | {gmap_path} -t 1 -n 1 -2"
            try: 
                result = subprocess.run(command,
                                    shell=True, check=True, text=True,
                                    stdout=subprocess.PIPE, 
                                    stderr=subprocess.DEVNULL)

                if verbose == True:
                    print(f"# INFO(align_sequence_to_ranges): gmap result for {seqname}:\n{result.stdout}")

                paths = re.search(r"Paths \([123456789]+", result.stdout)
                if paths:
                    rangematch = re.search(r"(chr\d[a-zA-Z])_sampleName=([^:]+):(\d+)-(\d+)", seqname)
                    if rangematch:
                        range_str = f'{rangematch.group(1)}@{rangematch.group(2)}'
                        range_start = int(rangematch.group(3))
                        range_end = int(rangematch.group(4))    

                    match = re.search(r"Genomic pos: ([\d,]+)..([\d,]+) \(([+-])", result.stdout)
                    if match:
                        gmap_start = int(match.group(1).replace(',',''))
                        gmap_end = int(match.group(2).replace(',',''))
                        gmap_strand = match.group(3)
                        if gmap_strand == '-':
                            range_strand = '-'
                            gmap_start, gmap_end = gmap_end, gmap_start 
                        
                        range_start += gmap_start - 1 
                        range_end = range_start + (gmap_end - gmap_start)
                        range_str = f'{range_str}:{range_start}-{range_end}({range_strand})'
                    
                    aligned_ranges.append(range_str)
                    os.remove(temp_file_name)

                else:
                    os.remove(temp_file_name)
                    if verbose == True:
                        print(f"# INFO(align_sequence_to_ranges): no valid alignment for {seqname}")
                        
            except subprocess.CalledProcessError as e:
                print(f'# ERROR(align_sequence_to_ranges): {e.cmd} failed: {e.stderr}')

    # Sort aligned ranges using the global get_rank function
    if ranked_genomes:
        aligned_ranges.sort(key=lambda x: get_rank(x, ranked_genomes))

    return ";".join(aligned_ranges)


# %%
def get_overlap_ranges_reference(gmap_match,hapIDranges,genomes,bed_folder_path,
                                coverage=0.75,all_graph_matches=False,
                                aligner_tool='gmap',
                                bedtools_path='bedtools',grep_path='grep',
                                agc_path='agc', agc_db_path='', gmap_path='gmap',
                                ranked_genomes=None,
                                verbose=False):
    """Retrieves PHG keys for ranges overlapping gmap match in reference genome.
    Passed coverage is used to intersect ranges and match. Overlap does not consider strandness. 
    Note: agc & gmap only used when all_graph_matches=True to confirm range matches.
    Returns: string with matched coords in TSV format.
    Column order in TSV: ref_chr, ref_start, ref_end, ref_strand (. if absent in ref), 
    match_genome, match_chr, match_start, match_end, match_strand, 
    match_identity,match_coverage,other_matches (Yes/No),graph_ranges."""

    keys = {}   # Allows storing multiple keys 
    match_tsv = ''
    mult_mappings = 'No'
    all_ranges = '.'

    chrom = gmap_match['chrom']
    genome = gmap_match['genome']
    start = gmap_match['start']
    end = gmap_match['end']
    strand = gmap_match['strand']
    ident = gmap_match['ident']
    cover = gmap_match['cover']
    if(gmap_match['matches'] > 1):
        mult_mappings = 'Yes'

    if verbose == True:
        print(f"# Checking match for {chrom}:{start}-{end} within reference")

    # prepare bedtools intersect command to find overlapping range, no strand check,
    # assume input coords are sorted to avoid getBin errors with chroms > 512MB,
    # see https://bedtools.readthedocs.io/en/latest
    command = f"{bedtools_path} intersect -a {hapIDranges} -b stdin -nonamecheck -e -F {coverage} -f {coverage} -sorted"             

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
            match_tsv = f'.\t.\t.\t.\t{genome}\t{chrom}\t{start}\t{end}\t{strand}\t{ident}\t{cover}\t{mult_mappings}\t'
            return match_tsv + all_ranges

        elif len(intersections) > 1:
            if(verbose == True):
                print(f"# WARN(get_overlap_ranges_reference): several overlaps, take 1st: {intersections}")
            intersections = intersections[:1]

        for feature in intersections:
            feature = str(feature).split("\t")
            feature[-1] = feature[-1].strip() 
            match_tsv = (
                f'{feature[0]}\t{feature[1]}\t{feature[2]}\t{strand}\t'
                f'{genome}\t{chrom}\t{start}\t{end}\t{strand}\t{ident}\t{cover}\t{mult_mappings}\t')
        
            if all_graph_matches == True:
                for c in range(0,len(genomes)):
                    k = feature[c+3]
                    if k == ".":
                        continue
                    else:
                        clean_k = k[1:-1] #remove <>
                        if clean_k not in keys:
                            keys[clean_k] = set()
                        keys[clean_k].add(genomes[c])

    except subprocess.CalledProcessError as e:
        print(f'# ERROR(get_overlap_ranges_reference): {e.cmd} failed: {e.stderr}')

    # retrieve genomic ranges matching these keys,
    # multiple matches allowed to support 1toN mappings (CNVs)
    if all_graph_matches == True:
        agc_ranges = []
        for k in keys:
            for genome in keys[k]:
                range_bedfile = f"{bed_folder_path}/{genome}.h.bed"
                command = f"{grep_path} {k} {range_bedfile}"
                if verbose == True:
                    print(f"# INFO(get_overlap_ranges_reference): grep for key {k} in {range_bedfile}")
                try:
                    result = subprocess.run(command,
                            shell=True,check=True,text=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

                    num_bed_lines = 0
                    for b in result.stdout.splitlines():
                        bed_data = b.split("\t")
                        if len(bed_data) > 4:
                            # Calculate 5% padding on each side. On genomic regions, UTRs from other genomes may be shorter than the sequence in other ranges, so adding padding on both sides to avoid cropping the end/start of a gene when add_ranges
                            start = int(bed_data[1])
                            end = int(bed_data[2])
                            block_length = end - start
                            padding = int(block_length * 0.05)
                            padding = max(padding, 1000)  # Use at least 500bp padding
                            # For short ranges 5% is not enough, so ensure always to use at least 500bp
                            padded_start = max(0, start - padding)  # Ensure non-negative
                            padded_end = end + padding
                            agc_range = f'{bed_data[0]}@{genome}:{padded_start}-{padded_end}' # agc format with 5% padding (>500bp)
                            agc_ranges.append(agc_range)
                            if verbose == True:
                                print(f"# INFO(get_overlap_ranges_reference): key {k} -> {bed_data[0]}:{start}-{end} + padding ({padding}bp) = {padded_start}-{padded_end}")
                        num_bed_lines += 1

                except subprocess.CalledProcessError as e:
                    print(f'# ERROR(get_overlap_ranges_reference): {e.cmd} failed: {e.stderr}')

        # Remove duplicates)
        agc_ranges = (set(agc_ranges))
        aligned_ranges = align_sequence_to_ranges(agc_path, agc_db_path, gmap_path,
                                                    gmap_match['sequence'], agc_ranges,
                                                    ranked_genomes=ranked_genomes,
                                                    aligner_tool=aligner_tool,
                                                    verbose=verbose)
        all_ranges = aligned_ranges 


    return match_tsv + all_ranges


# %%
def sort_genomes_by_range_number(hap_table_file, verbose=False):
    """Sorts the genomes in graph by number of ranges after parsing hap_table_file (max to min). 
    Returns list with genome names (without extension)."""

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
    
    # sort genomes
    if verbose == True:
        print("\n# genomes sorted by number of ranges:")

    pangenome_genomes = []
    for g in sorted(hapnum, key=hapnum.get, reverse=True):
        pangenome_genomes.append(g)
        if verbose == True:
            print(f'# {g} {hapnum[g]}')

    return pangenome_genomes


def genomes_from_graph(hapIDranges_filename):
    """Reads hapIDranges file and returns list of genomes in pangenome graph.
    Assumes that the 1st line starts with # and contains genome names."""

    genomes = []
    with open(hapIDranges_filename) as f:
        for line in f:
            if line.startswith('#'):
                genomes = line.split("\t")
                genomes = genomes[3:]  # skip first 3 columns
                genomes[-1] = genomes[-1].strip()  # remove newline
                break
    f.close()

    return genomes


# %%
def run_gmap_genomes(pangenome_genomes, gmap_path, gmap_db, fasta_filename, 
                         min_identity, min_coverage,
                         cores=4, prefix='temp', path='/tmp/', 
                         genomic=False,verbose=False):
    """Iteratively gmaps input FASTA file against list of genomes.
    Returns dictionary with GMAP matches with sequence names as 1ary keys.
    For each input sequence the following 2ary keys are created: 
    i) integer with number of matches (default 0),    
    ii) string with matched genome name (default ''), 
    iii to viii) strings with chromosome, start, end, strand, identity, cover
    """

    # regexes to parse gmap stderr
    pattern0 = r'Problem with sequence (\S+)'
    pattern1 = r'SIGSEGV'
    pattern2 = r'For big genomes of more than'
    pattern3 = r'This is a large genome of more than'

    gmap_matches = {}

    # parse sequences and init dictionary of matches
    seqnames, sequences = parse_fasta_file(fasta_filename)
    for seqname in seqnames:
        gmap_matches[seqname] = {}       
        gmap_matches[seqname]['sequence'] = ''
        gmap_matches[seqname]['matches'] = 0
        gmap_matches[seqname]['genome'] = ''
     
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
        if genomic == True:
            gmap_command = (
                f"{gmap_path} --nosplicing -D {gmap_db} -d {genome} -t {cores} {g_fasta_filename} -f gff3_gene > {g_gff_filename}") 
        if verbose == True:
            print(f'# {gmap_command}')

        try:
            result = subprocess.run(gmap_command, shell=True, check=True, 
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if verbose == True:
                print(result.stdout.decode())

        except subprocess.CalledProcessError as e:

            # long reads ie m84096_250523_032626_s1/247992835/ccs might fail unless splicing is disabled,
            # sometimes this warning is never printed, only SIGSEGV occurs
            match0 = re.search(pattern0, e.stderr.decode())
            if match0:
                print(f"# WARN: {match0.group(0)} ({genome}), consider using --genomic to disable splicing") 
            elif re.search(pattern1, e.stderr.decode()):
                print(f"# WARN: SIGSEGV detected ({genome}), consider using --genomic to disable splicing")

            else:
                # regexes to decide whether gmapl is needed (genomes > 4Gbp)
                match2 = re.search(pattern2, e.stderr.decode())
                match3 = re.search(pattern3, e.stderr.decode())

                # try gmapl instead
                if match2 or match3:
                    if verbose == True:
                        print(f"# Running gmapl for {genome}")

                    gmap_command = (
                        f"{gmap_path}l -D {gmap_db} -d {genome} -t {cores} {g_fasta_filename} -f gff3_gene > {g_gff_filename}")
                    if genomic == True:
                        gmap_command = (
                            f"{gmap_path}l --nosplicing -D {gmap_db} -d {genome} -t {cores} {g_fasta_filename} -f gff3_gene > {g_gff_filename}")
                    
                    try:
                        result = subprocess.run(gmap_command, shell=True, check=True, 
                                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        if verbose == True:
                            print(result.stdout.decode())

                    except subprocess.CalledProcessError as e:
                        print(f"\nERROR(run_gmap_genomes): '{e.cmd}' (gmapl) returned non-zero exit status {e.returncode}.")

                else:
                    print(f"\nERROR(run_gmap_genomes): '{e.cmd}' (gmap) returned non-zero exit status {e.returncode}.")

        genome_matches = valid_matches(g_gff_filename,min_identity,min_coverage,verbose=verbose)
        
        for seqname in genome_matches:
            gmap_matches[seqname]['genome'] = genome
            gmap_matches[seqname]['matches'] = genome_matches[seqname]['matches']
            gmap_matches[seqname]['chrom'] = genome_matches[seqname]['chrom']
            gmap_matches[seqname]['start'] = genome_matches[seqname]['start']
            gmap_matches[seqname]['end'] = genome_matches[seqname]['end']
            gmap_matches[seqname]['strand'] = genome_matches[seqname]['strand']
            gmap_matches[seqname]['ident'] = genome_matches[seqname]['ident']
            gmap_matches[seqname]['cover'] = genome_matches[seqname]['cover']
            gmap_matches[seqname]['sequence'] = sequences[seqname]

        # clean up temp files
        os.remove(g_gff_filename)
        os.remove(g_fasta_filename)
        
    return gmap_matches

# %%
def get_overlap_ranges_pangenome(gmap_match,hapIDranges,genomes,bedfile,bed_folder_path,
                                coverage=0.75,all_graph_matches=False,
                                aligner_tool='gmap',
                                bedtools_path='bedtools',grep_path='grep',
                                agc_path='agc', agc_db_path='', gmap_path='gmap',
                                ranked_genomes=None,
                                verbose=False):
    """Retrieves PHG keys for ranges overlapping gmap match in 1st matched pangenome assembly.
    BED file is usually a .h.bed file with sorted ranges extracted from PHG .h.vcf.gz files.
    Passed coverage is used to intersect ranges and match. Overlap does not consider strandness.
    Note: If a sequence is more than once mapped in a genome, only the 1st match is used.
    Note: agc & gmap only used when all_graph_matches=True to confirm range matches. # Now this should not be true
    Returns: string with matched coords in TSV format.
    Column order in TSV: ref_chr, ref_start, ref_end, ref_strand (. if absent in ref), 
    match_genome, match_chr, match_start, match_end, match_strand, 
    match_identity,match_coverage,other_matches (Yes/No),graph_ranges."""

    keys = {}
    match_tsv = ''
    graph_key = ''
    mult_mappings = 'No'
    aligned_ranges = '.'
    
    chrom = gmap_match['chrom']
    genome = gmap_match['genome']
    start = gmap_match['start']
    end = gmap_match['end']
    strand = gmap_match['strand']
    ident = gmap_match['ident']
    cover = gmap_match['cover']
    if(gmap_match['matches'] > 1):
        mult_mappings = 'Yes'

    if verbose == True:
        print(f"# Checking match for {chrom}:{start}-{end} at {bedfile}")

    # prepare bedtools intersect command to find overlapping range, no strand check,
    # see also https://github.com/arq5x/bedtools2/issues/679,
    # bedfile should contain lines like this:
    # chr1H	975076	975705	+	3e9613a955ec342330fd2138f652fddc	HOR_13942	chr1H	76744	77373	35451ad1b65dd5ff4d940f914c21dbe1
    command = f"{bedtools_path} intersect -a {bedfile} -b stdin -nonamecheck -e -F {coverage} -f {coverage} -sorted "

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
            match_tsv = f'.\t.\t.\t.\t{genome}\t{chrom}\t{start}\t{end}\t{strand}\t{ident}\t{cover}\t{mult_mappings}\t'
            return match_tsv + aligned_ranges

        elif len(intersections) > 1:
            if(verbose == True):
                print(f"# WARN(get_overlap_ranges_pangenome): > several overlaps, take 1st: {intersections}")
            intersections = intersections[:1]
                    
        for feature in intersections:
            feature = str(feature).split("\t")
            feature[-1] = feature[-1].replace("\n", "") 
            match_tsv = ( # strand unknown as match is missing from reference genome
                f'{feature[6]}\t{feature[7]}\t{feature[8]}\t.\t'
                f'{genome}\t{chrom}\t{start}\t{end}\t{strand}\t{ident}\t{cover}\t{mult_mappings}\t')
            graph_key = feature[4]

        if verbose == True:
            print(f"# INFO(get_overlap_ranges_pangenome): graph key: {graph_key}")

        # look for this key within graph ranges (grep)
        if all_graph_matches == True:
            command = f"{grep_path} {graph_key} {hapIDranges}"
            try:
                result = subprocess.run(command,
                                shell=True,check=True,text=True,input=match_interval,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE)

                graph_data = result.stdout.splitlines()

                if len(graph_data) > 0:

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

                else:
                    #if verbose == True:
                    print(f"# WARN(get_overlap_ranges_pangenome): no graph matches for {gmap_match} {graph_key}")

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

                num_bed_lines = 0
                for b in result.stdout.splitlines():
                    bed_data = b.split("\t")
                    if len(bed_data) > 4:
                        # Calculate 5% padding on each side
                        # At some cases, matches from other genomes may be shorter than the sequence in other ranges, so adding padding on 
                        # both sides to avoid cropping the end/start of a gene when add_ranges
                        start = int(bed_data[1])
                        end = int(bed_data[2])
                        block_length = end - start
                        padding = int(block_length * 0.05)
                        padded_start = max(0, start - padding)  # Ensure non-negative
                        padded_end = end + padding
                        padded_range = f'{bed_data[0]}@{keys[k]}:{padded_start}-{padded_end}' # agc format with 5% padding
                        if num_bed_lines == 0:
                            keys[k] = padded_range
                        else:
                            keys[k] += f' {padded_range}' # with a space separator to distinguish multiple ranges
                    num_bed_lines += 1

                    if verbose == True:
                        print(f"# INFO(get_overlap_ranges_pangenome): found range for key {k}: {bed_data[0]}:{start}-{end} with padding {padded_start}-{padded_end}")

            except subprocess.CalledProcessError as e:
                print(f'# ERROR(get_overlap_ranges_pangenome): {e.cmd} failed: {e.stderr}')
        
        #aligned_ranges = ';'.join(sorted(keys.values())) # strand unknown, unchecked
        aligned_ranges = align_sequence_to_ranges(agc_path, agc_db_path, gmap_path,
                                                    gmap_match['sequence'], sorted(keys.values()),
                                                    ranked_genomes=ranked_genomes,
                                                    aligner_tool=aligner_tool,
                                                    verbose=verbose)

    return match_tsv + aligned_ranges


# %%
def main():

    parser = argparse.ArgumentParser(
        description="Map sequences within pangenome graph.\n",
        epilog="Citation: see https://github.com/eead-csic-compbio/barleygraph\n",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "fasta_file",
        help="path to FASTA file with sequences to map (required)"
    )

    parser.add_argument(
        "--graph_yaml",
        default="graph.yaml",
        help="path to YAML file with pangenome config, dafault: graph.yaml"
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
        "--agc_exe",
        default='agc',
        help=f"path to agc executable, default: agc"
    )

    parser.add_argument(
        "--cor", 
        default=4, 
        help="number of cores for gmap, default: 4"
    )

    parser.add_argument(
        "--minident",
        default=98.0,
        help="min %%identity of gmap matches, default: 98.0"
    )

    parser.add_argument(
        "--mincover",
        default=95.0,
        help="min %%coverage of gmap matches, default: 95.0"
    )

    parser.add_argument(
        "--mincover_range",
        default=75.0,
        help="min %%coverage of gmap matches and pangenome ranges, default: 75.0"
    ) 

    parser.add_argument(
        "--single_genome",
        default='',
        help="selected genome to be scanned with GMAP, must be part of graph, default: all genomes. Note that --add_ranges may not work properly with this option.",
    )   # If used together with add_ranges a warning is printed

    parser.add_argument(
        '--verb', 
        action='store_true', 
        help='Increase verbosity in output'
    )
    
    parser.add_argument(
        '--genomic', 
        action='store_true', 
        help='Input sequences are genomic, turn off splicing'
    )
    
    # MODIFIED: Now is not a flag, but a choice of tool
    parser.add_argument('--add_ranges', 
        type=str,
        choices=['gmap', 'minimap2'],
        help='Add all pangenome ranges matching input sequences using specified tool (gmap or minimap2)'
    )

    args = parser.parse_args()

    # Sequential execution for now, parallelization could be implemented later if needed
    process_sequences_serial(args)


def process_sequences_serial(args):
    """Main processing logic for sequences (called by main or parallel workers)"""
    
    grep_exe = 'grep' # supposed to be available in PATH
    
    # required params
    fasta_file = args.fasta_file

    # parse graph YAML file
    chr_syns = {}
    try:
        f = open(args.graph_yaml, "r")
    except FileNotFoundError:
        print(f'# ERROR: cannot read file {args.graph_yaml}, please correct --graph_yaml')
        sys.exit(-1)
    else:
        with f:        
            config = yaml.load(f, Loader=yaml.FullLoader)

            hvcf_bed = config['hvcf_bed']
            agc_db = config['agc_assemblies']            
            gmap_db = config['gmap_db']
            reference_name = config['reference_name']
            hapIDtable  = config['hapIDtable']
            hapIDranges = config['hapIDranges']
            gmap_exe = config['gmap_exe']

            if 'chr_syns' in config:
                # if chr_syns is present, it is a dictionary with chromosome synonyms
                # e.g. {'chr1H': 'chr1H_LR890096.1', 'chr2H': 'chr2H_LR890097.1', ...}
                chr_syns = config['chr_syns']

    # get other optional params
    bedtools_exe  = args.bedtools_exe
    agc_exe       = args.agc_exe
    ncores        = int(args.cor)
    min_identity  = float(args.minident)
    min_coverage  = float(args.mincover)
    min_coverage_range = float(args.mincover_range)
    temp_path     = args.tmp_path
    verbose_out   = args.verb
    genomic       = args.genomic
    
    # Logic for add_ranges
    aligner_tool  = args.add_ranges
    do_add_ranges = (aligner_tool is not None)
    
    single_genome = args.single_genome

    ######################################################

    gmap_version = check_gmap_version(gmap_exe)

    print(f"# GMAP version: {gmap_version}")
    print(f"# config_file: {args.graph_yaml}")
    print(f"# fasta_file: {fasta_file}")
    print(f"# minimum identity %: {min_identity}")
    print(f"# minimum coverage %: {min_coverage}")
    print(f"# minimum coverage range %: {min_coverage_range}")
    print(f"# add_ranges: {do_add_ranges} (Tool: {aligner_tool if do_add_ranges else 'None'})")
    print(f"# genomic: {genomic}")

    if single_genome != '' and do_add_ranges:
        print(f"# WARNING: --add_ranges may not work properly when --single_genome is used")
        print(f"# because only ranges with exactly same haplotype are considered, not surfing through whole pangenome\n")

    # order of genes to be hierarchically scanned with GMAP
    ranked_pangenome_genomes = sort_genomes_by_range_number(
        hapIDtable, 
        verbose=verbose_out) 

    if single_genome == '':
        print(f"# ranked pangenome genomes: {', '.join(ranked_pangenome_genomes)}\n")    
    elif single_genome in ranked_pangenome_genomes: 
        ranked_pangenome_genomes = [ single_genome ]
        print(f"# single genome to be scanned with GMAP: {single_genome}\n") 

    # check GMAP indices are in place
    for genome in ranked_pangenome_genomes:
        gmap_index_file = f"{gmap_db}/{genome}/{genome}.ref153positions"        
        if not os.path.isfile(gmap_index_file):
            print(f"# ERROR: missing GMAP index for {genome}, please index")
            sys.exit(-2)

    graph_pangenome_genomes = genomes_from_graph(hapIDranges)   
    
    # define prefix for temp & output files
    temp_prefix = uuid.uuid4().hex
        
    # match all sequences in one batch    
    gmap_matches = run_gmap_genomes(
        ranked_pangenome_genomes, 
        gmap_exe, gmap_db, 
        fasta_file, 
        min_identity, min_coverage,
        ncores, 
        temp_prefix, temp_path,
        genomic,
        verbose=verbose_out)
    

    # print header
    print(f'#query\tref_chr\t\tref_start\tref_end\t\tref_strand\tgenome\t'
          'chr\tstart\tend\tstrand\tperc_ident\tperc_cover\tmultmaps\tgraph_ranges')
    
    # compute graph coordinates for matched sequences
    for seqname in gmap_matches:
        if gmap_matches[seqname]['matches'] > 0:

            if(gmap_matches[seqname]['genome'] == reference_name):

                matched_coords = get_overlap_ranges_reference(
                    gmap_matches[seqname], 
                    hapIDranges,
                    graph_pangenome_genomes, 
                    hvcf_bed, # --- CHANGED --- 
                    coverage=min_coverage_range/100,
                    all_graph_matches=do_add_ranges, 
                    aligner_tool=aligner_tool,       
                    bedtools_path=bedtools_exe, 
                    grep_path=grep_exe, 
                    agc_path=agc_exe,
                    agc_db_path=agc_db,
                    gmap_path=gmap_exe,
                    ranked_genomes=ranked_pangenome_genomes,
                    verbose=verbose_out)
                
            else:  
                matched_coords = get_overlap_ranges_pangenome(
                    gmap_matches[seqname],
                    hapIDranges,
                    graph_pangenome_genomes,
                    f"{hvcf_bed}/{gmap_matches[seqname]['genome']}.h.bed", # --- CHANGED ---
                    hvcf_bed, # --- CHANGED ---
                    coverage=min_coverage_range/100,
                    all_graph_matches=do_add_ranges, 
                    aligner_tool=aligner_tool,       
                    bedtools_path=bedtools_exe,
                    grep_path=grep_exe,
                    agc_path=agc_exe,
                    agc_db_path=agc_db,
                    gmap_path=gmap_exe,
                    ranked_genomes=ranked_pangenome_genomes, 
                    verbose=verbose_out)

            # print output coordinates and update chr names if needed
            outTSV = matched_coords.split("\t")
            matched_coords = "\t".join(outTSV[1:])

            if outTSV[0] in chr_syns:
                ref_chr = chr_syns[outTSV[0]]
            else:
                ref_chr = outTSV[0]

            print(f"{seqname}\t{ref_chr}\t{matched_coords}")


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
    import time
    import yaml


    start_time = time.time()
    main()
