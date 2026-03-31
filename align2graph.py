#!/usr/bin/env python3

# Aligns arbitrary DNA sequences to a gmap+PHG pangenome in up to 2 steps:
# i) hierarchically gmap input sequences to included genomes, starting from reference,
# ii) queries PHG haplotypes to obtain precomputed matches in other included genomes
# Returns 1-based reference-based coordinates of matched sequence
#
# J Sarria, B Contreras-Moreira
# Copyright [2024-26] Estacion Experimental de Aula Dei-CSIC

# %%
MINIMAP2_N_DEFAULT = 20


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
def valid_matches(gff_file, genome_name, query_sequences, min_identity, min_coverage, verbose=False):
    """Checks if GFF3 contains gene & mRNA features with required identity and coverage.
    Assumes 1st match is best, but counts all matches satisfying identity and coverage.
    Returns dictionary of lists with sequence names as 1ary key and the following 2ary keys:
    i) genome (string),
    ii) sequences (dictionary),
    iii) ident (float), 
    iv) cover (float), 
    v) chrom (string),
    vi) start (int),
    vii) end (int),
    viii) strand (string)"""

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
                                    matches[seqname] = []
                                    
                                matches[seqname].append({
                                    'genome':genome_name,
                                    'sequence':query_sequences.get(seqname, ''),
                                    'ident':identity,
                                    'cover':coverage,
                                    'chrom':fields[0],
                                    'start':int(fields[3]),
                                    'end':int(fields[4]),
                                    'strand':fields[6]
                                })
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


def parse_graph_range_string(range_str):
    match = re.search(r"^(\S+)@(\S+):(\d+)-(\d+)\(([+-])\)$", range_str)
    if not match:
        return None
    return {
        'chrom': match.group(1),
        'genome': match.group(2),
        'start': int(match.group(3)),
        'end': int(match.group(4)),
        'strand': match.group(5)
    }


def merge_wider_range(range_a, range_b):
    parsed_a = parse_graph_range_string(range_a)
    parsed_b = parse_graph_range_string(range_b)
    if not parsed_a or not parsed_b:
        return range_a

    if (parsed_a['chrom'] == parsed_b['chrom'] and
        parsed_a['genome'] == parsed_b['genome'] and
        parsed_a['strand'] == parsed_b['strand']):

        wider_start = min(parsed_a['start'], parsed_b['start'])
        wider_end = max(parsed_a['end'], parsed_b['end'])
        return f"{parsed_a['chrom']}@{parsed_a['genome']}:{wider_start}-{wider_end}({parsed_a['strand']})"

    len_a = parsed_a['end'] - parsed_a['start']
    len_b = parsed_b['end'] - parsed_b['start']
    return range_a if len_a >= len_b else range_b


# %%
def align_sequence_to_ranges(agc_path, agc_db_path, gmap_path, 
                             sequence, agc_ranges, ranked_genomes=None, 
                             aligner_tool='gmap', minimap_path='minimap2', 
                             min_identity=0.0,
                             verbose=False):
    """
    Calls agc to cut genome fragments corresponding to list of ranges in appropriate agc format.
    Then aligns input sequence to these ranges using either GMAP (pairwise) or Minimap2 (splice).
    If aligner_tool is GMAP and min_identity > 0, ranges below this identity are filtered out.
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

    def keep_non_overlapping_ranges(range_results):
        """Keeps non-overlapping graph ranges; preserves input order."""
        kept = []
        kept_parsed = []

        for range_result in range_results:
            parsed = parse_graph_range_string(range_result)
            if not parsed:
                continue

            overlaps = False
            for kp in kept_parsed:
                if parsed['chrom'] != kp['chrom']:
                    continue
                if parsed['genome'] != kp['genome']:
                    continue

                # overlap if intervals intersect
                if not (parsed['end'] < kp['start'] or parsed['start'] > kp['end']):
                    overlaps = True
                    break

            if not overlaps:
                kept.append(range_result)
                kept_parsed.append(parsed)

        return kept

    def align_single_range_minimap(seqname, query_temp_file_name):
        ref_temp_file_name = write_temp_fasta_file([seqname], {seqname:range_sequences[seqname]})
        range_results = []
        seen_ranges = set()

        command = [minimap_path, "-ax", "splice", "-N", str(MINIMAP2_N_DEFAULT), ref_temp_file_name, query_temp_file_name]
        try:
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, check=True)

            if verbose:
                print(f"# INFO(align_sequence_to_ranges): minimap2 finished for {seqname}")

            for line in result.stdout.splitlines():
                if line.startswith("@"):
                    continue

                cols = line.split("\t")
                if len(cols) < 6:
                    continue

                flag = int(cols[1])
                if flag & 4:
                    continue

                pos = int(cols[3])
                cigar = cols[5]

                aligned_bases = 0
                for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
                    if op in "MIDX=":
                        aligned_bases += int(length)

                minimap_identity = None
                nm_tag = None
                for field in cols[11:]:
                    if field.startswith("NM:i:"):
                        try:
                            nm_tag = int(field.split(":")[-1])
                        except ValueError:
                            nm_tag = None
                        break

                if aligned_bases > 0 and nm_tag is not None:
                    minimap_identity = (1.0 - (nm_tag / aligned_bases)) * 100.0

                if minimap_identity is None and min_identity > 0:
                    if verbose:
                        print(f"# INFO(align_sequence_to_ranges): skip {seqname}, cannot parse minimap identity and min_identity={min_identity}")
                    continue

                if minimap_identity is not None and minimap_identity < min_identity:
                    if verbose:
                        print(f"# INFO(align_sequence_to_ranges): skip {seqname}, minimap identity {minimap_identity:.2f} < min_identity {min_identity}")
                    continue

                ref_len = 0
                for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
                    if op in "MDN=X":
                        ref_len += int(length)

                rangematch = re.search(r"(chr\d[a-zA-Z])_sampleName=([^:]+):(\d+)-(\d+)", seqname)
                if rangematch:
                    contig = rangematch.group(1)
                    genome = rangematch.group(2)
                    chunk_start = int(rangematch.group(3))

                    final_start = chunk_start + pos
                    final_end = final_start + ref_len - 1
                    strand = '-' if (flag & 16) else '+'
                    range_result = f'{contig}@{genome}:{final_start}-{final_end}({strand})'
                    if range_result not in seen_ranges:
                        seen_ranges.add(range_result)
                        range_results.append(range_result)

            raw_hits = len(range_results)
            range_results = keep_non_overlapping_ranges(range_results)
            if verbose:
                print(f"# INFO(align_sequence_to_ranges): {seqname} minimap kept {len(range_results)} non-overlapping hit(s) in this range (raw={raw_hits})")

        except subprocess.CalledProcessError:
            print(f'# ERROR(align_sequence_to_ranges): minimap2 failed')

        finally:
            os.remove(ref_temp_file_name)

        return range_results

    def align_single_range_gmap(seqname):
        range_results = []
        seen_ranges = set()
        range_start = 0
        range_end = 0
        range_strand = '+'

        temp_file_name = write_temp_fasta_file([seqname,'query'],
                                               {seqname:range_sequences[seqname], 'query':sequence})

        command = f"cat {temp_file_name} | {gmap_path} -t 1 -n 1 -2"
        try:
            result = subprocess.run(command,
                                    shell=True, check=False, text=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)

            if result.returncode != 0:
                if verbose:
                    stderr_txt = (result.stderr or '').strip()
                    if stderr_txt:
                        print(f"# WARN(align_sequence_to_ranges): gmap pairwise exited {result.returncode} for {seqname}: {stderr_txt}")
                    else:
                        print(f"# WARN(align_sequence_to_ranges): gmap pairwise exited {result.returncode} for {seqname} (no stderr)")
                return []

            if verbose == True:
                print(f"# INFO(align_sequence_to_ranges): gmap result for {seqname}:\n{result.stdout}")

            paths = re.search(r"Paths \([123456789]+", result.stdout)
            if paths:
                gmap_identity = None
                identity_match = re.search(r"Percent identity:\s*([\d.]+)", result.stdout)
                if not identity_match:
                    identity_match = re.search(r"Identity:\s*([\d.]+)%", result.stdout)
                if identity_match:
                    gmap_identity = float(identity_match.group(1))

                if gmap_identity is None and min_identity > 0:
                    if verbose == True:
                        print(f"# INFO(align_sequence_to_ranges): skip {seqname}, cannot parse GMAP identity and min_identity={min_identity}")
                    return []

                if gmap_identity is not None and gmap_identity < min_identity:
                    if verbose == True:
                        print(f"# INFO(align_sequence_to_ranges): skip {seqname}, identity {gmap_identity:.2f} < min_identity {min_identity}")
                    return []

                rangematch = re.search(r"(chr\d[a-zA-Z])_sampleName=([^:]+):(\d+)-(\d+)", seqname)
                if rangematch:
                    range_prefix = f'{rangematch.group(1)}@{rangematch.group(2)}'
                    range_start = int(rangematch.group(3))
                    range_end = int(rangematch.group(4))
                else:
                    return []

                pos_matches = re.findall(r"Genomic pos: ([\d,]+)..([\d,]+) \(([+-])", result.stdout)
                for pos_match in pos_matches:
                    gmap_start = int(pos_match[0].replace(',',''))
                    gmap_end = int(pos_match[1].replace(',',''))
                    gmap_strand = pos_match[2]
                    range_strand = '+'
                    if gmap_strand == '-':
                        range_strand = '-'
                        gmap_start, gmap_end = gmap_end, gmap_start

                    final_start = range_start + gmap_start
                    final_end = final_start + (gmap_end - gmap_start)
                    range_result = f'{range_prefix}:{final_start}-{final_end}({range_strand})'
                    if range_result not in seen_ranges:
                        seen_ranges.add(range_result)
                        range_results.append(range_result)

                raw_hits = len(range_results)
                range_results = keep_non_overlapping_ranges(range_results)
                if verbose:
                    print(f"# INFO(align_sequence_to_ranges): {seqname} gmap kept {len(range_results)} non-overlapping hit(s) in this range (raw={raw_hits})")
                return range_results

            if verbose == True:
                print(f"# INFO(align_sequence_to_ranges): no valid alignment for {seqname}")

        except Exception as e:
            print(f'# ERROR(align_sequence_to_ranges): gmap pairwise exception for {seqname}: {e}')

        finally:
            os.remove(temp_file_name)

        return []

    # --- Step 2: Align using selected tool ---
    if aligner_tool in ('minimap', 'both'):
        query_temp_file_name = write_temp_fasta_file(['query'], {'query':sequence})
    else:
        query_temp_file_name = ''

    for seqname in range_seqnames:
        minimap_ranges = []
        gmap_ranges = []

        if aligner_tool in ('minimap', 'both'):
            minimap_ranges = align_single_range_minimap(seqname, query_temp_file_name)

        if aligner_tool in ('gmap', 'both'):
            gmap_ranges = align_single_range_gmap(seqname)

        if aligner_tool == 'both':
            # If both tools produced at most one hit, keep old behavior: merge to wider range
            if len(gmap_ranges) <= 1 and len(minimap_ranges) <= 1:
                if gmap_ranges and minimap_ranges:
                    aligned_ranges.append(merge_wider_range(gmap_ranges[0], minimap_ranges[0]))
                elif gmap_ranges:
                    aligned_ranges.extend(gmap_ranges)
                elif minimap_ranges:
                    aligned_ranges.extend(minimap_ranges)
            # If any tool produced multiple hits, keep only tool with higher number of hits
            elif len(gmap_ranges) >= len(minimap_ranges):
                if verbose:
                    print(f"# INFO(align_sequence_to_ranges): {seqname} multi-hit decision in BOTH mode -> keep GMAP ({len(gmap_ranges)} hits) over minimap ({len(minimap_ranges)} hits)")
                aligned_ranges.extend(gmap_ranges)
            else:
                if verbose:
                    print(f"# INFO(align_sequence_to_ranges): {seqname} multi-hit decision in BOTH mode -> keep minimap ({len(minimap_ranges)} hits) over GMAP ({len(gmap_ranges)} hits)")
                aligned_ranges.extend(minimap_ranges)
        elif aligner_tool == 'gmap':
            if gmap_ranges:
                aligned_ranges.extend(gmap_ranges)
        else:
            if minimap_ranges:
                aligned_ranges.extend(minimap_ranges)

    if query_temp_file_name:
        os.remove(query_temp_file_name)

    # Sort aligned ranges using the global get_rank function
    if ranked_genomes:
        aligned_ranges.sort(key=lambda x: get_rank(x, ranked_genomes))

    return ";".join(aligned_ranges)


# %%
def filter_nested_agc_ranges(agc_ranges):
    """Filter out AGC ranges that are nested inside other ranges for the same genome.
    Keeps only non-nested (containing) ranges to avoid redundant overlapping alignments.
    Args: list of agc_range strings like 'chr4H@Du_Li_Huang:605310496-605441755'
    Returns: filtered list with only non-nested ranges
    """
    if len(agc_ranges) <= 1:
        return agc_ranges
    
    # Parse AGC ranges by genome
    ranges_by_genome = {}
    for agc_range in agc_ranges:
        try:
            # Format: chr_locus@genome:start-end
            parts = agc_range.split('@')
            if len(parts) != 2:
                continue
            genome = parts[1].split(':')[0]
            coords = parts[1].split(':')[1]
            start, end = map(int, coords.split('-'))
            
            if genome not in ranges_by_genome:
                ranges_by_genome[genome] = []
            ranges_by_genome[genome].append((start, end, agc_range))
        except (ValueError, IndexError):
            # If parsing fails, keep the range as-is
            pass
    
    # Filter nested ranges within each genome
    filtered = []
    for genome, ranges in ranges_by_genome.items():
        non_nested = []
        for i, (start_i, end_i, range_i) in enumerate(ranges):
            is_nested = False
            for j, (start_j, end_j, _) in enumerate(ranges):
                if i != j:
                    # Check if range i is nested inside range j
                    if start_i >= start_j and end_i <= end_j:
                        is_nested = True
                        break
            if not is_nested:
                non_nested.append(range_i)
        filtered.extend(non_nested)
    
    return filtered if filtered else agc_ranges


# %%
def get_overlap_ranges_reference(gmap_match,hapIDranges,genomes,bed_folder_path,
                                coverage=0.75,mult_mappings='No',all_graph_matches=False,
                                aligner_tool='gmap',
                                bedtools_path='bedtools',grep_path='grep',
                                agc_path='agc', agc_db_path='', 
                                gmap_path='gmap', minimap_path='minimap',
                                ranked_genomes=None,
                                min_identity=0.0,
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
    all_ranges = '.'
    best_match = None
    best_num_ranges = 0
    best_ref_coords = None
    best_range_width = 0

    chrom = gmap_match['chrom']
    genome = gmap_match['genome']
    start = gmap_match['start']
    end = gmap_match['end']
    strand = gmap_match['strand']
    ident = gmap_match['ident']
    cover = gmap_match['cover']

    if verbose == True:
        print(f"# Checking match for {chrom}:{start}-{end} within reference")

    # Clean up path to avoid double slashes
    bed_folder_path = bed_folder_path.rstrip('/')

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
            print(f"# WARN(get_overlap_ranges_reference): several overlaps. Checking ranges (even if you didn't ask for) to pick most accurate reference coords.")

        # Only evaluate when: (1) add_ranges requested OR (2) multiple intersections to choose between
        should_evaluate = aligner_tool or len(intersections) > 1
        
        # If only 1 intersection and no add_ranges, just use it without evaluation
        if not should_evaluate and len(intersections) == 1:
            feature = str(intersections[0]).split("\t")
            feature[-1] = feature[-1].strip()
            ref_chr = feature[0]
            ref_start = int(feature[1])
            ref_end = int(feature[2])
            best_ref_coords = (ref_chr, ref_start, ref_end)
            best_match = ""
        else:
            # Evaluate each bedtools result to pick the best one
            for intersection_idx, feature_str in enumerate(intersections):
                feature = str(feature_str).split("\t")
                feature[-1] = feature[-1].strip()
                ref_chr = feature[0]
                ref_start = int(feature[1])
                ref_end = int(feature[2])
                range_width = ref_end - ref_start
                
                # Extract keys from this intersection
                keys_for_this = {}
                for c in range(0, len(genomes)):
                    k = feature[c+3]
                    if k == ".":
                        continue
                    else:
                        clean_k = k[1:-1]  # remove <>
                        if clean_k not in keys_for_this:
                            keys_for_this[clean_k] = set()
                        keys_for_this[clean_k].add(genomes[c])
                
                # Retrieve ranges for this intersection
                agc_ranges = []
                for k in keys_for_this:
                    for gen in keys_for_this[k]:
                        range_bedfile = f"{bed_folder_path}/{gen}.h.bed"
                        command = f"{grep_path} {k} {range_bedfile}"
                        try:
                            result = subprocess.run(command,
                                    shell=True, check=True, text=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
                            for b in result.stdout.splitlines():
                                bed_data = b.split("\t")
                                if len(bed_data) > 4:
                                    bstart = int(bed_data[1])
                                    bend = int(bed_data[2])
                                    block_length = bend - bstart
                                    padding = int(block_length * 0.05)
                                    padding = max(padding, 1000)
                                    padded_start = max(0, bstart - padding)
                                    padded_end = bend + padding
                                    agc_range = f'{bed_data[0]}@{gen}:{padded_start}-{padded_end}'
                                    agc_ranges.append(agc_range)
                        except subprocess.CalledProcessError:
                            pass
                
                # Evaluate this intersection
                if agc_ranges:
                    agc_ranges = list(set(agc_ranges))
                    agc_ranges = filter_nested_agc_ranges(agc_ranges)
                    # Always use gmap for evaluation, regardless of add_ranges setting
                    eval_aligner_tool = aligner_tool if aligner_tool else 'gmap'
                    aligned_ranges = align_sequence_to_ranges(agc_path, agc_db_path, gmap_path,
                                                                    gmap_match['sequence'], agc_ranges,
                                                                    ranked_genomes=ranked_genomes,
                                                                    aligner_tool=eval_aligner_tool,
                                                                    minimap_path=minimap_path,
                                                                    min_identity=min_identity,
                                                                    verbose=verbose)
                    
                    # Count matches in aligned_ranges
                    num_ranges = len([x for x in aligned_ranges.split(';') if x.strip()]) if aligned_ranges.strip() else 0
                    
                    # Pick best: most matches, or wider if tied
                    if num_ranges > best_num_ranges or (num_ranges == best_num_ranges and range_width > best_range_width):
                        best_num_ranges = num_ranges
                        best_match = aligned_ranges
                        best_ref_coords = (ref_chr, ref_start, ref_end)
                        best_range_width = range_width
                        best_keys = keys_for_this
                        if verbose:
                            print(f"# INFO(get_overlap_ranges_reference): intersection {intersection_idx} has {num_ranges} matched ranges (width={range_width}bp)")
        
        # Use the best match for output
        if best_ref_coords:
            ref_chr, ref_start, ref_end = best_ref_coords
            match_tsv = (
                f'{ref_chr}\t{ref_start}\t{ref_end}\t{strand}\t'
                f'{genome}\t{chrom}\t{start}\t{end}\t{strand}\t{ident}\t{cover}\t{mult_mappings}\t')
            # Only return graph ranges if add_ranges was requested
            all_ranges = best_match if (best_match and aligner_tool) else ""
        else:
            match_tsv = f'.\t.\t.\t.\t{genome}\t{chrom}\t{start}\t{end}\t{strand}\t{ident}\t{cover}\t{mult_mappings}\t'
            all_ranges = ""

    except subprocess.CalledProcessError as e:
        print(f'# ERROR(get_overlap_ranges_reference): {e.cmd} failed: {e.stderr}')

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
def has_graph_overlap(match, reference_name, hapIDranges, hvcf_bed,
                      bedtools_path='bedtools', coverage=0.75):
    """Returns True if match overlaps any graph haplotype range."""

    chrom = match['chrom']
    start = match['start']
    end = match['end']
    genome = match['genome']
    match_interval = f'{chrom}\t{start}\t{end}'

    if genome == reference_name:
        bed_source = hapIDranges
    else:
        bed_source = f"{hvcf_bed}/{genome}.h.bed"

    command = f"{bedtools_path} intersect -a {bed_source} -b stdin -nonamecheck -e -F {coverage} -f {coverage} -sorted"
    try:
        result = subprocess.run(command,
                                shell=True, check=True, text=True,
                                input=match_interval,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        return len(result.stdout.splitlines()) > 0
    except subprocess.CalledProcessError:
        return False


# %%
def run_gmap_genomes(pangenome_genomes, gmap_path, gmap_db, fasta_filename, 
                         min_identity, min_coverage,
                         cores=4, prefix='temp', path='/tmp/', 
                         genomic=False,
                         reference_name='', hapIDranges='', hvcf_bed='',
                         bedtools_path='bedtools', min_coverage_range=75.0,
                         verbose=False):
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
        gmap_matches[seqname] = [] 
     
    pending_seqnames = set(seqnames)

    # loop over genomes hierarchically
    for genome in pangenome_genomes:

        # only gmap sequences that still have no in-graph hit
        g_seqnames = [seqname for seqname in seqnames if seqname in pending_seqnames]
        g_sequences = {seqname: sequences[seqname] for seqname in g_seqnames}

        if len(g_seqnames) == 0:
            if verbose == True:
                print(f"# INFO(run_gmap_genomes): all sequences already have at least one in-graph hit\n")
            break

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

        genome_matches = valid_matches(g_gff_filename,genome,sequences,min_identity,min_coverage,verbose=verbose)
        for seqname in genome_matches:
            gmap_matches[seqname].extend(genome_matches[seqname])

            # stop iterating this sequence only when at least one match overlaps graph ranges
            for match in genome_matches[seqname]:
                has_overlap = has_graph_overlap(
                    match,
                    reference_name,
                    hapIDranges,
                    hvcf_bed,
                    bedtools_path=bedtools_path,
                    coverage=min_coverage_range/100)
                if has_overlap:
                    pending_seqnames.discard(seqname)
                    break

        # clean up temp files
        os.remove(g_gff_filename)
        os.remove(g_fasta_filename)
        
    return gmap_matches

# %%
def get_overlap_ranges_pangenome(gmap_match,hapIDranges,genomes,bedfile,bed_folder_path,
                                coverage=0.75,mult_mappings='No',all_graph_matches=False,
                                aligner_tool='gmap',
                                bedtools_path='bedtools',grep_path='grep',
                                agc_path='agc', agc_db_path='', 
                                gmap_path='gmap', minimap_path='minimap',
                                ranked_genomes=None,
                                min_identity=0.0,
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
    graph_keys = []
    aligned_ranges = '.'
    best_match = ""
    best_num_ranges = 0
    best_ref_coords = None
    best_range_width = 0
    
    # Clean up path to avoid double slashes
    bed_folder_path = bed_folder_path.rstrip('/')
    
    chrom = gmap_match['chrom']
    genome = gmap_match['genome']
    start = gmap_match['start']
    end = gmap_match['end']
    strand = gmap_match['strand']
    ident = gmap_match['ident']
    cover = gmap_match['cover']

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
            print(f"# INFO(get_overlap_ranges_pangenome):\n# chrom\tstart\tend\tstrand\tchecksum\tgenome\tref_chr\tref_start\tref_end\tref_checksum\n{result.stdout}")


        if(len(intersections) == 0):
            match_tsv = f'.\t.\t.\t.\t{genome}\t{chrom}\t{start}\t{end}\t{strand}\t{ident}\t{cover}\t{mult_mappings}\t'
            return match_tsv + aligned_ranges

        elif len(intersections) > 1:
            print(f"# WARN(get_overlap_ranges_pangenome): > several overlaps.  Checking ranges (even if you didn't ask for) to pick most accurate reference coords.")

        # Only evaluate when: (1) add_ranges requested OR (2) multiple intersections to choose between
        should_evaluate = aligner_tool or len(intersections) > 1
        
        # If only 1 intersection and no add_ranges, just use it without evaluation
        if not should_evaluate and len(intersections) == 1:
            feature = str(intersections[0]).split("\t")
            feature[-1] = feature[-1].replace("\n", "")
            ref_chr = feature[6]
            ref_start = int(feature[7])
            ref_end = int(feature[8])
            best_ref_coords = (ref_chr, ref_start, ref_end)
            best_match = ""
        else:
            # Evaluate each bedtools result to pick the best one
            for intersection_idx, feature_str in enumerate(intersections):
                feature = str(feature_str).split("\t")
                feature[-1] = feature[-1].replace("\n", "")
                ref_chr = feature[6]
                ref_start = int(feature[7])
                ref_end = int(feature[8])
                range_width = ref_end - ref_start
                
                # Collect keys from this intersection
                graph_keys_for_this = [feature[4]]
                
                # Look up graph ranges from hapIDranges
                keys_for_ranges = {}
                for graph_key in graph_keys_for_this:
                    command = f"{grep_path} {graph_key} {hapIDranges}"
                    try:
                        result = subprocess.run(command,
                                        shell=True, check=True, text=True, input=match_interval,
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        graph_data = result.stdout.splitlines()
                        
                        if len(graph_data) > 0:
                            for graph_line in graph_data:
                                gfeature = graph_line.split("\t")
                                gfeature[-1] = gfeature[-1].strip()
                                for c in range(0, len(genomes)):
                                    k = gfeature[c+3]
                                    if k != ".":
                                        clean_k = k[1:-1]
                                        if clean_k not in keys_for_ranges:
                                            keys_for_ranges[clean_k] = set()
                                        keys_for_ranges[clean_k].add(genomes[c])
                    except subprocess.CalledProcessError:
                        pass
                
                # Retrieve AGC ranges for this intersection
                agc_ranges = []
                for k in keys_for_ranges:
                    for key_genome in keys_for_ranges[k]:
                        range_bedfile = f"{bed_folder_path}/{key_genome}.h.bed"
                        command = f"{grep_path} {k} {range_bedfile}"
                        try:
                            result = subprocess.run(command,
                                            shell=True, check=True, text=True, input=match_interval,
                                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                            for b in result.stdout.splitlines():
                                bed_data = b.split("\t")
                                if len(bed_data) > 4:
                                    bstart = int(bed_data[1])
                                    bend = int(bed_data[2])
                                    block_length = bend - bstart
                                    padding = int(block_length * 0.05)
                                    padded_start = max(0, bstart - padding)
                                    padded_end = bend + padding
                                    padded_range = f'{bed_data[0]}@{key_genome}:{padded_start}-{padded_end}'
                                    agc_ranges.append(padded_range)
                        except subprocess.CalledProcessError:
                            pass
                
                # Evaluate this intersection
                if agc_ranges:
                    agc_ranges = sorted(set(agc_ranges))
                    agc_ranges = filter_nested_agc_ranges(agc_ranges)
                    # Always use gmap for evaluation, regardless of add_ranges setting
                    eval_aligner_tool = aligner_tool if aligner_tool else 'gmap'
                    aligned_ranges_result = align_sequence_to_ranges(agc_path, agc_db_path, gmap_path,
                                                                gmap_match['sequence'], agc_ranges,
                                                                ranked_genomes=ranked_genomes,
                                                                aligner_tool=eval_aligner_tool,
                                                                minimap_path=minimap_path,
                                                                min_identity=min_identity,
                                                                verbose=verbose)
                    
                    # Count matches 
                    num_ranges = len([x for x in aligned_ranges_result.split(';') if x.strip()]) if aligned_ranges_result.strip() else 0
                    
                    # Pick best: most matches, or wider if tied
                    if num_ranges > best_num_ranges or (num_ranges == best_num_ranges and range_width > best_range_width):
                        best_num_ranges = num_ranges
                        best_match = aligned_ranges_result
                        best_ref_coords = (ref_chr, ref_start, ref_end)
                        best_range_width = range_width
                        if verbose:
                            print(f"# INFO(get_overlap_ranges_pangenome): intersection {intersection_idx} has {num_ranges} matched ranges (width={range_width}bp)")
        
        # Use best match for output
        if best_ref_coords:
            ref_chr, ref_start, ref_end = best_ref_coords
            match_tsv = (
                f'{ref_chr}\t{ref_start}\t{ref_end}\t.\t'
                f'{genome}\t{chrom}\t{start}\t{end}\t{strand}\t{ident}\t{cover}\t{mult_mappings}\t')
        else:
            match_tsv = f'.\t.\t.\t.\t{genome}\t{chrom}\t{start}\t{end}\t{strand}\t{ident}\t{cover}\t{mult_mappings}\t'

    except subprocess.CalledProcessError as e:
        print(f'ERROR(get_overlap_ranges_pangenome): {e.cmd} failed: {e.stderr}')
    
    # Only return graph ranges if add_ranges was requested, otherwise return empty
    final_ranges = best_match if aligner_tool else '.'
    
    return match_tsv + final_ranges


# %%
def main():

    parser = argparse.ArgumentParser(
        description="Map sequences within pangenome graph, returns 1-based coordinates.\n",
        epilog=f"Version: {__version__} {__versiondate__}\n"+
            "Citation: see https://github.com/eead-csic-compbio/barleygraph\n",
        formatter_class=argparse.RawDescriptionHelpFormatter
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
        "--minimap_exe",
        default='minimap2',
        help=f"path to minimap executable, default: minimap2"
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
    
    parser.add_argument('--add_ranges', 
        type=str,
        choices=['gmap', 'minimap', 'both'],
        help='Add all pangenome ranges matching input sequences using specified tool (gmap, minimap, or both)'
    )
    
    parser.add_argument('--force_ranges',
        action='store_true',
        help='When no graph overlap is found, search across all genomes using gmap to find ranges'
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
    minimap_exe   = args.minimap_exe
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
    force_ranges  = args.force_ranges
    
    single_genome = args.single_genome

    ######################################################

    gmap_version = check_gmap_version(gmap_exe)

    print(f"# version: {__version__} {__versiondate__}")
    print(f"# GMAP version: {gmap_version}")
    print(f"# config_file: {args.graph_yaml}")
    print(f"# fasta_file: {fasta_file}")
    print(f"# minimum identity %: {min_identity}")
    print(f"# minimum coverage %: {min_coverage}")
    print(f"# minimum coverage range %: {min_coverage_range}")
    print(f"# add_ranges: {do_add_ranges} (Tool: {aligner_tool if do_add_ranges else 'None'})")
    print(f"# force_ranges: {force_ranges}")
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
        reference_name,
        hapIDranges,
        hvcf_bed,
        bedtools_exe,
        min_coverage_range,
        verbose=verbose_out)
    
    # Track which genomes were actually gmapped (have results in gmap_matches)
    # These should NOT be re-searched by force_ranges
    genomes_already_gmapped = set()
    for seqname in gmap_matches:
        for match in gmap_matches[seqname]:
            genomes_already_gmapped.add(match['genome'])
    


    # print header
    print(f'#query\tref_chr\t\tref_start\tref_end\t\tref_strand\tgenome\t'
          'chr\tstart\tend\tstrand\tperc_ident\tperc_cover\tmultmaps\tgraph_ranges')
    
    # Store main output line for each sequence
    seq_output_lines = {}
    
    # compute graph coordinates for matched sequences
    for seqname in gmap_matches:

        genome_match_counts = {}
        for match in gmap_matches[seqname]:
            gname = match['genome']
            if gname not in genome_match_counts:
                genome_match_counts[gname] = 0
            genome_match_counts[gname] += 1

        in_graph_rows = []
        out_graph_rows = []
        in_graph_genomes = set()
        out_graph_genomes = set()

        for match in gmap_matches[seqname]:

            # multmaps refers only to the genome shown in column 6
            mult_maps = 'Yes' if genome_match_counts.get(match['genome'], 0) > 1 else 'No'

            if(match['genome'] == reference_name):

                matched_coords = get_overlap_ranges_reference(
                    match,
                    hapIDranges,
                    graph_pangenome_genomes,
                    hvcf_bed,
                    coverage=min_coverage_range/100,
                    mult_mappings=mult_maps,
                    all_graph_matches=do_add_ranges,
                    aligner_tool=aligner_tool,
                    bedtools_path=bedtools_exe,
                    grep_path=grep_exe,
                    agc_path=agc_exe,
                    agc_db_path=agc_db,
                    gmap_path=gmap_exe,
                    minimap_path=minimap_exe,
                    ranked_genomes=ranked_pangenome_genomes,
                    min_identity=min_identity,
                    verbose=verbose_out)

            else:
                matched_coords = get_overlap_ranges_pangenome(
                    match,
                    hapIDranges,
                    graph_pangenome_genomes,
                    f"{hvcf_bed}/{match['genome']}.h.bed",
                    hvcf_bed,
                    coverage=min_coverage_range/100,
                    mult_mappings=mult_maps,
                    all_graph_matches=do_add_ranges,
                    aligner_tool=aligner_tool,
                    bedtools_path=bedtools_exe,
                    grep_path=grep_exe,
                    agc_path=agc_exe,
                    agc_db_path=agc_db,
                    gmap_path=gmap_exe,
                    minimap_path=minimap_exe,
                    ranked_genomes=ranked_pangenome_genomes,
                    min_identity=min_identity,
                    verbose=verbose_out)

            outTSV = matched_coords.split("\t")
            matched_tail = "\t".join(outTSV[1:])

            if outTSV[0] in chr_syns:
                ref_chr = chr_syns[outTSV[0]]
            else:
                ref_chr = outTSV[0]

            row = f"{seqname}\t{ref_chr}\t{matched_tail}"

            # outTSV[0] is ref_chr; dot means no overlap with graph ranges
            if outTSV[0] == '.':
                if verbose_out:
                    print(f"# INFO(process_sequences_serial): match has no graph coordinates for {seqname} in genome {match['genome']} ({match['chrom']}:{match['start']}-{match['end']} {match['strand']})")
                
                # For matches without graph overlap, output with simple chr@genome:start-end(strand) range
                simple_range = f"{match['chrom']}@{match['genome']}:{match['start']}-{match['end']}({match['strand']})"
                
                # Reconstruct row with simple range
                parts = row.split("\t")
                parts[-1] = simple_range
                row = "\t".join(parts)
                
                out_graph_rows.append(row)
                out_graph_genomes.add(match['genome'])
            else:
                in_graph_rows.append(row)
                in_graph_genomes.add(match['genome'])

        # Store main output line (prefer in_graph over out_graph)
        if len(in_graph_rows) > 0:
            seq_output_lines[seqname] = in_graph_rows[0]
        elif len(out_graph_rows) > 0:
            seq_output_lines[seqname] = out_graph_rows[0]

    # Post-processing: if force_ranges enabled, collect additional ranges with * marker
    force_ranges_results = {}  # {seqname: [range1*, range2*, ...]}
    
    if force_ranges and do_add_ranges and single_genome == '':
        # Extract all genomes that appear in graph_ranges columns of already-output lines
        genomes_in_graph_ranges = set()
        
        for seqname, output_line in seq_output_lines.items():
            parts = output_line.split("\t")
            if len(parts) > 13:  # Last column is graph_ranges
                ranges_str = parts[-1]
                if ranges_str and ranges_str != '.':
                    for range_item in ranges_str.split(';'):
                        if '@' in range_item:
                            # Extract genome from "chr4H@genome:start-end(strand)"
                            genome = range_item.split('@')[1].split(':')[0]
                            genomes_in_graph_ranges.add(genome)
        
        for seqname in gmap_matches:
            if not gmap_matches[seqname]:
                continue
            
            sequence = gmap_matches[seqname][0].get('sequence', '')
            if not sequence:
                continue
            
            force_ranges_results[seqname] = []
            
            # Case A: Handle genomes already searched by hierarchy but with no graph coordinates
            # These have matches in gmap_matches but don't appear in graph_ranges output
            already_searched_no_ranges = genomes_already_gmapped - genomes_in_graph_ranges
            
            for genome_to_check in sorted(already_searched_no_ranges):
                # Find if this genome has a match in gmap_matches for this sequence
                matching_result = None
                for match in gmap_matches[seqname]:
                    if match['genome'] == genome_to_check:
                        matching_result = match
                        break
                
                if matching_result:
                    # Try to find graph coordinates for this existing match (don't re-run gmap)
                    # This would require calling get_overlap_ranges_* functions
                    # For now, just output the basic match range with * marker (deferred graph lookup)
                    range_with_star = f"{matching_result['chrom']}@{genome_to_check}:{matching_result['start']}-{matching_result['end']}({matching_result['strand']})*"
                    force_ranges_results[seqname].append(range_with_star)
                    
                    if verbose_out:
                        print(f"# INFO(force_ranges): reusing match from {genome_to_check} for {seqname} (already searched, no graph coords)")
            
            # Case B: Search genomes never searched by hierarchy (not in genomes_already_gmapped)
            never_searched = set(graph_pangenome_genomes) - genomes_already_gmapped - genomes_in_graph_ranges
            
            if never_searched:
                if verbose_out:
                    print(f"# INFO(force_ranges): searching new genomes for {seqname}: {', '.join(sorted(never_searched))}")
                
                # For each never-searched genome, run gmap and add simple chr:start-end(strand)* range
                for missing_genome in sorted(never_searched):
                    try:
                        # Write temp FASTA
                        temp_fasta = os.path.join(temp_path, f"{temp_prefix}_{missing_genome}_{seqname}.fa")
                        with open(temp_fasta, 'w') as out_f:
                            out_f.write(f">{seqname}\n{sequence}\n")
                        
                        # Run gmap
                        gmap_cmd = [gmap_exe, '-d', missing_genome, '-D', gmap_db, '-f', 'gff3_gene', temp_fasta,
                                   '-n', '1', '--min-intronlength', '20', '-t', str(ncores)]
                        if genomic:
                            gmap_cmd.append('--nosplicing')
                        
                        result = subprocess.run(gmap_cmd, capture_output=True, text=True, timeout=300)
                        
                        if result.returncode == 0 and result.stdout:
                            # Parse GFF3 to extract match info
                            temp_gff = os.path.join(temp_path, f"{temp_prefix}_{missing_genome}_{seqname}.gff3")
                            with open(temp_gff, 'w') as gff_f:
                                gff_f.write(result.stdout)
                            
                            matches = valid_matches(temp_gff, missing_genome, {seqname: sequence}, min_identity, min_coverage, verbose=False)
                            
                            if matches and seqname in matches:
                                # Take first match, add range with * marker
                                m = matches[seqname][0]
                                range_with_star = f"{m['chrom']}@{missing_genome}:{m['start']}-{m['end']}({m['strand']})*"
                                force_ranges_results[seqname].append(range_with_star)
                                
                                if verbose_out:
                                    print(f"# INFO(force_ranges): found match in {missing_genome} for {seqname}")
                            
                            os.remove(temp_gff)
                    except subprocess.TimeoutExpired:
                        if verbose_out:
                            print(f"# WARN(force_ranges): gmap timeout for {missing_genome}")
                    except Exception as e:
                        if verbose_out:
                            print(f"# WARN(force_ranges): error for {missing_genome}: {e}")
                    finally:
                        try:
                            os.remove(temp_fasta)
                        except:
                            pass
    
    # Print output lines with force_ranges appended to graph_ranges column
    for seqname in seq_output_lines:
        output_line = seq_output_lines[seqname]
        
        # Append force_ranges results to graph_ranges column if any
        if seqname in force_ranges_results and force_ranges_results[seqname]:
            parts = output_line.split("\t")
            # Last column is graph_ranges
            if parts[-1] == '.':
                parts[-1] = ";".join(force_ranges_results[seqname])
            else:
                parts[-1] = parts[-1] + ";" + ";".join(force_ranges_results[seqname])
            output_line = "\t".join(parts)
        
        print(output_line)


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
    from version import __version__ , __versiondate__

    start_time = time.time()
    main()
