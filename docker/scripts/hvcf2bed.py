import sys
import gzip
import os
import re
import argparse
import glob
from tqdm import tqdm


def natural_chr_key(chrom):
    match = re.match(r'^(?:chr)?(\d+)([A-Za-z]*)$', chrom)
    if match:
        return (0, int(match.group(1)), match.group(2))
    return (1, chrom)


def create_query_ref_plot(output_sorted, genome_name, output_png, verbose=False):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    except Exception as e:
        print(f"WARNING {genome_name}: plotting skipped (matplotlib unavailable: {e})")
        return

    rows = []
    q_chr_max = {}
    r_chr_max = {}

    for line in output_sorted:
        cols = line.split("\t")
        if len(cols) < 10:
            continue

        q_chr = cols[0]
        q_start = int(cols[1])
        q_end = int(cols[2])
        r_chr = cols[6]
        r_start = int(cols[7])
        r_end = int(cols[8])

        rows.append((q_chr, q_start, q_end, r_chr, r_start, r_end))
        q_chr_max[q_chr] = max(q_chr_max.get(q_chr, 0), q_end)
        r_chr_max[r_chr] = max(r_chr_max.get(r_chr, 0), r_end)

    if len(rows) == 0:
        print(f"WARNING {genome_name}: no rows to plot")
        return

    # Build cumulative offsets for chaining chromosomes sequentially
    query_chrs = sorted(q_chr_max.keys(), key=natural_chr_key)
    ref_chrs = sorted(r_chr_max.keys(), key=natural_chr_key)

    q_offsets = {}
    q_offset = 0
    for chr_name in query_chrs:
        q_offsets[chr_name] = q_offset
        q_offset += q_chr_max[chr_name]

    r_offsets = {}
    r_offset = 0
    for chr_name in ref_chrs:
        r_offsets[chr_name] = r_offset
        r_offset += r_chr_max[chr_name]

    # Color map for chromosomes
    all_chrs_sorted = sorted(set(query_chrs + ref_chrs), key=natural_chr_key)
    chr_colors = {}
    colors = plt.cm.tab20([(i / max(len(all_chrs_sorted), 1)) for i in range(len(all_chrs_sorted))])
    for chr_name, color in zip(all_chrs_sorted, colors):
        chr_colors[chr_name] = color

    # Plot with points colored by query chromosome
    plt.figure(figsize=(14, 12), dpi=100)
    
    x_vals, y_vals, colors_list = [], [], []
    for q_chr, q_start, q_end, r_chr, r_start, r_end in rows:
        q_mid = (q_start + q_end) / 2.0
        r_mid = (r_start + r_end) / 2.0
        x_vals.append(q_offsets[q_chr] + q_mid)
        y_vals.append(r_offsets[r_chr] + r_mid)
        colors_list.append(chr_colors[q_chr])
    
    # Batch scatter plot for better performance
    plt.scatter(x_vals, y_vals, s=10, alpha=0.6, c=colors_list)

    # Add chromosome boundary lines on X-axis
    for chr_name in query_chrs[1:]:
        x_boundary = q_offsets[chr_name]
        plt.axvline(x=x_boundary, color='gray', linestyle='--', linewidth=0.5, alpha=0.3)

    # Add chromosome boundary lines on Y-axis
    for chr_name in ref_chrs[1:]:
        y_boundary = r_offsets[chr_name]
        plt.axhline(y=y_boundary, color='gray', linestyle='--', linewidth=0.5, alpha=0.3)

    # Create legend with chromosome colors
    legend_patches = [mpatches.Patch(color=chr_colors[chr_name], label=chr_name) 
                      for chr_name in all_chrs_sorted]
    plt.legend(handles=legend_patches, bbox_to_anchor=(1.05, 1), loc='upper left', 
               fontsize=8, title="Query chromosomes")

    plt.xlabel("Query chromosomes (chained, 0 within each chr)")
    plt.ylabel("Reference chromosomes (chained, 0 within each chr)")
    plt.title(f"{genome_name}: query vs reference haplotype block coordinates")
    plt.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.savefig(output_png, bbox_inches='tight')
    plt.close()

    if verbose:
        print(f"Plot written: {output_png} ({len(rows)} points)")



def analyze_overlap_type(prev_start, prev_end, curr_start, curr_end):
    """Classify overlap as NESTED or PARTIAL only if redundant region > 1 bp.
    Returns overlap type and size, or None if overlap is just touching (≤1 bp).
    """
    # Calculate actual overlap region
    overlap_start = max(prev_start, curr_start)
    overlap_end = min(prev_end, curr_end)
    overlap_size = overlap_end - overlap_start
    
    # Ignore trivial overlaps (just touching or zero-length)
    if overlap_size <= 1:
        return None
    
    if curr_start >= prev_start and curr_end <= prev_end:
        return "NESTED_IN_PREV"  # Current nested in previous
    elif prev_start >= curr_start and prev_end <= curr_end:
        return "NESTED_PREV_IN_CURR"  # Previous nested in current
    else:
        return "PARTIAL"  # Actual end-to-end overlap


def detect_overlaps(output_sorted):
    """Detect overlaps in sorted BED-like lines.
    Returns tuple: (query_overlaps, ref_overlaps)
    Each overlap is a dict with current and previous block details.
    """

    query_overlaps = []
    ref_overlaps = []

    prev_query_by_chr = {}
    prev_ref_by_chr = {}

    for line in output_sorted:
        cols = line.split("\t")
        if len(cols) < 10:
            continue

        q_chr = cols[0]
        q_start = int(cols[1])
        q_end = int(cols[2])
        q_checksum = cols[4]

        r_chr = cols[6]
        r_start = int(cols[7])
        r_end = int(cols[8])
        r_checksum = cols[9]

        # Query-space overlap: current starts before or at previous end
        if q_chr in prev_query_by_chr:
            prev = prev_query_by_chr[q_chr]
            if q_start <= prev['end']:
                overlap_type = analyze_overlap_type(prev['start'], prev['end'], q_start, q_end)
                if overlap_type is not None:  # Only store if redundant (> 1 bp overlap)
                    query_overlaps.append({
                        'chr': q_chr,
                        'prev_start': prev['start'],
                        'prev_end': prev['end'],
                        'prev_checksum': prev['checksum'],
                        'start': q_start,
                        'end': q_end,
                        'checksum': q_checksum,
                        'type': overlap_type,
                    })

        # Ref-space overlap: current starts before or at previous end
        if r_chr in prev_ref_by_chr:
            prev = prev_ref_by_chr[r_chr]
            if r_start <= prev['end']:
                overlap_type = analyze_overlap_type(prev['start'], prev['end'], r_start, r_end)
                if overlap_type is not None:  # Only store if redundant (> 1 bp overlap)
                    ref_overlaps.append({
                        'chr': r_chr,
                        'prev_start': prev['start'],
                        'prev_end': prev['end'],
                        'prev_checksum': prev['checksum'],
                        'start': r_start,
                        'end': r_end,
                        'checksum': r_checksum,
                        'type': overlap_type,
                    })

        prev_query_by_chr[q_chr] = {'start': q_start, 'end': q_end, 'checksum': q_checksum}
        prev_ref_by_chr[r_chr] = {'start': r_start, 'end': r_end, 'checksum': r_checksum}

    return query_overlaps, ref_overlaps

def parse_vcf_lines(lines, genome_name, verbose=False):
    viewed_checksums = set() 
    empty_ranges = 0
    skipped_checksums_bug = 0
    output_lines = []
    prev_checksum = None
    
    # Statistics tracking
    pattern_counts = {'multi_region': 0, 'standard': 0, 'swapped': 0, 'no_match': 0}

    # --- REGEX PATTERNS ---
    
    # 1. Multi-region format (quoted Regions)
    pat_multi = re.compile(
        r'SampleName=([^,>]+).*?Regions="([^"]+)".*?Checksum=([^,]+).*?RefChecksum=([^,>]+).*?RefRange=([^:]+):(\d+)-(\d+)'
    )
    
    # 2. Standard format (RefChecksum BEFORE RefRange)
    pat_std = re.compile(
        r'SampleName=([^,]+).*?Regions=([^:]+):(\d+)-(\d+)(?!.*Regions=).*?Checksum=([^,]+).*?RefChecksum=([^,]+).*?RefRange=([^:]+):(\d+)-(\d+)'
    )

    # 3. SWAPPED format (RefRange BEFORE RefChecksum) - This fixes your Ler0_1 issue
    pat_swapped = re.compile(
        r'SampleName=([^,]+).*?Regions=([^:]+):(\d+)-(\d+)(?!.*Regions=).*?Checksum=([^,]+).*?RefRange=([^:]+):(\d+)-(\d+).*?RefChecksum=([^,>]+)'
    )

    def clean_field(value):
        if value is None:
            return None
        return value.strip().strip('"')
        # Sometimes hvcf is buggy and has extra quotes or spaces, so we clean it up

    def process_entry(chr_, start, end, checksum, genome, ref_chr, ref_start, ref_end, ref_checksum, is_multi=False):
        nonlocal prev_checksum, empty_ranges, skipped_checksums_bug, verbose

        chr_ = clean_field(chr_)
        checksum = clean_field(checksum)
        genome = clean_field(genome)
        ref_chr = clean_field(ref_chr)
        ref_checksum = clean_field(ref_checksum)
        
        # Validate RefRange logic
        try:
            ref_start, ref_end = int(ref_start), int(ref_end)
        except ValueError:
            return

        # For multi-region entries, we rely on RefRange being valid
        if is_multi:
            ref_length = abs(ref_end - ref_start)
            if ref_length <= 0:
                empty_ranges += 1
                return
            
            output_lines.append(
                f"{ref_chr}\t{ref_start}\t{ref_end}\t+\t{checksum}\t{genome}\t{ref_chr}\t{ref_start}\t{ref_end}\t{ref_checksum}"
            )
            return
        
        # For standard/swapped formats
        start, end = int(start), int(end)
        length = end - start
        
        strand = "+"
        if start > end:
            strand = "-"
            start, end = end, start
            length = end - start
        
        if length < 0:
            empty_ranges += 1
            return 

        # Check for duplicates
        if checksum == prev_checksum:
            skipped_checksums_bug += 1
            return 
        elif checksum in viewed_checksums:
            pass # Keep going if seen before but not sequential

        output_lines.append(
            f"{chr_}\t{start}\t{end}\t{strand}\t{checksum}\t{genome}\t{ref_chr}\t{ref_start}\t{ref_end}\t{ref_checksum}"
        )
        
        prev_checksum = checksum
        viewed_checksums.add(checksum)

    # --- MAIN PARSING LOOP ---
    # Convert to list to get total line count for progress bar
    lines_list = list(lines)
    total_lines = len(lines_list)
    
    # Create progress bar
    pbar = tqdm(enumerate(lines_list, start=1), total=total_lines, 
                desc=f"Processing {genome_name}", unit=" lines", 
                disable=not verbose)
    
    for lineno, line in pbar:
        line = line.strip()
        
        if not line.startswith("##ALT"):
            continue

        # Try Pattern 1: Multi-region
        m2 = pat_multi.search(line)
        if m2:
            genome, regions, checksum, ref_checksum, ref_chr, ref_start, ref_end = m2.groups()
            process_entry(None, None, None, checksum, genome, ref_chr, ref_start, ref_end, ref_checksum, is_multi=True)
            prev_checksum = checksum
            pattern_counts['multi_region'] += 1
            continue

        # Try Pattern 2: Standard (RefChecksum first)
        m1 = pat_std.search(line)
        if m1:
            genome, chr_, start, end, checksum, ref_checksum, ref_chr, ref_start, ref_end = m1.groups()
            process_entry(chr_, start, end, checksum, genome, ref_chr, ref_start, ref_end, ref_checksum)
            pattern_counts['standard'] += 1
            continue

        # Try Pattern 3: Swapped (RefRange first)
        m3 = pat_swapped.search(line)
        if m3:
            genome, chr_, start, end, checksum, ref_chr, ref_start, ref_end, ref_checksum = m3.groups()
            process_entry(chr_, start, end, checksum, genome, ref_chr, ref_start, ref_end, ref_checksum)
            pattern_counts['swapped'] += 1
            continue

        pattern_counts['no_match'] += 1
        if verbose:
            pbar.write(f"WARNING (line {lineno}): No regex match for ALT line")
    
    pbar.close()

    if verbose:
        print(f"\n✓ Processed {len(output_lines)} haplotype blocks")
        print(f"  Pattern matches: Multi-region={pattern_counts['multi_region']}, Standard={pattern_counts['standard']}, Swapped={pattern_counts['swapped']}")
        print(f"  Issues: Empty ranges={empty_ranges}, Duplicate checksums={skipped_checksums_bug}, No match={pattern_counts['no_match']}")

    return output_lines

def main(args=None):
    parser = argparse.ArgumentParser(description='Convert HVCF files to BED format')
    parser.add_argument('vcf_folder', help='Path to folder containing h.vcf files')
    parser.add_argument('genome_name', nargs='?', default=None, help='Genome name (optional)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
    parser.add_argument('--no-plot', action='store_true', help='Skip query-vs-reference coordinate plot generation')
    
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    vcf_folder = args.vcf_folder
    genome_name = args.genome_name
    verbose = getattr(args, 'verbose', False)
    do_plot = not getattr(args, 'no_plot', False)

    if not os.path.isdir(vcf_folder):
        print(f"ERROR: Folder not found: {vcf_folder}")
        sys.exit(1)

    if genome_name is None:
        # Find all files
        files_gz = glob.glob(os.path.join(vcf_folder, "*.h.vcf.gz"))
        files_plain = glob.glob(os.path.join(vcf_folder, "*.h.vcf"))
        
        # Unique base names
        targets = set()
        for f in files_gz + files_plain:
            base = os.path.basename(f).replace('.h.vcf.gz', '').replace('.h.vcf', '')
            targets.add(base)
            
        for target in targets:
            process_single_file(vcf_folder, target, verbose, do_plot)
    else:
        process_single_file(vcf_folder, genome_name, verbose, do_plot)

def process_single_file(vcf_folder, genome_name, verbose=False, do_plot=True):
    vcf_file_gz = os.path.join(vcf_folder, f"{genome_name}.h.vcf.gz")
    vcf_file_plain = os.path.join(vcf_folder, f"{genome_name}.h.vcf")

    if os.path.exists(vcf_file_gz):
        vcf_file = vcf_file_gz
        file_type = "gzipped"
    elif os.path.exists(vcf_file_plain):
        vcf_file = vcf_file_plain
        file_type = "plain text"
    else:
        print(f"ERROR: VCF file not found for {genome_name}")
        return

    bed_file = os.path.join(vcf_folder, f"{genome_name}.h.bed")
    plot_file = os.path.join(vcf_folder, f"{genome_name}.h.coordplot.png")

    try:
        file_size = os.path.getsize(vcf_file) / (1024*1024)  # MB

        if verbose:
            print(f"\n{'='*60}")
            print(f"Processing: {genome_name}")
            print(f"File: {os.path.basename(vcf_file)} ({file_type})")
            print(f"Size: {file_size:.1f} MB")
            print(f"{'='*60}")

            if file_type == "plain text":
                alt_count = 0
                total_lines = 0
                with open(vcf_file, "r") as f:
                    for line in f:
                        total_lines += 1
                        if line.startswith("##ALT="):
                            alt_count += 1
                print(f"⚠️  UNCOMPRESSED FILE - Expected slower performance")
                print(f"   Total lines: {total_lines:,}")
                print(f"   ALT headers: {alt_count:,}")
                print(f"   → Large uncompressed files are slow due to disk I/O")
                print(f"   → Many ALT headers require regex parsing on long lines")

        if vcf_file.endswith(".gz"):
            with gzip.open(vcf_file, "rt") as f:
                output = parse_vcf_lines(f, genome_name, verbose=verbose)
        else:
            if verbose:
                print(f"Loading {file_size:.1f} MB into memory...")
            with open(vcf_file, "r") as f:
                output = parse_vcf_lines(f, genome_name, verbose=verbose)

        if verbose:
            print(f"Sorting {len(output)} records...")

        output_sorted = sorted(output, key=lambda x: (x.split("\t")[0], int(x.split("\t")[1])))

        query_overlaps, ref_overlaps = detect_overlaps(output_sorted)
        overlap_report = os.path.join(vcf_folder, f"{genome_name}.h.overlap_report.tsv")

        with open(overlap_report, "w") as rep:
            rep.write("space\tchr\tprev_start\tprev_end\tprev_checksum\tstart\tend\tchecksum\ttype\n")
            for ov in query_overlaps:
                rep.write(
                    f"query\t{ov['chr']}\t{ov['prev_start']}\t{ov['prev_end']}\t{ov['prev_checksum']}\t{ov['start']}\t{ov['end']}\t{ov['checksum']}\t{ov['type']}\n"
                )
            for ov in ref_overlaps:
                rep.write(
                    f"reference\t{ov['chr']}\t{ov['prev_start']}\t{ov['prev_end']}\t{ov['prev_checksum']}\t{ov['start']}\t{ov['end']}\t{ov['checksum']}\t{ov['type']}\n"
                )

        if query_overlaps or ref_overlaps:
            # Count overlap types
            query_nested_in = sum(1 for o in query_overlaps if o['type'] == 'NESTED_IN_PREV')
            query_nested_out = sum(1 for o in query_overlaps if o['type'] == 'NESTED_PREV_IN_CURR')
            query_partial = sum(1 for o in query_overlaps if o['type'] == 'PARTIAL')
            ref_nested_in = sum(1 for o in ref_overlaps if o['type'] == 'NESTED_IN_PREV')
            ref_nested_out = sum(1 for o in ref_overlaps if o['type'] == 'NESTED_PREV_IN_CURR')
            ref_partial = sum(1 for o in ref_overlaps if o['type'] == 'PARTIAL')
            
            print(
                f"WARNING {genome_name}: overlaps detected\n"
                f"  Query (total={len(query_overlaps)}): {query_nested_in} nested-in + {query_nested_out} nested-out + {query_partial} partial\n"
                f"  Reference (total={len(ref_overlaps)}): {ref_nested_in} nested-in + {ref_nested_out} nested-out + {ref_partial} partial\n"
                f"  Report: {overlap_report}"
            )
            if verbose:
                for ov in query_overlaps[:5]:
                    print(
                        f"  QUERY [{ov['type']}]: "
                        f"{ov['chr']}:{ov['prev_start']}-{ov['prev_end']} ({ov['prev_checksum']}) "
                        f"vs {ov['chr']}:{ov['start']}-{ov['end']} ({ov['checksum']})"
                    )
                for ov in ref_overlaps[:5]:
                    print(
                        f"  REF [{ov['type']}]: "
                        f"{ov['chr']}:{ov['prev_start']}-{ov['prev_end']} ({ov['prev_checksum']}) "
                        f"vs {ov['chr']}:{ov['start']}-{ov['end']} ({ov['checksum']})"
                    )

        if verbose:
            print(f"Writing to: {os.path.basename(bed_file)}")

        with open(bed_file, "w") as out:
            out.write("#chrom\tstart\tend\tstrand\tchecksum\tgenome\tref_chr\tref_start\tref_end\tref_checksum\n")
            out.write("\n".join(output_sorted) + "\n")

        if do_plot:
            create_query_ref_plot(output_sorted, genome_name, plot_file, verbose=verbose)

        print(f"DONE {genome_name}: {len(output_sorted)} haplotype blocks")

    except Exception as e:
        print(f"ERROR processing {genome_name}: {e}")

if __name__ == "__main__":
    main()