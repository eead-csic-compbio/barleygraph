import sys
import gzip
import os
import re
import argparse
import glob

def parse_vcf_lines(lines, genome_name, verbose=False):
    viewed_checksums = set() 
    empty_ranges = 0
    skipped_checksums_bug = 0
    output_lines = []
    prev_checksum = None

    # Pre-compile regex patterns (order matters!)
    # Multi-region format (quoted Regions field) - more flexible pattern
    pat_multi = re.compile(
        r'SampleName=([^,>]+).*?Regions="([^"]+)".*?Checksum=([^,]+).*?RefChecksum=([^,>]+).*?RefRange=([^:]+):(\d+)-(\d+)'
    )
    
    # Standard format (RefChecksum before RefRange)
    pat_std = re.compile(
        r'SampleName=([^,]+).*?Regions=([^:]+):(\d+)-(\d+)(?!.*Regions=).*?Checksum=([^,]+).*?RefChecksum=([^,]+).*?RefRange=([^:]+):(\d+)-(\d+)'
    )

    def process_entry(chr_, start, end, checksum, genome, ref_chr, ref_start, ref_end, ref_checksum, is_multi=False):
        nonlocal prev_checksum, empty_ranges, skipped_checksums_bug, verbose
        
        # For multi-region entries, we get the RefRange which is already valid
        if is_multi:
            ref_start, ref_end = int(ref_start), int(ref_end)
            ref_length = abs(ref_end - ref_start)
            
            if ref_length <= 0:
                empty_ranges += 1
                if verbose:
                    print(f"# Skipping {checksum}: RefRange length {ref_length} is invalid")
                return
            
            # For multi-region, output single entry using RefRange
            output_lines.append(
                f"{ref_chr}\t{ref_start}\t{ref_end}\t+\t{checksum}\t{genome}\t{ref_chr}\t{ref_start}\t{ref_end}\t{ref_checksum}"
            )
            if verbose:
                print(f"# Added multi-region entry: {checksum}")
            return
        
        # For standard format
        start, end = int(start), int(end)
        length = end - start
        
        strand = "+"
        if start > end:
            strand = "-"
            start, end = end, start
            length = end - start
        
        # Only filter if length is strictly negative
        if length < 0:
            empty_ranges += 1
            return 

        # Check for duplicates
        if checksum == prev_checksum:
            if verbose:
                print(f"# WARNING: duplicate checksum {checksum} skipped")
            skipped_checksums_bug += 1
            return 
        elif checksum in viewed_checksums:
            if verbose:
                print(f"# WARNING: checksum {checksum} already seen")

        output_lines.append(
            f"{chr_}\t{start}\t{end}\t{strand}\t{checksum}\t{genome}\t{ref_chr}\t{ref_start}\t{ref_end}\t{ref_checksum}"
        )
        
        prev_checksum = checksum
        viewed_checksums.add(checksum)

    # Main Loop
    for lineno, line in enumerate(lines, start=1):
        line = line.strip()
        
        if not line.startswith("##ALT"):
            continue

        # Check multi-region FIRST (has quoted regions with commas inside quotes)
        m2 = pat_multi.search(line)
        if m2:
            genome, regions, checksum, ref_checksum, ref_chr, ref_start, ref_end = m2.groups()
            
            # For multi-region, collect all segment lengths to see if any are valid
            segments = regions.split(",")
            valid_segment_lengths = []
            
            for segment in segments:
                seg_match = re.match(r'([^:]+):(\d+)-(\d+)', segment.strip())
                if seg_match:
                    s_chr, s_start, s_end = seg_match.groups()
                    seg_length = abs(int(s_end) - int(s_start))
                    valid_segment_lengths.append(seg_length)
            
            if verbose:
                non_zero = [l for l in valid_segment_lengths if l > 0]
                zero_count = len([l for l in valid_segment_lengths if l == 0])
                if zero_count > 0:
                    print(f"# Matched multi-region format for {checksum} ({len(segments)} segments: {zero_count} zero-length, {len(non_zero)} valid)")
                else:
                    print(f"# Matched multi-region format for {checksum} ({len(segments)} segments, all valid)")
            
            # Process using RefRange (which should be valid)
            process_entry(None, None, None, checksum, genome, ref_chr, ref_start, ref_end, ref_checksum, is_multi=True)
            prev_checksum = checksum
            continue

        # Check standard format (single region, no quotes)
        m1 = pat_std.search(line)
        if m1:
            genome, chr_, start, end, checksum, ref_checksum, ref_chr, ref_start, ref_end = m1.groups()
            process_entry(chr_, start, end, checksum, genome, ref_chr, ref_start, ref_end, ref_checksum)
            continue

        if verbose:
            print(f"# WARNING (line {lineno}): No match for ALT line.")
            print(f"# Line content (first 500 chars): {line[:500]}")

    if verbose:
        print(f"# Skipped {skipped_checksums_bug} duplicate checksums")
        print(f"# Skipped {empty_ranges} empty ranges")
        print(f"# Processed {len(output_lines)} haplotype blocks")

    return output_lines

def main(args=None):
    """
    Convert HVCF files to BED format.
    If no genome_name is provided, converts all .h.vcf and .h.vcf.gz files in the folder.
    """
    if args is None:
        parser = argparse.ArgumentParser(
            description='Convert HVCF files to BED format',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='''
Examples:
  # Convert a specific genome
  hvcf2bed /path/to/vcf/folder GenomeName
  
  # Convert all h.vcf files in folder
  hvcf2bed /path/to/vcf/folder
  
  # With verbose output
  hvcf2bed /path/to/vcf/folder -v
            '''
        )
        parser.add_argument('vcf_folder', help='Path to folder containing h.vcf files')
        parser.add_argument('genome_name', nargs='?', default=None, help='Genome name (optional). If not provided, converts all h.vcf files in folder')
        parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
        args = parser.parse_args()

    vcf_folder = args.vcf_folder
    genome_name = args.genome_name
    verbose = getattr(args, 'verbose', False)

    if not os.path.isdir(vcf_folder):
        print(f"❌ ERROR: Folder not found: {vcf_folder}")
        sys.exit(1)

    # If genome_name not provided, find all .h.vcf and .h.vcf.gz files
    if genome_name is None:
        vcf_files = []
        
        # Find all .h.vcf.gz files
        vcf_files.extend(glob.glob(os.path.join(vcf_folder, "*.h.vcf.gz")))
        
        # Find all .h.vcf files (but not .h.vcf.gz)
        vcf_files.extend(glob.glob(os.path.join(vcf_folder, "*.h.vcf")))
        
        # Remove duplicates and .h.vcf files if .h.vcf.gz version exists
        vcf_files_dict = {}
        for vcf_file in vcf_files:
            base_name = os.path.basename(vcf_file).replace('.h.vcf.gz', '').replace('.h.vcf', '')
            if base_name not in vcf_files_dict:
                vcf_files_dict[base_name] = vcf_file
            else:
                # Prefer .h.vcf.gz over .h.vcf
                if vcf_file.endswith('.h.vcf.gz'):
                    vcf_files_dict[base_name] = vcf_file
        
        vcf_files = list(vcf_files_dict.values())
        
        if not vcf_files:
            print(f"❌ ERROR: No .h.vcf or .h.vcf.gz files found in {vcf_folder}")
            sys.exit(1)
        
        if verbose:
            print(f"# Found {len(vcf_files)} h.vcf file(s) to process:")
            for vcf_file in vcf_files:
                print(f"#   {os.path.basename(vcf_file)}")
        
        # Process each file
        for vcf_file in vcf_files:
            base_name = os.path.basename(vcf_file).replace('.h.vcf.gz', '').replace('.h.vcf', '')
            process_single_file(vcf_folder, base_name, verbose)
    else:
        # Process single genome
        process_single_file(vcf_folder, genome_name, verbose)


def process_single_file(vcf_folder, genome_name, verbose=False):
    """Helper function to process a single HVCF file"""
    vcf_file_gz = os.path.join(vcf_folder, f"{genome_name}.h.vcf.gz")
    vcf_file_plain = os.path.join(vcf_folder, f"{genome_name}.h.vcf")
    
    if os.path.exists(vcf_file_gz):
        vcf_file = vcf_file_gz
    elif os.path.exists(vcf_file_plain):
        vcf_file = vcf_file_plain
    else:
        if verbose:
            print(f"# ERROR: VCF file not found for {genome_name}")
            print(f"# Looked for: {vcf_file_gz} or {vcf_file_plain}")
        else:
            print(f"❌ ERROR: VCF file not found for {genome_name}")
            print(f"  Looked for: {vcf_file_gz} or {vcf_file_plain}")
        return

    bed_file = os.path.join(vcf_folder, f"{genome_name}.h.bed")

    if verbose:
        print(f"# Processing: {os.path.basename(vcf_file)}")

    try:
        if vcf_file.endswith(".gz"):
            with gzip.open(vcf_file, "rt") as f:
                output = parse_vcf_lines(f, genome_name, verbose=verbose)
        else:
            with open(vcf_file, "r") as f:
                output = parse_vcf_lines(f, genome_name, verbose=verbose)

        # Sort output lines by chromosome and start position
        output_sorted = sorted(
            output,
            key=lambda x: (x.split("\t")[0], int(x.split("\t")[1]))
        )

        with open(bed_file, "w") as out:
            # Write header
            out.write("#chrom\tstart\tend\tstrand\tchecksum\tgenome\tref_chr\tref_start\tref_end\tref_checksum\n")
            out.write("\n".join(output_sorted) + "\n")

        if verbose:
            print(f"# Converted to: {os.path.basename(bed_file)} ({len(output_sorted)} entries)")
        else:
            print(f"✓ {genome_name}: {len(output_sorted)} haplotype blocks")

    except Exception as e:
        if verbose:
            print(f"# ERROR processing {genome_name}: {e}")
            import traceback
            traceback.print_exc()
        else:
            print(f"❌ ERROR processing {genome_name}: {e}")
        return


if __name__ == "__main__":
    main()
