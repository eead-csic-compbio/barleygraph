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

    # Pre-compile regex patterns - simplified approach
    # Extract individual fields using simple patterns
    pat_samplename = re.compile(r'SampleName=([^,>]+)')
    pat_regions = re.compile(r'Regions=([^,>]+)')
    pat_checksum = re.compile(r'Checksum=([^,>]+)')
    pat_refrange = re.compile(r'RefRange=([^,>]+)')
    pat_refchecksum = re.compile(r'RefChecksum=([^,>]+)')
    pat_regions_quoted = re.compile(r'Regions="([^"]+)"')

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

        # Extract all fields
        samplename_match = pat_samplename.search(line)
        checksum_match = pat_checksum.search(line)
        refrange_match = pat_refrange.search(line)
        refchecksum_match = pat_refchecksum.search(line)
        
        if not (samplename_match and checksum_match and refrange_match and refchecksum_match):
            if verbose:
                print(f"# WARNING (line {lineno}): Missing required fields in ALT line.")
                print(f"# Line content (first 500 chars): {line[:500]}")
            continue
        
        sample_name = samplename_match.group(1)
        checksum = checksum_match.group(1)
        refchecksum = refchecksum_match.group(1)
        refrange_str = refrange_match.group(1)
        
        # Parse RefRange (chr:start-end)
        refrange_parts = refrange_str.split(':')
        if len(refrange_parts) != 2:
            if verbose:
                print(f"# WARNING (line {lineno}): Invalid RefRange format: {refrange_str}")
            continue
        
        ref_chr = refrange_parts[0]
        ref_coords = refrange_parts[1].split('-')
        if len(ref_coords) != 2:
            if verbose:
                print(f"# WARNING (line {lineno}): Invalid RefRange coordinates: {refrange_parts[1]}")
            continue
        
        try:
            ref_start = int(ref_coords[0])
            ref_end = int(ref_coords[1])
        except ValueError:
            if verbose:
                print(f"# WARNING (line {lineno}): Non-numeric RefRange coordinates: {refrange_parts[1]}")
            continue
        
        # Check for quoted regions (multi-region)
        regions_quoted_match = pat_regions_quoted.search(line)
        if regions_quoted_match:
            # Multi-region format
            regions_str = regions_quoted_match.group(1)
            segments = regions_str.split(',')
            for segment in segments:
                segment = segment.strip()
                seg_parts = segment.split(':')
                if len(seg_parts) == 2:
                    s_chr = seg_parts[0]
                    s_coords = seg_parts[1].split('-')
                    if len(s_coords) == 2:
                        try:
                            s_start = int(s_coords[0])
                            s_end = int(s_coords[1])
                            process_entry(s_chr, s_start, s_end, checksum, sample_name, ref_chr, ref_start, ref_end, refchecksum, is_multi=True)
                        except ValueError:
                            pass
            prev_checksum = checksum
        else:
            # Single-region format
            regions_match = pat_regions.search(line)
            if regions_match:
                regions_str = regions_match.group(1)
                region_parts = regions_str.split(':')
                if len(region_parts) == 2:
                    chr_ = region_parts[0]
                    coords = region_parts[1].split('-')
                    if len(coords) == 2:
                        try:
                            start = int(coords[0])
                            end = int(coords[1])
                            process_entry(chr_, start, end, checksum, sample_name, ref_chr, ref_start, ref_end, refchecksum)
                        except ValueError:
                            if verbose:
                                print(f"# WARNING (line {lineno}): Non-numeric coordinates in Regions: {regions_str}")
                    else:
                        if verbose:
                            print(f"# WARNING (line {lineno}): Invalid Regions format: {regions_str}")
                else:
                    if verbose:
                        print(f"# WARNING (line {lineno}): Invalid Regions format: {regions_str}")
            else:
                if verbose:
                    print(f"# WARNING (line {lineno}): No Regions field found.")
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
