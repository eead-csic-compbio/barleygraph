# Date: 26.03.2025
# Works with PHGv2 
# citation: 10.1093/bioinformatics/btac410
# Used in version phg version 2.4.52.207

# Last update: 28/05/2025
# Time for a full run: Approx 3 days
# Space used was 500 gb with more new gmap versions. Now should be less

# Update: 28/07/2025
# Updated hvcf2bed.py script to addapt to barleymap structure
# Use gmap version installed into esplandian, using the version of gmap 2013-08-31 (less storage, and same that barleymap version)

# Update 04/08/2025
# hvcf2bed.py script now sorts the output by chromosome and start position

# Update 23/09/2025
# hvcf2bed.py script now warns if a checksum is duplicated in the same vcf file, and delete very short ranges (length <=1)
# Commented out build of old kmer index, because is being deprecated in the new versions of phg

###
# _________________________________________________________________________________________________
###

# Remember to run with proper conda environment
echo "Remember to set up properly the conda env: phgv2-conda"
echo 
conda activate phgv2-conda
phg --version

# Ensure proper gmap version to be used, out of conda env
gmap_build_command=/usr/local/bin/gmap_build

# Run gmap_build with no arguments to get version info
gmap_version_output=$($gmap_build_command 2>&1)
if ! echo "$gmap_version_output" | grep -q "GMAP package, version"; then
    echo "Error: Could not determine gmap_build version."
    exit 1
fi
if ! echo "$gmap_version_output" | grep -q "2013-08-31"; then
    echo "Error: gmap_build version is not 2013-08-31. Found:"
    echo "$gmap_version_output" | grep "GMAP package, version"
    exit 1
fi

gmap_build_version=$(echo "$gmap_version_output" | grep "GMAP package, version")
echo "Using gmap_build version: $gmap_build_version"

# Update the proper conda env:

# phg setup-environment
# Commented to avoid it during tests

# Generate necessary folders
mkdir -p Pan20/data
mkdir -p Pan20/output
mkdir -p Pan20/vcf_dbs
mkdir -p Pan20/output/alignment_files
mkdir -p Pan20/output/vcf_files
mkdir -p Pan20/output/read_mapping
mkdir -p Pan20/output/imputed_vcf_files
mkdir -p gmap_db

# Initialize the database
phg initdb --db-path Pan20/vcf_dbs/

# Get the genomes
# from the folder:
# https://galaxy-web.ipk-gatersleben.de/libraries/folders/F3f5830403180d620/page/1
# citation: https://www.nature.com/articles/s41586-020-2947-8

# Download the genomes
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=b8a0d6158b9961df&folder_ids=    # Akashinriki
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=0d16186aaff7cbfd&folder_ids=    # B1K-04-12
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=6d4c8d216c52b289&folder_ids=    # Barke
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=68013dab1c13fb37&folder_ids=    # Chiba
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=6d9affd96770ffb9&folder_ids=    # Du_Li_Huang
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=0998682e795f41c2&folder_ids=    # GoldenPromise
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=f30a35c999095ed7&folder_ids=    # Hockett
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=9b5c597dcbb59371&folder_ids=    # HOR_10350
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=c08f58209a47fe2d&folder_ids=    # HOR_13821         ##### Take a look, chr name no ENA
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=59606d2a36c77a56&folder_ids=    # HOR_13942         ##### Take a look, chr name no ENA
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=3f00d864a8cefbf4&folder_ids=    # HOR_21599
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=be8c8ac3dbd1dd54&folder_ids=    # HOR_3081 
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=b2486a20bc56b90f&folder_ids=    # HOR_3365
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=83e633fc7a5c928a&folder_ids=    # HOR_7552
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=795b4cb2a2107ba7&folder_ids=    # HOR_8148
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=c86c1b73aa7102dd&folder_ids=    # HOR_9043
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=1ac3b66e11b8fe50&folder_ids=    # Igri
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=e2132aef71b11dbf&folder_ids=    # OUN333
# https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=0d72ca01c763d02d&folder_ids=    # Planet

# https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_904849725.1?download=true&gzip=true   MorexV3

# Akashinriki
if [ ! -f Pan20/data/Akashinriki.fasta.gz ] && [ ! -f Pan20/data/Akashinriki.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=b8a0d6158b9961df&folder_ids=" -O Pan20/data/Akashinriki.fasta.gz
fi

# B1K-04-12
if [ ! -f Pan20/data/B1K-04-12.fasta.gz ] && [ ! -f Pan20/data/B1K-04-12.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=0d16186aaff7cbfd&folder_ids=" -O Pan20/data/B1K-04-12.fasta.gz
fi

# Barke
if [ ! -f Pan20/data/Barke.fasta.gz ] && [ ! -f Pan20/data/Barke.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=6d4c8d216c52b289&folder_ids=" -O Pan20/data/Barke.fasta.gz
fi

# Chiba
if [ ! -f Pan20/data/Chiba.fasta.gz ] && [ ! -f Pan20/data/Chiba.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=68013dab1c13fb37&folder_ids=" -O Pan20/data/Chiba.fasta.gz
fi

# Du_Li_Huang
if [ ! -f Pan20/data/Du_Li_Huang.fasta.gz ] && [ ! -f Pan20/data/Du_Li_Huang.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=6d9affd96770ffb9&folder_ids=" -O Pan20/data/Du_Li_Huang.fasta.gz
fi

# GoldenPromise
if [ ! -f Pan20/data/GoldenPromise.fasta.gz ] && [ ! -f Pan20/data/GoldenPromise.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=0998682e795f41c2&folder_ids=" -O Pan20/data/GoldenPromise.fasta.gz
fi

# Hockett
if [ ! -f Pan20/data/Hockett.fasta.gz ] && [ ! -f Pan20/data/Hockett.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=f30a35c999095ed7&folder_ids=" -O Pan20/data/Hockett.fasta.gz
fi

# HOR_10350
if [ ! -f Pan20/data/HOR_10350.fasta.gz ] && [ ! -f Pan20/data/HOR_10350.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=9b5c597dcbb59371&folder_ids=" -O Pan20/data/HOR_10350.fasta.gz
fi

# HOR_13821
if [ ! -f Pan20/data/HOR_13821.fasta.gz ] && [ ! -f Pan20/data/HOR_13821.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=c08f58209a47fe2d&folder_ids=" -O Pan20/data/HOR_13821.fasta.gz
fi

# HOR_13942
if [ ! -f Pan20/data/HOR_13942.fasta.gz ] && [ ! -f Pan20/data/HOR_13942.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=59606d2a36c77a56&folder_ids=" -O Pan20/data/HOR_13942.fasta.gz
fi

# HOR_21599
if [ ! -f Pan20/data/HOR_21599.fasta.gz ] && [ ! -f Pan20/data/HOR_21599.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=3f00d864a8cefbf4&folder_ids=" -O Pan20/data/HOR_21599.fasta.gz
fi

# HOR_3081
if [ ! -f Pan20/data/HOR_3081.fasta.gz ] && [ ! -f Pan20/data/HOR_3081.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=be8c8ac3dbd1dd54&folder_ids=" -O Pan20/data/HOR_3081.fasta.gz
fi

# HOR_3365
if [ ! -f Pan20/data/HOR_3365.fasta.gz ] && [ ! -f Pan20/data/HOR_3365.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=b2486a20bc56b90f&folder_ids=" -O Pan20/data/HOR_3365.fasta.gz
fi

# HOR_7552
if [ ! -f Pan20/data/HOR_7552.fasta.gz ] && [ ! -f Pan20/data/HOR_7552.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=83e633fc7a5c928a&folder_ids=" -O Pan20/data/HOR_7552.fasta.gz
fi

# HOR_8148
if [ ! -f Pan20/data/HOR_8148.fasta.gz ] && [ ! -f Pan20/data/HOR_8148.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=795b4cb2a2107ba7&folder_ids=" -O Pan20/data/HOR_8148.fasta.gz
fi

# HOR_9043
if [ ! -f Pan20/data/HOR_9043.fasta.gz ] && [ ! -f Pan20/data/HOR_9043.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=c86c1b73aa7102dd&folder_ids=" -O Pan20/data/HOR_9043.fasta.gz
fi

# Igri
if [ ! -f Pan20/data/Igri.fasta.gz ] && [ ! -f Pan20/data/Igri.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=1ac3b66e11b8fe50&folder_ids=" -O Pan20/data/Igri.fasta.gz
fi

# OUN333
if [ ! -f Pan20/data/OUN333.fasta.gz ] && [ ! -f Pan20/data/OUN333.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=e2132aef71b11dbf&folder_ids=" -O Pan20/data/OUN333.fasta.gz
fi

# Planet
if [ ! -f Pan20/data/Planet.fasta.gz ] && [ ! -f Pan20/data/Planet.fa ]; then
    wget "https://galaxy-web.ipk-gatersleben.de/api/libraries/datasets/download/uncompressed?ld_ids=0d72ca01c763d02d&folder_ids=" -O Pan20/data/Planet.fasta.gz
fi

# MorexV3
if [ ! -f Pan20/data/MorexV3.fasta.gz ] && [ ! -f Pan20/data/MorexV3.fa ]; then
    wget "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_904849725.1?download=true&gzip=true" -O Pan20/data/MorexV3.fasta.gz
fi



# Check if there are files ending in .fasta.gz in the /Pan20/data directory
if ls Pan20/data/*.fasta.gz 1> /dev/null 2>&1; then
    # Unzip and change the chromosome names for all files in /Pan20/data ending with .fasta.gz
    for genome_file in Pan20/data/*.fasta.gz; do
        output_file="Pan20/data/$(basename "$genome_file" .fasta.gz).fa"
        sample_name=$(basename "$genome_file" .fasta.gz)
        if [ ! -f "$output_file" ]; then
            zcat "$genome_file" | perl -lne 'if(/^>ENA\|[^|]+\|(\S+).*?chromosome: (\S+)/){ print ">chr$2 sampleName='"$sample_name"'" } else { last if(/^>/); print }' > "$output_file"
            # rm "$genome_file" # Commented for debugging purposes
        fi

        # Check if the file was created successfully and is not empty
        if [ -s "$output_file" ]; then
            echo "Processed and created: $output_file"
            continue
        else
        # Check it the names already start with "chr" and has not ENA
            zcat "$genome_file" | perl -lne 'if(/^>chr([0-9]+H)/){print ">chr$1 sampleName='"$sample_name"'"}elsif(/^>/){last}else{print}' > "$output_file"

            # rm "$genome_file" # Commented for debugging purposes
        fi
    done         




else
    echo "No .fasta.gz files found in the Pan20/data/ directory. Skipping this step."
fi

# # Execute prepare assemblies from PHG (but only if there is a file that should exists that does not)    ## Uses 32 threads
# for genome_file in "${genome_files[@]}"; do
#    output_file="Pan20/data/prepared_assemblies/$(basename "$genome_file" .fasta.gz).fa"
#    if [ ! -f "$output_file" ]; then
#        echo "phg prepare_assemblies --keyfile Pan20/data/prepare_assemblies.txt --output-dir Pan20/data/prepared_assemblies/ --threads 32"
#        phg prepare-assemblies --keyfile Pan20/data/prepare_assemblies.txt --output-dir Pan20/data/prepared_assemblies/ --threads 32
#        break
#    fi
#done
#########   EVERYTHING COMMENTED, BECAUSE THE FILES ARE ALREADY PREPARED AT THIS VERSION OF THE SCRIPT ######################


### Download Morex.gff

if [ ! -f Pan20/data/MorexV3.gff ]; then
    wget "https://doi.ipk-gatersleben.de/DOI/b2f47dfb-47ff-4114-89ae-bad8dcc515a1/5d16cc17-c37f-417f-855d-c5e72c721f6c/1/DOWNLOAD" -O Pan20/data/RAW_MorexV3.gff
fi
# Not necessary to modify the gff due chromosomes names: already match with the prepared fasta files
# But remove all is not in 7 chromosomes

if [ ! -f Pan20/data/MorexV3.gff ]; then
    perl -ne 'if(/^chr[1-7]H\t/){print}' Pan20/data/RAW_MorexV3.gff > Pan20/data/MorexV3.gff
    rm Pan20/data/RAW_MorexV3.gff
fi

# Prepare reference-ranges
if [ ! -f Pan20/output/ref_ranges.bed ]; then
    phg create-ranges --gff Pan20/data/MorexV3.gff --boundary gene -o Pan20/output/ref_ranges.bed --reference-file Pan20/data/MorexV3.fa
fi


# if exists the keyfile for the aligment, delete it
if [ -f Pan20/data/alignment_keyfile.txt ]; then
    rm Pan20/data/alignment_keyfile.txt
fi

# Create the keyfile for the alignment
# Check if the .fa files exist in the Pan20/data/ directory
# If they do, create the keyfile, but not "creating" a "*.fa" file
if ls Pan20/data/*.fa 1> /dev/null 2>&1; then
    for prepared_file in Pan20/data/*.fa; do
        if [[ "$(basename "$prepared_file")" != "MorexV3.fa" ]]; then
            echo "$prepared_file" >> Pan20/data/alignment_keyfile.txt
        fi
    done
else
    echo "No .fa files found in the Pan20/data/ directory."
fi

# Save in a variable the genome files
genome_files=(Pan20/data/*.fa)

# Align assemblies
# Execute align_assemblies from PHG (but only if there is a file that should exists that does not)    ## Uses 32 threads
for genome_file in "${genome_files[@]}"; do
    if [[ "$(basename "$genome_file" .fa)" == "MorexV3" ]]; then
        continue
    fi
    if [ ! -f "Pan20/output/alignment_files/$(basename "$genome_file" .fa).maf" ]; then
        phg align-assemblies --gff Pan20/data/MorexV3.gff -o Pan20/output/alignment_files/ --total-threads 32 --in-parallel 2 --assembly-file-list Pan20/data/alignment_keyfile.txt --output-dir Pan20/output/alignment_files/ --reference-file Pan20/data/MorexV3.fa
        break
    fi
done


# Checkpoint. Script only continues if the alignment files are created
for genome_file in "${genome_files[@]}"; do
    if [[ "$(basename "$genome_file" .fa)" == "MorexV3" ]]; then
        continue
    fi
    if [ ! -f "Pan20/output/alignment_files/$(basename "$genome_file" .fa).maf" ]; then
        echo "Alignment files not created for $(basename "$genome_file" .fa). Exiting script."
        exit 1
    fi
done

# Build the agc file to store the genomes
if [ ! -f Pan20/vcf_dbs/assemblies.agc ]; then
    phg agc-compress --db-path Pan20/vcf_dbs/ --fasta-list Pan20/data/alignment_keyfile.txt --reference-file Pan20/data/MorexV3.fa
fi

# Create the reference VCF
if [ ! -f Pan20/vcf_dbs/hvcf_files/MorexV3.h.vcf.gz ]; then
    phg create-ref-vcf --bed Pan20/output/ref_ranges.bed --reference-file Pan20/data/MorexV3.fa --reference-name MorexV3 --db-path Pan20/vcf_dbs/
fi

# Create from allignments the VCFs
# Execute create-maf-vcf from PHG (but only if there is a file that should exists that does not)    ## Uses 32 threads
for genome_file in "${genome_files[@]}"; do
    if [[ "$(basename "$genome_file" .fa)" == "MorexV3" ]]; then
        continue
    fi
    if [ ! -f "Pan20/output/vcf_files/$(basename "$genome_file" .fa).h.vcf.gz" ]; then
        phg create-maf-vcf --bed Pan20/output/ref_ranges.bed --reference-file Pan20/data/MorexV3.fa -o Pan20/output/vcf_files/ --db-path Pan20/vcf_dbs/ --skip-metrics --maf-dir Pan20/output/alignment_files/
    fi
done

# Set a checkpoint to see if the files are already created
for genome_file in "${genome_files[@]}"; do
    output_file="Pan20/data/$(basename "$genome_file" .fa).fa"

    if [ ! -f "Pan20/vcf_dbs/hvcf_files/$(basename "$genome_file" .fa).h.vcf.gz" ]; then
        phg load-vcf --db-path Pan20/vcf_dbs/ --vcf-dir Pan20/output/vcf_files/ --threads 32
        break
    fi
done


# Generate necessary files: hapIDranges:
if [ ! -f Pan20/output/hapIDranges.tsv ]; then
    phg sample-hapid-by-range --input-dir Pan20/vcf_dbs/hvcf_files/ --output-file Pan20/output/hapIDranges.tsv
fi

# Generate necessary files: hapIDranges:
if [ ! -f Pan20/output/hapIDtable.tsv ]; then
    phg hapid-sample-table --hvcf-dir Pan20/vcf_dbs/hvcf_files/ --output-file Pan20/output/hapIDtable.tsv
fi

# Generate the gmap database
for genome_file in "${genome_files[@]}"; do
    genome_name=$(basename "$genome_file" .fa)
    if [ ! -f "gmap_db/${genome_name}/${genome_name}.chromosome" ]; then
        $gmap_build_command -D gmap_db/ -d "${genome_name}" "Pan20/data/${genome_name}.fa"
    fi
done

############
# Skipping this part, is being deprecated in the new versions of phg
### In order to create the kmer index, we have had memmory issues in the past with java
### Ensure that there is enough memory
### In the past, 200gb where not enough, so its enough with 300gb
### Check how much memory is available

#if [[ "$JAVA_OPTS" != *"-Xmx300g"* ]]; then
#   export JAVA_OPTS="-Xmx300g"
#fi

#echo "Java memory set to:"
#echo $JAVA_OPTS

#if [ ! -f Pan20/output/kmerIndex.txt ]; then
#    phg build-kmer-index --db-path Pan20/vcf_dbs/ --index-file Pan20/output/kmerIndex.txt --hvcf-dir Pan20/vcf_dbs/hvcf_files/ --use-big-discard-set
#fi
############

### Preparing kmer mapping in BETA
if [ ! -f Pan20/output/ropeBWT_index.fmd ]; then
    phg rope-bwt-index --hvcf-dir Pan20/vcf_dbs/hvcf_files/ --db-path Pan20/vcf_dbs/ --output-dir Pan20/output/ --index-file-prefix ropeBWT_index --threads 32
fi

#######
#######
# If you want to align to the pangenome using barleymap, is necessary to have a bed file of each genome haplotypes:
#######
#######

# Create, if not exists, the python script that process the hvcf files into bed files

if [ ! -f hvcf2bed.py ]; then
    echo "Creating hvcf2bed.py script"

    hvcf2bed_script="hvcf2bed.py"

    cat << 'EOF' > "$hvcf2bed_script"
import sys
import gzip
import os
import re

def parse_vcf_lines(lines, genome_name):
    viewed_checksums = []
    empty_ranges = 0
    skipped_checksums_bug = 0
    output_lines = []
    prev_checksum = None  # Track previous checksum

    for lineno, line in enumerate(lines, start=1):
        line = line.strip()
        
        # Match single-segment ALT line
        match1 = re.search(
            r'SampleName=([^,]+),Regions=([^:]+):(\d+)-(\d+),Checksum=([^,]+),RefChecksum=([^,]+),RefRange=([^:]+):(\d+)-(\d+)',
            line
        )
        if match1:
            genome, chr_, start, end, checksum, ref_checksum, ref_chr, ref_start, ref_end = match1.groups()
            strand = "+"
            start, end = int(start), int(end)
            if start > end:
                strand = "-"
                start, end = end, start

            length = end - start

            # Check duplicate checksum & short length
            if checksum == prev_checksum and length <= 1:
                print(
                    f"WARNING (line {lineno}): duplicate checksum {checksum} with length {length} skipped",
                    file=sys.stderr
                )
                continue
            elif checksum in viewed_checksums: 
                print(
                    f"WARNING (line {lineno}): checksum {checksum} already seen, possible data issue",
                    file=sys.stderr
                )
            elif length <= 1:
                empty_ranges += 1
                continue
            else:
                output_lines.append(
                    f"{chr_}\t{start}\t{end}\t{strand}\t{checksum}\t{genome}\t{ref_chr}\t{ref_start}\t{ref_end}\t{ref_checksum}"
                )
                prev_checksum = checksum
                viewed_checksums.append(checksum)
                continue

        # Match multi-segment ALT line
        match2 = re.search(
            r'SampleName=([^,]+),Regions="([^"]+)",Checksum=([^,]+),RefChecksum=([^,]+),RefRange=([^:]+):(\d+)-(\d+)',
            line
        )
        if match2:
            genome, regions, checksum, ref_checksum, ref_chr, ref_start, ref_end = match2.groups()
            for segment in regions.split(","):
                m = re.match(r'([^:]+):(\d+)-(\d+)', segment)
                if m:
                    chr_, start, end = m.groups()
                    strand = "+"
                    start, end = int(start), int(end)
                    if start > end:
                        strand = "-"
                        start, end = end, start

                    length = end - start

                    # Check duplicate checksum & short length
                    if checksum == prev_checksum and length <= 1:
                        skipped_checksums_bug = skipped_checksums_bug + 1  
                        continue    #skip this duplicate

                    output_lines.append(
                        f"{chr_}\t{start}\t{end}\t{strand}\t{ref_checksum}\t{genome}\t{ref_chr}\t{ref_start}\t{ref_end}\t{checksum}"
                    )
                    prev_checksum = checksum
    print(f"Skipped {skipped_checksums_bug} duplicate checksums with length <= 1")
    print(f"Skipped {empty_ranges} empty ranges")

    return output_lines

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 hvcf2bed.py <vcf_folder> <genome_name>")
        sys.exit(1)

    vcf_folder = sys.argv[1]
    genome_name = sys.argv[2]
    vcf_file = os.path.join(vcf_folder, f"{genome_name}.h.vcf.gz")
    bed_file = os.path.join(vcf_folder, f"{genome_name}.h.bed")

    if not os.path.exists(vcf_file):
        print(f"VCF file not found: {vcf_file}")
        sys.exit(1)

    with gzip.open(vcf_file, "rt") as f:
        output = parse_vcf_lines(f, genome_name)

    # Sort output lines by chromosome and start position
    output_sorted = sorted(
        output,
        key=lambda x: (x.split("\t")[0], int(x.split("\t")[1]))
    )

    with open(bed_file, "w") as out:
        out.write("\n".join(output_sorted) + "\n")

    print(f"Converted {vcf_file} to {bed_file}")
EOF
fi

# Create the bed files for each genome
for genome_file in "${genome_files[@]}"; do
    genome_name=$(basename "$genome_file" .fa)
    if [ ! -f "Pan20/vcf_dbs/hvcf_files/${genome_name}.h.bed" ]; then
        echo "Creating bed file for $genome_name"
        python3 hvcf2bed.py Pan20/vcf_dbs/hvcf_files/ "$genome_name"
    fi
done

# Index the fasta files
for genome_file in "${genome_files[@]}"; do
    genome_name=$(basename "$genome_file" .fa)
    if [ ! -f "Pan20/data/${genome_name}.fa.fai" ]; then
        echo "Indexing fasta file for $genome_name"
        samtools faidx "Pan20/data/${genome_name}.fa"
    fi
done

# If align2graph.py is not present, download it
if [ ! -f align2graph.py ]; then
    echo "Downloading align2graph.py script"
    wget "https://raw.githubusercontent.com/eead-csic-compbio/eead-csic-compbio.github.io/master/scripts/align2graph.py" -O align2graph.py
fi
