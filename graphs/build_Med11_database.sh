# Date: 26.03.2025

# Works with PHGv2 
# citation: 10.1093/bioinformatics/btac410
# Used in version phg version 2.4.52.207

# Last update: 07.04.2025
# Time for a full run: Approx 3 days
# Space used arround 500gb nowadays

# Update: 28/07/2025
# Updated hvcf2bed.py script to addapt to barleymap structure
# Use gmap version installed into esplandian, using the version of gmap 2013-08-31 (less storage, and same that barleymap version)

# Update 31/07/2025
# Added GDB_136 sequenced and ensambled at IPK
# Increases java heap space to 350, explotes with GDB136 included

# Update 04/08/2025
# hvcf2bed.py script now sorts the output by chromosome and start position

# Update 23/09/2025
# hvcf2bed.py script now warns if a checksum is duplicated in the same vcf file, and delete very short ranges (length <=1)
# Commented out build of old kmer index, because is being deprecated in the new versions of phg

###
# _________________________________________________________________________________________________
###

# Remember to run with proper conda environment
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
mkdir -p Med11/data
mkdir -p Med11/output
mkdir -p Med11/vcf_dbs
mkdir -p Med11/output/alignment_files
mkdir -p Med11/output/vcf_files
mkdir -p Med11/output/read_mapping
mkdir -p Med11/output/imputed_vcf_files
mkdir -p gmap_db

# Initialize the database
phg initdb --db-path Med11/vcf_dbs/

# Get the genomes
# from the folder:
# https://galaxy-web.ipk-gatersleben.de/libraries/folders/Fd071e794759ab192
# citation: https://doi.org/10.1038/s41586-024-08187-1

# Download the genomes
# https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783365.1?download=true&gzip=true   HOR_14121
# https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783375.1?download=true&gzip=true   HOR_21595
# https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783475.1?download=true&gzip=true   HOR_3474
# https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783295.1?download=true&gzip=true   HOR_13942
# https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783435.1?download=true&gzip=true   HOR_21599
# https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783415.1?download=true&gzip=true   HOR_3365
# https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783545.1?download=true&gzip=true   HOR_2830
# https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783585.1?download=true&gzip=true   HOR_10892
# https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783485.1?download=true&gzip=true   HOR_2779
# https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783555.1?download=true&gzip=true   HOR_1168
# https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783505.1?download=true&gzip=true   HOR_12184
# https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_904849725.1?download=true&gzip=true   MorexV3

if [ ! -f Med11/data/HOR_14121.fasta.gz ] && [ ! -f Med11/data/HOR_14121.fa ]; then
    wget "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783365.1?download=true&gzip=true" -O Med11/data/HOR_14121.fasta.gz
fi
if [ ! -f Med11/data/HOR_21595.fasta.gz ] && [ ! -f Med11/data/HOR_21595.fa ]; then
    wget "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783375.1?download=true&gzip=true" -O Med11/data/HOR_21595.fasta.gz
fi
if [ ! -f Med11/data/HOR_3474.fasta.gz ] && [ ! -f Med11/data/HOR_3474.fa ]; then
    wget "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783475.1?download=true&gzip=true" -O Med11/data/HOR_3474.fasta.gz
fi
if [ ! -f Med11/data/HOR_13942.fasta.gz ] && [ ! -f Med11/data/HOR_13942.fa ]; then
    wget "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783295.1?download=true&gzip=true" -O Med11/data/HOR_13942.fasta.gz
fi
if [ ! -f Med11/data/HOR_21599.fasta.gz ] && [ ! -f Med11/data/HOR_21599.fa ]; then
    wget "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783435.1?download=true&gzip=true" -O Med11/data/HOR_21599.fasta.gz
fi
if [ ! -f Med11/data/HOR_3365.fasta.gz ] && [ ! -f Med11/data/HOR_3365.fa ]; then
    wget "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783415.1?download=true&gzip=true" -O Med11/data/HOR_3365.fasta.gz
fi
if [ ! -f Med11/data/HOR_2830.fasta.gz ] && [ ! -f Med11/data/HOR_2830.fa ]; then
    wget "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783545.1?download=true&gzip=true" -O Med11/data/HOR_2830.fasta.gz
fi
if [ ! -f Med11/data/HOR_10892.fasta.gz ] && [ ! -f Med11/data/HOR_10892.fa ]; then
    wget "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783585.1?download=true&gzip=true" -O Med11/data/HOR_10892.fasta.gz
fi
if [ ! -f Med11/data/HOR_2779.fasta.gz ] && [ ! -f Med11/data/HOR_2779.fa ]; then
    wget "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783485.1?download=true&gzip=true" -O Med11/data/HOR_2779.fasta.gz
fi
if [ ! -f Med11/data/HOR_1168.fasta.gz ] && [ ! -f Med11/data/HOR_1168.fa ]; then
    wget "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783555.1?download=true&gzip=true" -O Med11/data/HOR_1168.fasta.gz
fi
if [ ! -f Med11/data/HOR_12184.fasta.gz ] && [ ! -f Med11/data/HOR_12184.fa ]; then
    wget "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_949783505.1?download=true&gzip=true" -O Med11/data/HOR_12184.fasta.gz
fi
if [ ! -f Med11/data/MorexV3.fasta.gz ] && [ ! -f Med11/data/MorexV3.fa ]; then
    wget "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_904849725.1?download=true&gzip=true" -O Med11/data/MorexV3.fasta.gz
fi

## From the own assembly of GDB136, which is not yet published:

#    scp /scratch/GDB136/IPK/output/pseudomolecules_v1/pseudomolecules_v1.fasta Med11/data/GDB_136_raw.fa
#    cat Med11/data/GDB_136_raw.fa | perl -lne 'if (/^>GDB\d+_(chr\S+)/) { print ">$1" } elsif (/^>/) { last } else { print }' > Med11/data/GDB_136.fa
#    rm Med11/data/GDB_136_raw.fa
#    mv Med11/data/GDB_136.fa Med11/data/GDB_136_raw.fa
#    rm Med11/data/prepare_assemblies.txt
#    echo -e "Med11/data/GDB_136_raw.fa\tGDB_136" > Med11/data/prepare_assemblies.txt
#    phg prepare-assemblies --keyfile Med11/data/prepare_assemblies.txt --output-dir Med11/data/ --threads 32

# Check if there are files ending in .fasta.gz in the /data/ directory
if ls Med11/data/*.fasta.gz 1> /dev/null 2>&1; then
    # Unzip and change the chromosome names for all files in /data/ ending with .fasta.gz
    for genome_file in Med11/data/*.fasta.gz; do
        output_file="Med11/data/$(basename "$genome_file" .fasta.gz).fa"
        sample_name=$(basename "$genome_file" .fasta.gz)
        if [ ! -f "$output_file" ]; then
            zcat "$genome_file" | perl -lne 'if(/^>ENA\|[^|]+\|(\S+).*?chromosome: (\S+)/){ print ">chr$2 sampleName='"$sample_name"'" } else { last if(/^>/); print }' > "$output_file"
            # rm "$genome_file" # Commented for debugging purposes
        fi
    done
else
    echo "No .fasta.gz files found in the Med11/data/ directory. Skipping this step."
fi


# # Execute prepare assemblies from PHG (but only if there is a file that should exists that does not)    ## Uses 32 threads
# for genome_file in "${genome_files[@]}"; do
#    output_file="Med11/data/prepared_assemblies/$(basename "$genome_file" .fasta.gz).fa"
#    if [ ! -f "$output_file" ]; then
#        echo "phg prepare_assemblies --keyfile Med11/data/prepare_assemblies.txt --output-dir Med11/data/prepared_assemblies/ --threads 32"
#        phg prepare-assemblies --keyfile Med11/data/prepare_assemblies.txt --output-dir Med11/data/prepared_assemblies/ --threads 32
#        break
#    fi
#done
#########   EVERYTHING COMMENTED, BECAUSE THE FILES ARE ALREADY PREPARED AT THIS VERSION OF THE SCRIPT ######################


### Download Morex.gff

if [ ! -f Med11/data/MorexV3.gff ]; then
    wget "https://doi.ipk-gatersleben.de/DOI/b2f47dfb-47ff-4114-89ae-bad8dcc515a1/5d16cc17-c37f-417f-855d-c5e72c721f6c/1/DOWNLOAD" -O Med11/data/RAW_MorexV3.gff
fi
# Not necessary to modify the gff due chromosomes names: already match with the prepared fasta files
# But remove all is not in 7 chromosomes

if [ ! -f Med11/data/MorexV3.gff ]; then
    perl -ne 'if(/^chr[1-7]H\t/){print}' Med11/data/RAW_MorexV3.gff > Med11/data/MorexV3.gff
    rm Med11/data/RAW_MorexV3.gff
fi

# Prepare reference-ranges
if [ ! -f Med11/output/ref_ranges.bed ]; then
    phg create-ranges --gff Med11/data/MorexV3.gff --boundary gene -o Med11/output/ref_ranges.bed --reference-file Med11/data/MorexV3.fa
fi


# if exists the keyfile for the aligment, delete it
if [ -f Med11/data/alignment_keyfile.txt ]; then
    rm Med11/data/alignment_keyfile.txt
fi

# Create the keyfile for the alignment
# Check if the .fa files exist in the Med11/data/ directory
# If they do, create the keyfile, but not "creating" a "*.fa" file
if ls Med11/data/*.fa 1> /dev/null 2>&1; then
    for prepared_file in Med11/data/*.fa; do
        if [[ "$(basename "$prepared_file")" != "MorexV3.fa" ]]; then
            echo "$prepared_file" >> Med11/data/alignment_keyfile.txt
        fi
    done
else
    echo "No .fa files found in the Med11/data/ directory."
fi

# Save in a variable the genome files
genome_files=(Med11/data/*.fa)

# Align assemblies
# Execute align_assemblies from PHG (but only if there is a file that should exists that does not)    ## Uses 32 threads
for genome_file in "${genome_files[@]}"; do
    if [[ "$(basename "$genome_file" .fa)" == "MorexV3" ]]; then
        continue
    fi
    if [ ! -f "Med11/output/alignment_files/$(basename "$genome_file" .fa).maf" ]; then
        phg align-assemblies --gff Med11/data/MorexV3.gff -o Med11/output/alignment_files/ --total-threads 32 --in-parallel 2 --assembly-file-list Med11/data/alignment_keyfile.txt --output-dir Med11/output/alignment_files/ --reference-file Med11/data/MorexV3.fa
        break
    fi
done


# Checkpoint. Script only continues if the alignment files are created
for genome_file in "${genome_files[@]}"; do
    if [[ "$(basename "$genome_file" .fa)" == "MorexV3" ]]; then
        continue
    fi
    if [ ! -f "Med11/output/alignment_files/$(basename "$genome_file" .fa).maf" ]; then
        echo "Alignment files not created for $(basename "$genome_file" .fa). Exiting script."
        exit 1
    fi
done

# Build the agc file to store the genomes
if [ ! -f Med11/vcf_dbs/assemblies.agc ]; then
    phg agc-compress --db-path Med11/vcf_dbs/ --fasta-list Med11/data/alignment_keyfile.txt --reference-file Med11/data/MorexV3.fa
fi

# Create the reference VCF
if [ ! -f Med11/vcf_dbs/hvcf_files/MorexV3.h.vcf.gz ]; then
    phg create-ref-vcf --bed Med11/output/ref_ranges.bed --reference-file Med11/data/MorexV3.fa --reference-name MorexV3 --db-path Med11/vcf_dbs/
fi

# Create from allignments the VCFs
# Execute create-maf-vcf from PHG (but only if there is a file that should exists that does not)    ## Uses 32 threads
for genome_file in "${genome_files[@]}"; do
    if [[ "$(basename "$genome_file" .fa)" == "MorexV3" ]]; then
        continue
    fi
    if [ ! -f "Med11/output/vcf_files/$(basename "$genome_file" .fa).h.vcf.gz" ]; then
        phg create-maf-vcf --bed Med11/output/ref_ranges.bed --reference-file Med11/data/MorexV3.fa -o Med11/output/vcf_files/ --db-path Med11/vcf_dbs/ --skip-metrics --maf-dir Med11/output/alignment_files/
    fi
done

# Set a checkpoint to see if the files are already created
for genome_file in "${genome_files[@]}"; do
    output_file="Med11/data/$(basename "$genome_file" .fa).fa"

    if [ ! -f "Med11/vcf_dbs/hvcf_files/$(basename "$genome_file" .fa).h.vcf.gz" ]; then
        phg load-vcf --db-path Med11/vcf_dbs/ --vcf-dir Med11/output/vcf_files/ --threads 32
        break
    fi
done


# Generate necessary files: hapIDranges:
if [ ! -f Med11/output/hapIDranges.tsv ]; then
    phg sample-hapid-by-range --input-dir Med11/vcf_dbs/hvcf_files/ --output-file Med11/output/hapIDranges.tsv
fi

# Generate necessary files: hapIDranges:
if [ ! -f Med11/output/hapIDtable.tsv ]; then
    phg hapid-sample-table --hvcf-dir Med11/vcf_dbs/hvcf_files/ --output-file Med11/output/hapIDtable.tsv
fi

# Generate the gmap database
for genome_file in "${genome_files[@]}"; do
    genome_name=$(basename "$genome_file" .fa)
    if [ ! -f "gmap_db/${genome_name}/${genome_name}.chromosome" ]; then
        $gmap_build_command -D gmap_db/ -d "${genome_name}" "Med11/data/${genome_name}.fa"
    fi
done

########
# Skipping this part, is being deprecated in the new versions of phg
### In order to create the kmer index, we have had memmory issues in the past with java
### Ensure that there is enough memory
### In the past, 200gb where not enough, so its enough with 300gb
### Check how much memory is available

#if [[ "$JAVA_OPTS" != *"-Xmx350g"* ]]; then
#    export JAVA_OPTS="-Xmx350g"
#fi

#echo "Java memory set to:"
#echo $JAVA_OPTS

#if [ ! -f Med11/outputkmerIndex.txt ]; then
#    phg build-kmer-index --db-path Med11/vcf_dbs/ --index-file Med11/outputkmerIndex.txt --hvcf-dir Med11/vcf_dbs/hvcf_files/ --use-big-discard-set
#fi
########

### Preparing kmer mapping in BETA
if [ ! -f Med11/output/ropeBWT_index.fmd ]; then
    phg rope-bwt-index --hvcf-dir Med11/vcf_dbs/hvcf_files/ --db-path Med11/vcf_dbs/ --output-dir Med11/output --index-file-prefix ropeBWT_index --threads 32
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
    output_lines = []

    for line in lines:
        line = line.strip()
        
        # Match single-segment ALT line
        match1 = re.search(r'SampleName=([^,]+),Regions=([^:]+):(\d+)-(\d+),Checksum=([^,]+),RefChecksum=([^,]+),RefRange=([^:]+):(\d+)-(\d+)', line)
        if match1:
            genome, chr_, start, end, checksum, ref_checksum, ref_chr, ref_start, ref_end = match1.groups()
            strand = "+"
            start, end = int(start), int(end)
            if start > end:
                strand = "-"
                start, end = end, start

            output_lines.append(f"{chr_}\t{start}\t{end}\t{strand}\t{checksum}\t{genome}\t{ref_chr}\t{ref_start}\t{ref_end}\t{ref_checksum}")
            continue

        # Match multi-segment ALT line
        match2 = re.search(r'SampleName=([^,]+),Regions="([^"]+)",Checksum=([^,]+),RefChecksum=([^,]+),RefRange=([^:]+):(\d+)-(\d+)', line)
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

                    output_lines.append(f"{chr_}\t{start}\t{end}\t{strand}\t{ref_checksum}\t{genome}\t{ref_chr}\t{ref_start}\t{ref_end}\t{checksum}")

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
EOF
fi

# Create the bed files for each genome
for genome_file in "${genome_files[@]}"; do
    genome_name=$(basename "$genome_file" .fa)
    if [ ! -f "Med11/vcf_dbs/hvcf_files/${genome_name}.h.bed" ]; then
        echo "Creating bed file for $genome_name"
        python3 hvcf2bed.py Med11/vcf_dbs/hvcf_files/ "$genome_name"
    fi
done

# Index the fasta files
for genome_file in "${genome_files[@]}"; do
    genome_name=$(basename "$genome_file" .fa)
    if [ ! -f "Med11/data/${genome_name}.fa.fai" ]; then
        echo "Indexing fasta file for $genome_name"
        samtools faidx "Med11/data/${genome_name}.fa"
    fi
done

# Prepare to work with align2graph

# If align2graph.py is not present, download it
if [ ! -f align2graph.py ]; then
    echo "Downloading align2graph.py script"
    wget "https://raw.githubusercontent.com/eead-csic-compbio/eead-csic-compbio.github.io/master/scripts/align2graph.py" -O align2graph.py
fi

