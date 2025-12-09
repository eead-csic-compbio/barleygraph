import pandas as pd
import os
import time
import sys
import gzip
import argparse
import subprocess
from tqdm import tqdm
import warnings

# Visualization Imports
import matplotlib
matplotlib.use('Agg') # Key for running on servers without a display
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#######################################################################################################
#                                      HELPER FUNCTIONS
#######################################################################################################

def extract_chr_list(bed_file, genotype_group_sort_hash):
    chr_list = []
    reference_genotype = next(genotype for genotype, info in genotype_group_sort_hash.items() if info['Group'] == 'Reference')
    
    # Handle filename variations
    base_name = os.path.basename(bed_file)
    dir_name = os.path.dirname(bed_file)
    
    # Construct reference filename pattern
    # Assuming bed_file format is {genotype}.h.bed
    current_genotype = next(iter(genotype_group_sort_hash))
    reference_bed_file = bed_file.replace(current_genotype, reference_genotype)

    if not os.path.exists(reference_bed_file):
        # Fallback: try to construct it manually if replace failed
        reference_bed_file = os.path.join(dir_name, f"{reference_genotype}.h.bed")
        if not os.path.exists(reference_bed_file):
             raise FileNotFoundError(f"Reference BED file not found: {reference_bed_file}")
    
    with open(reference_bed_file, 'r') as bed_f:
        header_skipped = False
        for line in bed_f:
            if line.startswith('#'): continue
            if not header_skipped:
                header_skipped = True
                continue
            
            chr_name = line.split('\t')[0]
            if chr_name not in chr_list:
                chr_list.append(chr_name)
            
    return chr_list

def read_and_validate_grouped_file(grouped_file, verbose=False):
    try:
        df_grouped_validate = pd.read_csv(grouped_file, sep=r'\s+', engine='python')
        # Strip whitespace
        df_grouped_validate = df_grouped_validate.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
        
        required_cols = ['Sort', 'Genotype', 'Group']
        if not all(col in df_grouped_validate.columns for col in required_cols):
            raise ValueError(f"Grouped file must contain columns: {required_cols}")
        
        duplicates = df_grouped_validate[df_grouped_validate['Genotype'].duplicated(keep=False)]
        if not duplicates.empty:
            raise ValueError(f"Duplicate genotypes found: {duplicates['Genotype'].unique().tolist()}")
        
        reference_count = (df_grouped_validate['Group'] == 'Reference').sum()
        if reference_count == 0:
            raise ValueError("Grouped file must contain at least one sample with Group='Reference'")
        
        if verbose:
            print(f"Found {len(df_grouped_validate)} unique genotypes")
        
        genotype_group_sort_hash = df_grouped_validate.set_index('Genotype')[['Group', 'Sort']].to_dict(orient='index')

    except Exception as e:
        print(f"\nERROR reading grouped file: {grouped_file}")
        print("Expected format: Sort <tab> Genotype <tab> Group")
        raise ValueError(f"Error reading grouped file: {e}")
    
    return genotype_group_sort_hash

def add_colors_to_hapIDs(genotype_group_sort_hash):
    base_colors = [
        '#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffd92f', 
        '#a65628', '#f781bf', '#999999', '#00ced1', '#000075', '#ffb300', 
        '#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', 
        '#e5c494', '#b3b3b3'
    ]
    
    hapID_color_map = {}
    color_idx = 0
    for genotype, info in genotype_group_sort_hash.items():
        if info['Group'] == 'Reference':
            hapID_color_map[genotype] = '#000000'
        elif info['Group'] in ['Pangenome', 'Imputed']:
            hapID_color_map[genotype] = base_colors[color_idx % len(base_colors)]
            color_idx += 1

    return hapID_color_map

def find_hvcf_files(vcf_folder, genotype_group_sort_hash):
    hvcf_files = {}
    for genotype in genotype_group_sort_hash.keys():
        vcf_file_gz = os.path.join(vcf_folder, f"{genotype}.h.vcf.gz")
        vcf_file_plain = os.path.join(vcf_folder, f"{genotype}.h.vcf")
        
        if os.path.exists(vcf_file_gz):
            hvcf_files[genotype] = vcf_file_gz
        elif os.path.exists(vcf_file_plain):
            hvcf_files[genotype] = vcf_file_plain
        else:
            raise FileNotFoundError(f"VCF file not found for genotype: {genotype} in {vcf_folder}")
    
    return hvcf_files

def hvcf_to_bed(genotype_group_sort_hash, hvcf_files, verbose=False):
    genotypes = genotype_group_sort_hash.keys()
    bed_files = {}
    
    for genotype in genotypes:
        vcf_file = hvcf_files[genotype]
        # Handle both .gz and plain for output name
        if vcf_file.endswith('.h.vcf.gz'):
            bed_file = vcf_file.replace('.h.vcf.gz', '.h.bed')
        else:
            bed_file = vcf_file.replace('.h.vcf', '.h.bed')
        
        # Check if BED file already exists and is valid
        if os.path.exists(bed_file):
            try:
                open_func = gzip.open if bed_file.endswith('.gz') else open
                mode = 'rt' if bed_file.endswith('.gz') else 'r'
                
                with open_func(bed_file, mode) as f:
                    valid_format = True
                    lines_checked = 0
                    for line in f:
                        if line.startswith('#'): continue
                        parts = line.strip().split('\t')
                        if len(parts) != 10:
                            valid_format = False
                            break
                        lines_checked += 1
                        if lines_checked >= 3: break
                    
                    if valid_format and lines_checked >= 1: # At least one line is fine
                        if verbose: print(f"INFO: Reusing existing BED file: {bed_file}")
                        bed_files[genotype] = bed_file
                        continue
            except Exception:
                pass # Fall through to regeneration
            
            if verbose: print(f"Regenerating invalid BED: {bed_file}")
            os.remove(bed_file)

        if verbose: print(f"Converting {vcf_file} to BED...")
        
        # Call hvcf2bed using subprocess
        try:
            folder = os.path.dirname(vcf_file)
            cmd = ['hvcf2bed', folder, genotype]
            if verbose:
                cmd.append('-v')
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            if verbose:
                print(result.stdout)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Error converting {genotype}: {e.stderr}")

        bed_files[genotype] = bed_file

    return bed_files

def generate_df(bed_files, hapID_color_map, genotype_group_sort_hash, hvcfs_folder, verbose=False):
    data_rows = []
    
    donor_color_map = {}
    for genotype, info in genotype_group_sort_hash.items():
        if info['Group'] in ['Pangenome', 'Reference']:
            donor_color_map[genotype] = hapID_color_map[genotype]
    
    for genotype, bed_file in tqdm(bed_files.items(), desc="Processing BED files", disable=not verbose):
        if not os.path.exists(bed_file): continue
        
        group_status = genotype_group_sort_hash[genotype]['Group']
        open_func = gzip.open if bed_file.endswith('.gz') else open
        mode = 'rt' if bed_file.endswith('.gz') else 'r'
        
        with open_func(bed_file, mode) as f:
            header_skipped = False
            for line in f:
                if line.startswith('#'): continue
                if not header_skipped:
                    header_skipped = True
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) != 10: continue
                
                chrom = parts[6]
                ref_start = int(parts[7])
                ref_end = int(parts[8])
                donor_genome = parts[5]
                
                if group_status == 'Reference':
                    color = hapID_color_map[genotype]
                elif group_status == 'Pangenome':
                    color = hapID_color_map[genotype]
                elif group_status == 'Imputed':
                    color = donor_color_map.get(donor_genome, '#CCCCCC')
                else:
                    color = '#CCCCCC'
                
                data_rows.append({
                    'Genotype': genotype,
                    'chr': chrom,
                    'ref_start': ref_start,
                    'ref_end': ref_end,
                    'donor_genome': donor_genome,
                    'color': color
                })
    
    df = pd.DataFrame(data_rows)
    
    if verbose:
        output_dir = os.path.join(hvcfs_folder, "plots")
        os.makedirs(output_dir, exist_ok=True)
        debug_tsv = os.path.join(output_dir, "haplotype_painting_debug_dataframe.tsv")
        df.to_csv(debug_tsv, sep='\t', index=False)
        print(f"Debug dataframe saved to: {debug_tsv}")
    
    return df

def plot_haplotype_painting(hvcfs_folder, chr_to_plot, chr_list, hapID_color_map, genotype_group_sort_hash, df, region_start, region_end, plot_pangenomes, verbose=False):
    base_output_dir = os.path.join(hvcfs_folder, "plots")
    os.makedirs(base_output_dir, exist_ok=True)

    chromosomes_to_process = chr_to_plot if chr_to_plot else chr_list

    for chrom in chromosomes_to_process:
        if region_start and region_end:
            chrom_output_png = os.path.join(base_output_dir, f"{chrom}_{region_start}-{region_end}_haplotype_painting.png")
        else:
            chrom_output_png = os.path.join(base_output_dir, f"{chrom}_FULL_haplotype_painting.png")

        start_time_chr = time.time()
        cdf = df[df['chr'] == chrom].copy()

        if region_start and region_end:
            cdf = cdf[(cdf['ref_start'] >= region_start) & (cdf['ref_end'] <= region_end)].copy()
        
        if not plot_pangenomes:
            genotypes_to_exclude = [g for g, info in genotype_group_sort_hash.items() if info['Group'] == 'Pangenome']
            cdf = cdf[~cdf['Genotype'].isin(genotypes_to_exclude)].copy()
        
        if cdf.empty:
            print(f"No data for {chrom}, skipping.")
            continue
            
        cdf['ref_start_Mbp'] = cdf['ref_start'] / 1e6
        cdf['ref_end_Mbp'] = cdf['ref_end'] / 1e6
        
        # Categorical sorting
        cdf['Genotype'] = pd.Categorical(cdf['Genotype'], categories=sorted(hapID_color_map.keys(), key=lambda x: genotype_group_sort_hash[x]['Sort']), ordered=True)
        cdf = cdf.sort_values('Genotype')
        
        num_genotypes = len(cdf['Genotype'].unique())
        
        # Dynamic figure sizing with better margins for small plots
        height = max(num_genotypes * 0.5, 5)  # Increased minimum
        chr_length_mbp = cdf['ref_end'].max() / 1e6
        width = max(chr_length_mbp / 100 + 5, 12)
        width = min(width, 30)
        
        # Extra space for small sample counts (to fit title/legend)
        if num_genotypes <= 3:
            height += 1.5
        
        fig, ax = plt.subplots(figsize=(width, height))
        yticks = []
        yticklabels = []

        genotypes_to_plot = [g for g in cdf['Genotype'].unique()]

        if verbose: print(f"Plotting chromosome {chrom} with {num_genotypes} sample(s)...")

        for idx, genotype in enumerate(tqdm(genotypes_to_plot, desc=f"Plotting {chrom}", unit="genotype", disable=not verbose)):
            y_pos = len(genotypes_to_plot) - idx - 1
            yticks.append(y_pos)
            yticklabels.append(genotype)
                    
            genotype_rows = cdf[cdf['Genotype'] == genotype]
            for _, row in genotype_rows.iterrows():
                ax.add_patch(mpatches.Rectangle(
                    (row['ref_start_Mbp'], y_pos - 0.4),
                    row['ref_end_Mbp'] - row['ref_start_Mbp'],
                    0.8,
                    facecolor=row['color'],
                    edgecolor='none'
                ))

        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels, fontsize=10)
        ax.set_xlabel("Position (Mbp)", fontsize=12)
        ax.set_title(f"{chrom} Haplotype Blocks", fontsize=14, pad=20)
        
        # Add padding to y-axis limits
        ax.set_ylim(-0.5, len(genotypes_to_plot) - 0.5)

        # Dynamic X limits with margins
        if not cdf.empty:
            x_min = cdf['ref_start_Mbp'].min()
            x_max = cdf['ref_end_Mbp'].max()
            x_range = x_max - x_min
            if x_range < 10: 
                margin = max(x_range * 0.05, 0.5)
            elif x_range < 100: 
                margin = x_range * 0.02
            else: 
                margin = max(x_range * 0.01, 5.0)
            ax.set_xlim(x_min - margin, x_max + margin)

        # Legend
        legend_patches = []
        for genotype, color in hapID_color_map.items():
            group = genotype_group_sort_hash[genotype]['Group']
            if group in ['Reference', 'Pangenome']:
                legend_patches.append(mpatches.Patch(color=color, label=genotype))

        ax.legend(
            handles=legend_patches, 
            bbox_to_anchor=(1.02, 1), 
            loc='upper left', 
            borderaxespad=0.,
            fontsize=10
        )
        
        # Layout with warning suppression
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='.*Tight layout.*')
            try:
                plt.tight_layout(rect=[0.05, 0.05, 0.85, 0.95])
            except:
                plt.subplots_adjust(left=0.1, right=0.85, top=0.93, bottom=0.08)

        plt.savefig(chrom_output_png, dpi=300, bbox_inches='tight', pad_inches=0.3)
        plt.close()
        
        print(f"Plot saved: {chrom_output_png}")


#######################################################################################################
#                                          MAIN FUNCTION
#######################################################################################################

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        description='Generate haplotype painting plots from h.vcf files',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--hvcf-folder', required=True, help='Path to folder containing h.vcf files')
    parser.add_argument('--samples-list', required=True, help='Path to grouped samples file (TSV with Sort, Genotype, Group columns)')
    parser.add_argument('-c', '--chromosome', nargs='+', default=None, help='Chromosome(s) to plot (e.g., chr1H chr2H)')
    parser.add_argument('-r', '--region', type=str, default=None, help='Region to plot in format START-END (e.g. 1000-2000)')
    parser.add_argument('--plot-pangenome-references', action='store_true', default=False, help='Include pangenome samples in plots')
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help='Enable verbose output')

    parsed_args = parser.parse_args(args)

    # Extract args
    hvcfs_folder = parsed_args.hvcf_folder
    grouped_file = parsed_args.samples_list
    chr_to_plot = parsed_args.chromosome
    region = parsed_args.region
    plot_pangenomes = parsed_args.plot_pangenome_references
    verbose = parsed_args.verbose

    # --- Validation Logic Moved Here ---
    if not os.path.exists(grouped_file):
        print(f"\nERROR: Samples list file not found: {grouped_file}")
        sys.exit(1)

    region_start = None
    region_end = None   
    if region:
        try:
            region_parts = region.split('-')
            if len(region_parts) != 2: raise ValueError("Format error")
            region_start = int(region_parts[0])
            region_end = int(region_parts[1])
        except:
            print("Region must be START-END")
            sys.exit(1)

    # --- Execution Pipeline ---
    try:
        genotype_group_sort_hash = read_and_validate_grouped_file(grouped_file, verbose=verbose)
        hapID_color_map = add_colors_to_hapIDs(genotype_group_sort_hash)
        hvcf_files = find_hvcf_files(hvcfs_folder, genotype_group_sort_hash)
        
        # Call the internal tool logic (replacing os.system)
        bed_files = hvcf_to_bed(genotype_group_sort_hash, hvcf_files, verbose=verbose)
        
        chr_list = extract_chr_list(next(iter(bed_files.values())), genotype_group_sort_hash)
        
        df = generate_df(bed_files, hapID_color_map, genotype_group_sort_hash, hvcfs_folder, verbose=verbose)
        
        plot_haplotype_painting(
            hvcfs_folder, chr_to_plot, chr_list, hapID_color_map, 
            genotype_group_sort_hash, df, region_start, region_end, 
            plot_pangenomes, verbose=verbose
        )

    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
