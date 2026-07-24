## barleygraph

BARLEYGRAPH provides [PHG](https://github.com/maize-genetics/phg_v2)-based barley pangenome graphs 
for **sequence mapping** and **haplotype analysis**. This software is to be used from a 
[container](https://github.com/eead-csic-compbio?tab=packages&repo_name=barleygraph)
shipping with prebuilt PHG graphs of barley pangenomes and tools.

![PHG_database](https://github.com/eead-csic-compbio/barleygraph/blob/main/miscellaneous/Esquema_PHG.png)

>    It is based on 2.4 phg version.
>    
>    Find here the [2.4.75.230 release](https://github.com/maize-genetics/phg_v2/releases/tag/2.4.75.230).

### Introduction

Inspired by [BARLEYMAP](https://barleymap.eead.csic.es),
sequence alignments are performed with [GMAP](http://research-pub.gene.com/gmap), which supports both 
genomic sequences and transcripts. The series of genome sequences are scanned hierarchically; 
the scan stops with the first match. GMAP matches and precomputed graph ranges are intersected with
[BEDTOOLS](https://bedtools.readthedocs.io/en/latest). Genome assemblies compression and manage are done with [AGC 3.1](https://github.com/refresh-bio/agc).

> If your aiming is to align sequences and spot them in the several genomes of the graph, you may try the [Barleymap graph mode](https://barleymap.eead.csic.es/barleymap/graph/) on the website.
> It is user-friendly, quick and alignments are computed remotely. If you are looking for mapping a big dataset, add more genomes to the graph, or do haplotype analysis,
> you will need to work from the docker package, which requires terminal-commands knowledge and a bioinformatic server.

#### Available Graphs (Pan20 Release)

This release includes the **Pan20** graph. Note that the reference is **MorexV3** as annotated 
at [IPK](https://galaxy-web.ipk-gatersleben.de/libraries/folders/Fa676e8f07209a3be/dataset/78efbc10d9dd2218):

|graph|notes|genome names and scan order|
|:----|:----|:-----------|
|Pan20|Barley pangenome V1|MorexV3, Barke, HOR_9043, HOR_10350, HOR_3081, HOR_3365, Planet, HOR_7552, Akashinriki, OUN333, HOR_13942, HOR_13821, HOR_21599, Igri, Chiba, B1K-04-12, Du_Li_Huang, HOR_8148, GoldenPromise, Hockett|

---

### Quick Start Guide

**1. Pull the Docker Image**
Download the specific Pan20 release snapshot from the GitHub Container Registry:
```bash
docker pull ghcr.io/eead-csic-compbio/barleygraph:Pan20-20260721
```

Once downloaded, get into the docker image. Important to mark: You will need to specify a path for the gmap database. It will also be useful to move in and out files from your image to your local computer.
```bash
docker run --rm -v /path/to/local_gmap_db/:/gmap_db/ -it ghcr.io/eead-csic-compbio/barleygraph:Pan20-20260721 /bin/bash
```

> Since gmap index are bulky, those will be stored outside the container. More info below.

Several steps are required to run BARLEYGRAPH, depending if you want to do [mapping](https://github.com/eead-csic-compbio/barleygraph/edit/main/README.md#mapping) and/or [haplotype](https://github.com/eead-csic-compbio/barleygraph/edit/main/README.md#haplotype-analysis) analysis:

# Mapping

<img src="https://github.com/eead-csic-compbio/barleygraph/blob/main/miscellaneous/Align2graph_esquema.png"  width="400">

The figure below shows the full `align2graph` workflow:

<img src="https://github.com/eead-csic-compbio/barleygraph/blob/main/miscellaneous/align2graph_workflow.png" width="900">

### 1. Create local folder for GMAP indices [mapping]

This is done outside the container, as indices are bulky; in Linux you could do it as follows:

    mkdir -m 777 -p /path/to/local_gmap_db

### 2. Launch container and build GMAP indices [mapping]

The container will first be downloaded if not found locally.
The `index` command takes hours, over 8GB RAM and up to 163GB of disk, but it's only required the first time:

    docker run --rm -v /path/to/local_gmap_db/:/gmap_db/ -it ghcr.io/eead-csic-compbio/barleygraph:Pan20-20260721 index

| Genotype | Gmap index size (GB) |
| :--- | :--- |
| MorexV3 | 7.1 |
| Akashinriki | 7.1 |
| B1K-04-12 | 13 |
| Barke | 7.1 |
| Chiba | 7.2 |
| Du_Li_Huang | 7.1 |
| GoldenPromise | 7.1 |
| HOR_10350 | 7 |
| HOR_13821 | 7.1 |
| HOR_13942 | 7.1 |
| HOR_21599 | 7.1 |
| HOR_3081 | 7.1 |
| HOR_3365 | 7.1 |
| HOR_7552 | 7.1 |
| HOR_8148 | 7.1 |
| HOR_9043 | 7.1 |
| Hockett | 7.1 |
| Igri | 7.1 |
| OUN333 | 7.1 |
| Planet | 7.1 |

 
### 3. Map sequences in FASTA files [mapping]

Run this command line to find out about available optional flags in terminal:

    docker run --rm -v ~/local_gmap_db/:/gmap_db/ -it ghcr.io/eead-csic-compbio/barleygraph:Pan20-20260721 align2graph --help

Required arguments are ```--graph_yaml <mmap_pro/gmap_geno.yaml>``` and ```<input_fasta>```.

    docker run --rm -v ~/local_gmap_db/:/gmap_db/ -it ghcr.io/eead-csic-compbio/barleygraph:Pan20-20260721 align2graph --graph_yaml Pan20_mmap-pro.yaml sequences.fa
    docker run --rm -v ~/local_gmap_db/:/gmap_db/ -it ghcr.io/eead-csic-compbio/barleygraph:Pan20-20260721 align2graph --graph_yaml Pan20_gmap_geno.yaml sequences.fa

Optional arguments:

| Argument | Description |
| :--- | :--- |
| `--tmp_path TMP_PATH` | Path to writable folder for temporary files, default: /tmp |
| `--bedtools_exe BEDTOOLS_EXE` | Path to bedtools executable, default: bedtools |
| `--agc_exe AGC_EXE` | Path to agc executable, default: agc |
| `--minimap_exe MINIMAP_EXE` | Path to minimap executable, default: minimap2 |
| `--cor COR` | Number of cores for gmap, default: 4 |
| `--minident MINIDENT` | Min %identity of gmap matches, default: 98.0 |
| `--mincover MINCOVER` | Min %coverage of gmap matches, default: 95.0 |
| `--mincover_range MINCOVER_RANGE` | Min %coverage of gmap matches and pangenome ranges, default: 75.0 |
| `--single_genome SINGLE_GENOME` | Selected genome to be scanned with GMAP, must be part of graph, default: all genomes. Note that --add_ranges may not work properly with this option. |
| `--verb` | Increase verbosity in output |
| `--genomic` | Input sequences are genomic, turn off splicing |
| `--add_ranges {gmap,minimap,both}` | Add all pangenome ranges matching input sequences using specified tool (gmap, minimap, or both) |
| `--force_ranges` | When no graph overlap is found, search across all genomes using gmap to find ranges |


If ```--add_ranges <mode>```is on, output will contain each genome coordinates where your query sequence is found. You may use either GMAP, minimap2 or both at once, increasing required time but also accuracy.
You can read a better explanation at _Mapping pangenes_ section of [barleygraph paper](https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcRrTA6qsTDar992Z8SD7orh3o0ReV5Ng7c5lKsqFddk2g&s=10).

    docker run --rm -v ~/local_gmap_db/:/gmap_db/ -it ghcr.io/eead-csic-compbio/barleygraph:Pan20-20260721 align2graph --graph_yaml <config.yaml> --add_ranges <mode> sequences.fa


### 4. Mapping output

The parameters of a mapping run are included in the header as # comments.
The mapping results are in TSV format with the following columns:

|column name|explanations|
|:----------|:-----------|
|query|name of query sequence|
|ref_chr|name of chromosome in graph, taken from reference genome (MorexV3)|
|ref_start|1-based start coordinate in graph of range containing match, trimmed to GMAP coordinates if genome is reference MorexV3|
|ref_end|1-based end coordinate in graph of range containing match, trimmed to GMAP coordinates if genome is reference MorexV3|
|ref_strand|strand of graph range containing match, `.` if absent in reference genome|
|genome|name of genome of first GMAP match|
|chr|name of chromosome of first GMAP match|
|start|1-based start coordinate of first GMAP match|
|end|1-based end coordinate of first GMAP match|
|strand|strand of first GMAP match|
|perc_ident|% sequence identity of first GMAP match|
|perc_cover|% sequence cover of first GMAP match|
|multmaps|other GMAP matches (Yes/No)|
|graph_ranges|graph ranges of all genomes containing matching, requires flag `--add_ranges`|

Example output after mapping the VRN2 nucleotide sequence.

    # GMAP version: 2013-08-31
    # config_file: graph.yaml
    # fasta_file: old_bruno/VRN2.fa
    # minimum identity %: 98.0
    # minimum coverage %: 95.0
    # minimum coverage range %: 75.0
    # genomic: False
    # ranked pangenome genomes: MorexV3, HOR_2830, HOR_1168, HOR_14121, GDB_136, HOR_3365, HOR_3474, HOR_13942, HOR_21599, HOR_12184, HOR_2779, HOR_10892, HOR_21595

    #query	ref_chr		ref_start	ref_end		ref_strand	genome	chr	start	end	strand	perc_ident	perc_cover	multmaps	graph_ranges
    VRN2	chr4H_LR890099.1	604188191	604197211	.	HOR_2830	chr4H	602386783	602388450	-	98.8	100.0	No	.

# Haplotype analysis
kmer-based alignment of input FASTQ ﬁles from low-deep whole genome sequencing data. Infer haplotype blocks from the pangenome, achieving a pseudoassembly. It applies the [Viterbi algorithm](https://en.wikipedia.org/wiki/Viterbi_algorithm) to preserve local haplotype context, improving both ancestry inference and imputation accuracy

<img src="https://github.com/eead-csic-compbio/barleygraph/blob/main/miscellaneous/Imputation_esquema.png" width="300">

### Build kmer index

If you are looking for imputate sequencing data against our released databases, kmer index are precomputed; you don't need to run this step. 
If you have modified or increased a database, It is mandatory for afterwards imputate queries. 

> Requires high ammount of resources, >12 RAM GB, >100 GB of free -temporary- disk space, and up to 32 CPU hours. It depends on your machine and threads used (expect 1 hour if you provide 32 CPU cores).

    build_imputation_index
    build_imputation index -t <int>

 >    Note: Default threads if not specified = 8

It generates a [ropebwt3](https://github.com/lh3/ropebwt3) index file, with _.fmd_ extension.

Running the commnd ```build_imputation_index``` will allow you to choose which database are you interested on (mmap_pro or gmap_geno).
Index is moved to gmap_db external folder while remain accesible inside the docker image.

### Imputation

From fastq files, maps and get back the most likely haplotype paths. The resulting file is actually a [h.vcf](https://phg.maizegenetics.net/hvcf_specifications/) file.
Take as imput raw fastq reads (two files, if there are pair-ended reads). Files must end with ```.fastq```,```.fq```,```.fastq.gz``` or ```fq.gz```.

    imputation <R1_fastq> 
    imputation <R1_fastq> <R2_fastq>

>    Note: You can specify a number of threads to commit
> 
>    If empty field,  default value is 8 threads

    imputation <R1_fastq> -t <int>
    imputation <R1_fastq> <R2_fastq> -t <int>

You will be asked to choose which database you want to use (mmap_pro or gmap_geno). 

The output file is a **h**aplotype **v**ariant **c**all **f** or hVCF. Find more details at the [official specification documents](https://phg.maizegenetics.net/convenience_commands/#create-a-gff-file-from-an-imputed-hvcf-file). It is basically a gapless pseudoassembly based on the inference of haplotype blocks, where each line correspond to an individual block or range.


### Imputed haplotypes processing (optional)

Optional step: Convert h.vcf (haplotype VCF) files to BED format for easier analysis and visualization:

    hvcf2bed <vcf_folder> [genome_name]
    hvcf2bed <vcf_folder> [genome_name] -v

**Arguments:**
- `vcf_folder`: Path to folder containing h.vcf or h.vcf.gz files
- `genome_name` (optional): Specific genome to convert. If not provided, converts all h.vcf files in folder
- `-v, --verbose`: Enable verbose output for debugging

**Output:** Creates a BED file for each h.vcf file with columns: chrom, start, end, strand, checksum, genome, ref_chr, ref_start, ref_end, ref_checksum

**Example:**
```bash
hvcf2bed Med13/vcf_dbs/hvcf_files/ MorexV3
hvcf2bed Med13/vcf_dbs/hvcf_files/ -v  # Convert all files with verbose output
```

### Haplotype path painting

Generate visual plots of haplotype blocks from h.vcf files showing how different samples' genomes are composed of pangenome haplotypes:

> NOTE: You need to have all both pangenome and imputed .h.vcf files in the same folder. You may either move the imputed file to the folder with other pangenome files, the oposite (move the pangenome files into the output imputed files) or just create a new folder with soft links to group them up.

![Haplopainting_example](https://github.com/eead-csic-compbio/barleygraph/blob/main/miscellaneous/chr4H_FULL_haplotype_painting.png)

    haplopainting --hvcf-folder <path> --samples-list <file> [options]

**Arguments:**
- `--hvcf-folder` (required): Path to folder containing h.vcf files
- `--samples-list` (required): TSV file with columns: Sort, Genotype, Group (Reference/Pangenome/Imputed)
- `-c, --chromosome` (optional): Specific chromosome(s) to plot (e.g., chr1H chr2H). If not specified, all are plotted
- `-r, --region` (optional): Region to plot in format START-END (e.g., 1000000-2000000)
- `--plot-pangenome-references`: Include pangenome samples in plots (default: exclude)
- `-v, --verbose`: Enable verbose output showing progress

**Samples list format:**
```
Sort	Genotype	Group
1	MorexV3	Reference
2	HOR_2830	Pangenome
3	HOR_1168	Pangenome
...
15  YourSample  Imputed
...
```

>    Note: You can use `samplelist.tsv` as a template and add your imputed samples there.

**Output:** PNG plots saved in `<hvcf_folder>/plots/` directory

**Example:**
```bash
haplopainting --hvcf-folder Med13/vcf_dbs/hvcf_files/ --samples-list samplelist.tsv -v
haplopainting --hvcf-folder Med13/vcf_dbs/hvcf_files/ --samples-list samplelist.tsv -c chr1H chr2H --plot-pangenome-references
haplopainting --hvcf-folder Med13/vcf_dbs/hvcf_files/ --samples-list samplelist.tsv -c chr1H -r 1000000-2000000 -v
```

## Current status:

Pan20 at this moment is made up with mmap_pro or gmap_geno modes (more info in ongoing paper publication). You also have an Example dataset made of few base pairs of two arabidopsis chromosomes, useful to test some utilities and build new ones.

## Modifying or new databases

PHG databases are scalable and easy to modify. Adding new ```.h.vcf``` and ```.bed``` to a folder is enough to increase/decrease genomes to work on imputation and mapping. 
Docker image is fully equiped for it. To add new genomes, you will have to align sequences to the reference genome, and construct ```h.vcf``` files. Follow [PHG](https://github.com/maize-genetics/phg_v2) guidelines 
for an easy default (mmap_pro) alignment, or our in-terminal command indications in [barleygraph paper](https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcRrTA6qsTDar992Z8SD7orh3o0ReV5Ng7c5lKsqFddk2g&s=10) for other architectures.

If you want to build from scratch a new database, you may find interesting our [build PHG database script](graphs/build_PHG_database.sh). Just by modifing the [config file](graphs/PHG_database.example.config) and running the command like:
```bash
./build_PHG_database.sh --config <path/to/database.config>
```
All will be done automatically, for default _mmap\_pro_ architecture. For others, like _gmap\_geno_ you will needto substitute the ```phg align_assemblies``` command by two individual steps, on your convenience.

## For developers

Debuging the container may be though, its convenient to use a example toyset. Running:
```bash
docker pull ghcr.io/eead-csic-compbio/barleygraph:Example_Ara-20260721
```
You are getting a set of few arabidopsis genomes croped in some kb for chr1 & chr2 that are convinient to use.
To build an image using the [docker file](docker/Dockerfile), you only need a local file Example_Ara.tgz, which can not be uploaded here, but that you can export from the pulled docker. You might try:
```bash
docker run --rm   -v /scratch/PHG_barleymap/barleygraph/graphs/Ara_toyset/gmap_db/:/gmap_db/   -it barleygraph:example_ara   /bin/bash
tar -czvf Example_Ara .
```
Store the file in same folder as the Docker file and run it:
```bash
docker build \
  --build-arg graph=Example_Ara \
  -t <image_name:tag> \
  -f docker/Dockerfile .
```
To check how has this dataset being built: [Example_Ara config file](graphs/PHG_Example_Ara_database.config).

## Troubleshooting

If the `docker` commands above fail with an error similar to 

    permission denied while trying to connect to the Docker daemon socket

please check the instructions at https://docs.docker.com/engine/install/linux-postinstall

## References

See the files at [graphs/](https://github.com/eead-csic-compbio/barleygraph/tree/main/graphs) for the source of genome sequences and the MorexV3 gene annotation.

* Cantalapiedra CP, Boudiar R, Casas AM et al (2015) BARLEYMAP: physical and genetic mapping of nucleotide sequences and annotation of surrounding loci in barley. Mol Breeding 35:13. https://doi.org/10.1007/s11032-015-0253-1

* Bradbury PJ, Casstevens T, Jensen SE et al (2022) The Practical Haplotype Graph, a platform for storing and using pangenomes for imputation. Bioinformatics 38(15):3698-370. https://doi.org/10.1093/bioinformatics/btac410

* Wu TD, Watanabe CK (2005) GMAP: a genomic mapping and alignment program for mRNA and EST sequences. Bioinformatics 21(9):1859-1875. https://doi.org/10.1093/bioinformatics/bti310

* Quinlan AR, Hall IM (2010) BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26(6):841-842. https://doi.org/10.1093/bioinformatics/btq033

* Jayakodi M, Padmarasu S, Haberer G et al (2020) The barley pan-genome reveals the hidden legacy of mutation breeding. Nature 588:284-289. https://doi.org/10.1038/s41586-020-2947-8

* Mascher M, Wicker T, Jenkins J, et al (2021) Long-read sequence assembly: a technical evaluation in barley. The Plant Cell 33(6):1888-1906. https://doi.org/10.1093/plcell/koab077

## Citation


## Funding 

[PRIMA GENDIBAR, PCI2019-103526] supported by Horizon 2020 

[A08_23R] funded by Aragón goverment 

[FAS2022_052, INFRA24018, conexión BCB] funded by CSIC

[PID2022-142116OB-I00] by MICIU/AEI/10.13039/501100011033

![AEI](./miscellaneous/AEI.jpg)

