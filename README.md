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
genomic sequences and transcripts. The genome sequences making up a pangenome graph are scanned hierarchically; 
the scan stops with the first match. GMAP matches and precomputed graph ranges are intersected with
[BEDTOOLS](https://bedtools.readthedocs.io/en/latest). Genome assembly compression and management are done with 
[AGC 3.1](https://github.com/refresh-bio/agc).

> If your aim is to align barley sequences and locate them in the individual genomes of the graph, you may want to try the [graph mode](https://barleymap.eead.csic.es/barleymap/graph/) in the Barleymap Web app. It is user-friendly, quick and alignments are computed remotely. 


>If you need to align a large dataset or carry out haplotype analysis, you will need to work with the Docker image, which requires typing commands on the terminal and disk space. Check the guide below.

#### Available pangenome graphs

Currently this repository distributes flavours of the **Pan20** graph. 
Note that the reference is **MorexV3** as annotated 
at [IPK](https://galaxy-web.ipk-gatersleben.de/libraries/folders/Fa676e8f07209a3be/dataset/78efbc10d9dd2218), HC genes only:

|graph|notes|genome names and scan order|
|:----|:----|:-----------|
|Pan20|Barley pangenome V1|MorexV3, Barke, HOR_9043, HOR_10350, HOR_3081, HOR_3365, Planet, HOR_7552, Akashinriki, OUN333, HOR_13942, HOR_13821, HOR_21599, Igri, Chiba, B1K-04-12, Du_Li_Huang, HOR_8148, GoldenPromise, Hockett|

---

### Quick start Docker guide

**1. Pull the Docker Image**. Download an image from the GitHub Container Registry:

    docker pull ghcr.io/eead-csic-compbio/barleygraph:20260807

**2. Create local folders for GMAP indices, graphs and results**. This is done in the host computer, outside the container. This is required to keep the graph data separated from the code (Docker), and also to keep persistent copies of your results so that you can access them even when the Docker container is not running. You will need abundant disk space for the data. For instance, the downloadable `Pan20-mmap-pro` graph takes up 20GB and supports haplotype analysis only. You would need another 150GB should you build the GMAP indices required for align2grap. For instance, in Linux you could create the following folders in your home:

    mkdir -m 777 -p ${HOME}/graph_db          #required
    mkdir -m 777 -p ${HOME}/graph_db/input    #to place input FASTQ/FASTA files
    mkdir -m 777 -p ${HOME}/graph_db/results  #to store output files
    mkdir -m 777 -p ${HOME}/gmap_db           #only if you plant to run align2graph

**3. Run the image binding the local folders**. Bindings look like /full/path/local:/container. The following command will launch a container on your terminal, the promot should be similar to you@8ee9e89ed09c:/barleygraph$ :

    docker run -it -v ${HOME}/gmap_db/:/gmap_db -v ${HOME}/graph_db:/graph_db barleygraph:2026-08-07

**4. Setup a graph**. At the container terminal type and run:

    setup_graph -l                     #to see currently supported graphs
    setup_graph -G Pan20-mmap-pro
    setup_graph -G Pan20-mmap-pro -g   #to additionally make GMAP indices; can take hours

**5. Imputate and call haplotypes** requires 1 single-end or 2 pair-end FASTQ files, which might be compressed:

    imputation -G Pan20-mmap-pro -1 /graph_db/input/HOR_10096_GBS.fq -o /graph_db/results/

This command should a hVCF output file `/graph_db/results/HOR_10096_GBS_1.h.vcf` and a folder `/graph_db/results/HOR_10096_GBS_1.hvcfdir/` that we can use in the next setp.

**6. Paint haplotypes**:

    haplopainting -h  #checkout options

    haplopainting --hvcf-folder /graph_db/results//HOR_10096_GBS_1.hvcfdir/ --samples-list /graph_db/results//HOR_10096_GBS_1.hvcfdir/Pan20_samplelist.tsv -f pdf --plot-pangenome-references

In addition to BED files converted from the original vVCF, this command will produce haplotype plots in folder `/graph_db/results/HOR_10096_GBS_1.hvcfdir/plots/`, one per chromosome, with graph genomes on top and sample below:

<img src="https://github.com/eead-csic-compbio/barleygraph/blob/main/miscellaneous/chr4H_FULL_haplotype_painting.png"  width="400">

**7. Mapping sequences in FASTA files**.

    align2graph -h

    align2graph --graph_yaml /graph_db/Pan20/Pan20-mmap-pro/Pan20-mmap-pro.yaml /graph_db/input/Vrn2.fna 

### More details

The next figures describe the `align2graph` algorithm:

<img src="https://github.com/eead-csic-compbio/barleygraph/blob/main/miscellaneous/Align2graph_esquema.png"  width="400">

<img src="https://github.com/eead-csic-compbio/barleygraph/blob/main/miscellaneous/align2graph_workflow.png" width="900">

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

#### Imputation and haplotype analysis

<img src="https://github.com/eead-csic-compbio/barleygraph/blob/main/miscellaneous/Imputation_esquema.png" width="300">

The resulting file is a [h.vcf](https://phg.maizegenetics.net/hvcf_specifications/) file, a **h**aplotype **v**ariant **c**all **f** or hVCF. Find more details at the [official specification documents](https://phg.maizegenetics.net/convenience_commands/#create-a-gff-file-from-an-imputed-hvcf-file). It is essentially a gapless pseudoassembly based on the inference of haplotype blocks, where each line corresponds to an individual block or range.

Generate visual plots of haplotype blocks from h.vcf files showing how different samples' genomes are composed of pangenome haplotypes.


### Troubleshooting

If the `docker` commands above fail with an error similar to 

    permission denied while trying to connect to the Docker daemon socket

please check the instructions at https://docs.docker.com/engine/install/linux-postinstall

### References

See the files at [graphs/](https://github.com/eead-csic-compbio/barleygraph/tree/main/graphs) for the source of genome sequences and the MorexV3 gene annotation.

* Cantalapiedra CP, Boudiar R, Casas AM et al (2015) BARLEYMAP: physical and genetic mapping of nucleotide sequences and annotation of surrounding loci in barley. Mol Breeding 35:13. https://doi.org/10.1007/s11032-015-0253-1

* Bradbury PJ, Casstevens T, Jensen SE et al (2022) The Practical Haplotype Graph, a platform for storing and using pangenomes for imputation. Bioinformatics 38(15):3698-370. https://doi.org/10.1093/bioinformatics/btac410

* Wu TD, Watanabe CK (2005) GMAP: a genomic mapping and alignment program for mRNA and EST sequences. Bioinformatics 21(9):1859-1875. https://doi.org/10.1093/bioinformatics/bti310

* Quinlan AR, Hall IM (2010) BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26(6):841-842. https://doi.org/10.1093/bioinformatics/btq033

* Jayakodi M, Padmarasu S, Haberer G et al (2020) The barley pan-genome reveals the hidden legacy of mutation breeding. Nature 588:284-289. https://doi.org/10.1038/s41586-020-2947-8

* Mascher M, Wicker T, Jenkins J, et al (2021) Long-read sequence assembly: a technical evaluation in barley. The Plant Cell 33(6):1888-1906. https://doi.org/10.1093/plcell/koab077

### Citation


##3 Funding 

This work was supported by AEI/10.13039/501100011033/FEDER/UE [PID2022-142116OB-I00 and predoctoral contract PREP2022_EEAD_52 to JSA], Horizon 2020 PRIMA [PCI2019-103526] and SusCrop ERA-NET Recobar [771134], Government of Aragon [A08_23R] and CSIC [FAS2022_052, INFRA24018].

![AEI](./miscellaneous/AEI.jpg)
