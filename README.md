## barleygraph

BARLEYGRAPH provides [PHG](https://github.com/maize-genetics/phg_v2)-based barley pangenome graphs 
for **sequence mapping** and **haplotype analysis**. This software is to be used from a 
[container](https://github.com/eead-csic-compbio?tab=packages&repo_name=barleygraph)
shipping with prebuilt graphs and tools.
For convenience, sequence mapping can also be done on the Barleymap Web app, [graph mode](https://barleymap.eead.csic.es/barleymap/graph/).

![PHG_database](https://github.com/eead-csic-compbio/barleygraph/blob/main/miscellaneous/Esquema_PHG.png)

>    It is based on PHG v2.4.
>    
>    Find here the [2.4.75.230 release](https://github.com/maize-genetics/phg_v2/releases/tag/2.4.75.230).

### Introduction

Inspired by [BARLEYMAP](https://barleymap.eead.csic.es),
sequence alignments are performed with [GMAP](http://research-pub.gene.com/gmap), which supports both 
genomic sequences and transcripts. The genome sequences making up a pangenome graph are scanned hierarchically starting with MorexV3; 
the scan stops with the first match. GMAP matches and precomputed graph ranges are intersected with
[BEDTOOLS](https://bedtools.readthedocs.io/en/latest). Genome assembly compression and management are done with 
[AGC](https://github.com/refresh-bio/agc).

> If your aim is to align barley sequences and locate them in the individual genomes of the graph, you may want to try the [graph mode](https://barleymap.eead.csic.es/barleymap/graph/) in the Barleymap Web app. It is user-friendly, quick and alignments are computed remotely. 


>If you need to align a large dataset or carry out haplotype analysis, you will need to work with the Docker image, which requires typing commands on the terminal and disk space. Check the guide below.

#### Available pangenome graphs

Currently this repository distributes flavours of the **Pan20** graph. 
Note that the reference is **MorexV3** as annotated at 
[IPK](https://galaxy-web.ipk-gatersleben.de/libraries/folders/Fa676e8f07209a3be/dataset/78efbc10d9dd2218), HC genes only:

|graph|genomes sorted by scan order/contributed ranges|
|:----|:----------------------------------------------|
|Pan20-mmap-pro|MorexV3, Hockett, Igri, Du_Li_Huang, Planet, HOR_9043, Barke, HOR_13821, HOR_21599, HOR_3365, HOR_8148, HOR_13942, GoldenPromise, B1K-04-12, HOR_10350, HOR_7552, OUN333, Chiba, Akashinriki, HOR_3081|
|Pan20-gmap-geno|MorexV3, HOR_21599, GoldenPromise, Du_Li_Huang, Barke, HOR_7552, HOR_13821, Hockett, OUN333, Igri, HOR_8148, HOR_3081, Planet, HOR_9043, HOR_13942, Akashinriki, HOR_3365, Chiba, HOR_10350, B1K-04-12|

### Quick start Docker guide

**1. Pull the Docker Image**. Download an image from the GitHub Container Registry:

    docker pull ghcr.io/eead-csic-compbio/barleygraph:2026-09-04

**2. Create local persistent folders for graphs & GMAP indices**. This is done in the host computer, outside the container. This is required to keep persistent graph data separated from the code (Docker). You will need abundant disk space for the data. For instance, the downloadable `Pan20-mmap-pro` graph takes up to 20GB and supports haplotype analysis only. You would need another 150GB should you build the GMAP indices required for align2grap. For instance, in Linux you could create the following folders in your home:

    mkdir -p ${HOME}/graph_db     #required
    mkdir -p ${HOME}/gmap_db      #required to run align2graph

These folders will be bound by the Docker container at runtime. Binding arguments look like this: `/full/path/local:/container`.

**3. Create a local folder for results (optional)**. If you run the `imputation` and `haplopainting` scripts you will need also a writable folder to store results which you can review even when the container is off:

    mkdir -p ${HOME}/results

**4. Binding the local input folder**. You will need to add another argument to bind the folder containing your input FASTA and FASTQ files. For instance this could be a folder in your home which will be referrod to within the container as `/user`:

    -v ${HOME}/datafiles/:/user

See examples below on how to actually analyze your own input.

**5. Check and setup graphs**. Try the following commands, the argument `-u $(id -u):$(id -g)` is required to make sure any written files belong to the user running the container:

    # list graphs that can be downloaded from this site
    docker run -it -v ${HOME}/graph_db:/graph_db ghcr.io/eead-csic-compbio/barleygraph:latest setup_graph -l

    # see currently locally installed graphs
    docker run -it -v ${HOME}/graph_db:/graph_db ghcr.io/eead-csic-compbio/barleygraph:latest setup_graph -I

    # download and install a listed graph, this will take an hour
    docker run -it -u $(id -u):$(id -g) -v ${HOME}/graph_db:/graph_db ghcr.io/eead-csic-compbio/barleygraph:latest setup_graph -G Pan20-mmap-pro

    # optionally make GMAP indices; required to run align2graph, will take more time 
    docker run -it -u $(id -u):$(id -g) -v ${HOME}/gmap_db/:/gmap_db -v ${HOME}/graph_db:/graph_db ghcr.io/eead-csic-compbio/barleygraph:latest setup_graph -G Pan20-mmap-pro -g

**6. Imputate and call haplotypes** requires 1 single-end or 2 pair-end FASTQ files, which might be compressed. This step requires over 24GB RAM:

    # check options
    docker run -it ghcr.io/eead-csic-compbio/barleygraph:latest imputation

    # run test FASTQ file, binding args are stored in variable for convenience
    BINDS="-v ${HOME}/results:/results -v ${HOME}/graph_db:/graph_db"
    docker run -it -u $(id -u):$(id -g) ${BINDS} ghcr.io/eead-csic-compbio/barleygraph:latest imputation -G Pan20-mmap-pro -1 test.fq -o /results/

    # example with input data provided by user, see step 4 and BINDS below
    BINDS="-v ${HOME}/datafiles/:/user -v ${HOME}/results:/results -v ${HOME}/graph_db:/graph_db"
    docker run -it -u $(id -u):$(id -g) ${BINDS} ghcr.io/eead-csic-compbio/barleygraph:latest imputation -G Pan20-mmap-pro -1 /user/example.fq -o /results/

Thess command produce a hVCF output file `..._1.h.vcf` and a folder `..._1.hvcfdir/` that we can use in the next step.
**Note**: If you want to imputate with several graphs make sure you use different output folders to store the results.

**7. Paint haplotypes**. This requires results produced in the previous step, which we add to shell variable PREVRES for convenience:

    # check options
    docker run ghcr.io/eead-csic-compbio/barleygraph:latest haplopainting -h

    # run with previous imputation results
    BINDS="-v ${HOME}/results:/results -v ${HOME}/graph_db:/graph_db"
    PREVRES="--hvcf-folder /results/..._1.hvcfdir/ --samples-list /results/..._1.hvcfdir/Pan20_samplelist.tsv"
    docker run -it -u $(id -u):$(id -g) ${BINDS} ghcr.io/eead-csic-compbio/barleygraph:latest haplopainting ${PREVRES} -f pdf --plot-pangenome-references

In addition to BED files converted from the original hVCF, this command will produce haplotype plots in folder `/results/..._1.hvcfdir/plots/`, one per chromosome, with graph genomes on top and sample below:

<img src="https://github.com/eead-csic-compbio/barleygraph/blob/main/miscellaneous/chr4H_FULL_haplotype_painting.png"  width="400">

**8. Mapping sequences in FASTA files**. This requires GMAP indices, see step 5:

    # check options 
    docker run ghcr.io/eead-csic-compbio/barleygraph:latest align2graph -h

    # run test FASTA file, YAML config file for relevant graph is required
    BINDS="-v ${HOME}/gmap_db:/gmap_db -v ${HOME}/graph_db:/graph_db"
    YAML="--graph_yaml /graph_db/Pan20/Pan20-mmap-pro/Pan20-mmap-pro.yaml"
    docker run -it ${BINDS} ghcr.io/eead-csic-compbio/barleygraph:latest align2graph ${YAML} test.fa
    
    # request aligned segments in all matched graph genomes
    docker run -it ${BINDS} ghcr.io/eead-csic-compbio/barleygraph:latest align2graph ${YAML} test.fa --add_ranges both

    # example with input data provided by user
    BINDS="-v ${HOME}/datafiles/:/user -v ${HOME}/results:/results -v ${HOME}/graph_db:/graph_db"
    docker run -it ${BINDS} ghcr.io/eead-csic-compbio/barleygraph:latest align2graph ${YAML} /user/example.fa --add_ranges both

**9. Get citation**.

    docker run -it ghcr.io/eead-csic-compbio/barleygraph:latest citation


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

    # version: 8bbfb5f 2026-05-28
    # GMAP version: 2013-08-31
    # config_file: /graph_db/Pan20/Pan20-mmap-pro/Pan20-mmap-pro.yaml
    # fasta_file: /graph_db/input/Vrn2.fna
    # minimum identity %: 98.0
    # minimum coverage %: 95.0
    # minimum coverage range %: 75.0
    # add_ranges: False (Tool: None)
    # force_ranges: False
    # genomic: False
    # ranked pangenome genomes: MorexV3, Hockett, Igri, Du_Li_Huang, Planet, HOR_9043, Barke, HOR_13821, HOR_21599, HOR_3365, HOR_8148, HOR_13942, GoldenPromise, B1K-04-12, HOR_10350, HOR_7552, OUN333, Chiba, Akashinriki, HOR_3081

    #query  ref_chr         ref_start       ref_end         ref_strand      genome  chr     start   end     strand  perc_ident       perc_cover      multmaps        graph_ranges
    Horvu_13942_4H01G516500.1      chr4H_LR890099.1        604188191       604202141       .       Igri     chr4H   602527414       602529082       -       98.0    100.0   No      .

Note that pangenome genomes are ranked by contributed ranges; those number will change across graphs.

#### Imputation and haplotype analysis

<img src="https://github.com/eead-csic-compbio/barleygraph/blob/main/miscellaneous/Imputation_esquema.png" width="300">

The resulting file is a [h.vcf](https://phg.maizegenetics.net/hvcf_specifications/) file, a **h**aplotype **v**ariant **c**all **f** or hVCF. Find more details at the [official specification documents](https://phg.maizegenetics.net/convenience_commands/#create-a-gff-file-from-an-imputed-hvcf-file). It is essentially a gapless pseudoassembly based on the inference of haplotype blocks, where each line corresponds to an individual block or range.

Generate visual plots of haplotype blocks from h.vcf files showing how different samples' genomes are composed of pangenome haplotypes.


### Troubleshooting

<!-- -u $(id -u):$(id -g) -->

* If docker fails to download container due to no disk space left you might need to [change the location](https://evodify.com/change-docker-storage-location) of images, will require admin.

* If the `docker` commands above fail with an error similar to `permission denied while trying to connect to the Docker daemon socket` please check the instructions at https://docs.docker.com/engine/install/linux-postinstall

* On Windows [WSL](https://learn.microsoft.com/es-es/windows/wsl/install) the container often requires [sudo](https://stackoverflow.com/questions/64710480/docker-client-under-wsl2-doesnt-work-without-sudo) to run.

* Downloading a graph (step 5) might fail due to network issues. In this case a solution is to run the container interactively with `docker run -it -v ${HOME}/graph_db:/graph_db ghcr.io/eead-csic-compbio/barleygraph:latest` and then run the setup command on the prompt, for instance `setup_graph -G Pan20-mmap-pro`. You might have to repeat the setup command several times if the connection is poor, but this way all the successfully downloaded graph parts are skipped.

* Building GMAP indices (step 5) requires some RAM; if you run out of memory you'll get an error message like `Killed agc getset`.

* Imputation with Pan20 graphs requires over 24GB RAM; if you run out of memory you'll get an error message like `Exception in thread "main" java.lang.OutOfMemoryError: Java heap space`.


### References

See the files at [graphs/](https://github.com/eead-csic-compbio/barleygraph/tree/main/graphs) for the source of genome sequences and the MorexV3 gene annotation.

* Cantalapiedra CP, Boudiar R, Casas AM et al (2015) BARLEYMAP: physical and genetic mapping of nucleotide sequences and annotation of surrounding loci in barley. Mol Breeding 35:13. https://doi.org/10.1007/s11032-015-0253-1

* Bradbury PJ, Casstevens T, Jensen SE et al (2022) The Practical Haplotype Graph, a platform for storing and using pangenomes for imputation. Bioinformatics 38(15):3698-370. https://doi.org/10.1093/bioinformatics/btac410

* Wu TD, Watanabe CK (2005) GMAP: a genomic mapping and alignment program for mRNA and EST sequences. Bioinformatics 21(9):1859-1875. https://doi.org/10.1093/bioinformatics/bti310

* Quinlan AR, Hall IM (2010) BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26(6):841-842. https://doi.org/10.1093/bioinformatics/btq033

* Jayakodi M, Padmarasu S, Haberer G et al (2020) The barley pan-genome reveals the hidden legacy of mutation breeding. Nature 588:284-289. https://doi.org/10.1038/s41586-020-2947-8

* Mascher M, Wicker T, Jenkins J, et al (2021) Long-read sequence assembly: a technical evaluation in barley. The Plant Cell 33(6):1888-1906. https://doi.org/10.1093/plcell/koab077

### Citation

Sarria J, Amhal H, Ramírez CJ, Igartua E, Casas AM, Contreras-Moreira B (2026) A pangenome-graph approach for mapping and imputing barley sequences. bioRxiv 2026.08.06.741139; doi: https://doi.org/10.64898/2026.08.06.741139

## Funding 

This work was supported by AEI/10.13039/501100011033/FEDER/UE [PID2022-142116OB-I00 and predoctoral contract PREP2022_EEAD_52 to JSA], Horizon 2020 PRIMA [PCI2019-103526] and SusCrop ERA-NET Recobar [771134], Government of Aragon [A08_23R] and CSIC [FAS2022_052, INFRA24018].

![AEI](./miscellaneous/AEI.jpg)
