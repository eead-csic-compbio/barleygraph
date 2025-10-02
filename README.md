## barleygraph

BARLEYGRAPH provides [PHG](https://github.com/maize-genetics/phg_v2)-based barley pangenome graphs 
for sequence mapping and haplotype analysis. This software is to be used from a Docker container 
shipping with prebuilt PHG graphs of barley pangenomes and tools. 

Inspired by [BARLEYMAP](https://barleymap.eead.csic.es),
sequence alignments are performed with [GMAP](http://research-pub.gene.com/gmap), which supports both 
genomic sequences and transcripts. GMAP matches and precomputed graph ranges are intersected with
[BEDTOOLS](https://bedtools.readthedocs.io/en/latest).

Two graphs are available, which include a series of genome sequences which are scanned hierarchically;
the scan stops with the first match. Note that in both cases the reference is **MorexV3** as annotated 
at IPK:

|graph|notes|genome names and scan order|
|:----|:----|:-----------|
|Med13|Merges 13 genomes of Mediterranean barley landraces|MorexV3, HOR_2830, HOR_1168, HOR_14121, GDB_136, HOR_3365, HOR_3474, HOR_13942, HOR_21599, HOR_12184, HOR_2779, HOR_10892, HOR_21595|
|Pan20|Barley pangenome V1|MorexV3, Barke, HOR_9043, HOR_10350, HOR_3081, HOR_3365, Planet, HOR_7552, Akashinriki, OUN333, HOR_13942, HOR_13821, HOR_21599, Igri, Chiba, B1K-04-12, Du_Li_Huang, HOR_8148, GoldenPromise, Hockett|

## howto to run

Several steps are required to run it, depending if you want to do [mapping] and/or [haplotype} analysis:

1. Create local folder for GMAP indices, outside the container, as these are bulky. In Linux you can do:

    mkdir -m 777 -p /path/to/local_gmap_db

2. Launch container and build GMAP indices, takes hours and up to 223GB of disk **[mapping]**

    docker run --rm -v /path/to/local_gmap_db/:/gmap_db/ -it eeadcsiccompbio/barleygraph:Med13-latest index

    or

    docker run --rm -v /path/to/local_gmap_db/:/gmap_db/ -it eeadcsiccompbio/barleygraph:Pan20-latest index
    
3. Map sequences in FASTA files **[mapping]**

    docker run --rm -v ~/local_gmap_db/:/gmap_db/ -it eeadcsiccompbio/barleygraph:Med13-latest align2graph # see optional params

    docker run --rm -v ~/local_gmap_db/:/gmap_db/ -it eeadcsiccompbio/barleygraph:Med13-latest align2graph sequences.fa

    docker run --rm -v ~/local_gmap_db/:/gmap_db/ -it eeadcsiccompbio/barleygraph:Med13-latest align2graph --add_ranges sequences.fa


## mapping output

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



## References

See the files at [graphs/](https://github.com/eead-csic-compbio/barleygraph/tree/main/graphs) for the source of genome sequences and the MorexV3 gene annotation.

* Cantalapiedra CP, Boudiar R, Casas AM et al (2015) BARLEYMAP: physical and genetic mapping of nucleotide sequences and annotation of surrounding loci in barley. Mol Breeding 35:13. https://doi.org/10.1007/s11032-015-0253-1

* Bradbury PJ, Casstevens T, Jensen SE et al (2022) The Practical Haplotype Graph, a platform for storing and using pangenomes for imputation. Bioinformatics 38(15):3698-370. https://doi.org/10.1093/bioinformatics/btac410

* Wu TD, Watanabe CK (2005) GMAP: a genomic mapping and alignment program for mRNA and EST sequences. Bioinformatics 21(9):1859-1875. https://doi.org/10.1093/bioinformatics/bti310

* Quinlan AR, Hall IM (2010) BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26(6):841-842. https://doi.org/10.1093/bioinformatics/btq033

* Jayakodi M, Padmarasu S, Haberer G et al (2020) The barley pan-genome reveals the hidden legacy of mutation breeding. Nature 588:284-289. https://doi.org/10.1038/s41586-020-2947-8

* Mascher M, Wicker T, Jenkins J, et al (2021) Long-read sequence assembly: a technical evaluation in barley. The Plant Cell 33(6):1888-1906. https://doi.org/10.1093/plcell/koab077


## Funding 

[PRIMA GENDIBAR, PCI2019-103526] supported by Horizon 2020 

[A08_23R] funded by Aragón goverment 

[FAS2022_052, INFRA24018, conexión BCB] funded by CSIC

[PID2022-142116OB-I00] by MICIU/AEI/10.13039/501100011033

![AEI](./AEI.jpg)

