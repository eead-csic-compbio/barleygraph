## barleygraph

[PHG](https://github.com/maize-genetics/phg_v2)-based barley pangenome graphs for sequence mapping and haplotype analysis.

This software is to be used from a container shipping with prebuilt PHG graphs of barley pangenomes and tools. Two graphs are available:

|graph|notes|genome names|
|:----|:----|:-----------|
|Med13|Merges 13 genomes of Mediterranean barley landraces|MorexV3 HOR_10892 HOR_1168 HOR_12184 HOR_13942 HOR_14121 HOR_21595 HOR_21599 HOR_2779 HOR_2830 HOR_3365 HOR_3474 GDB_136|
|Pan20|Barley pangenome V1|MorexV3 Akashinriki B1K-04-12 Barke Chiba Du_Li_Huang GoldenPromise HOR_10350 HOR_13821 HOR_13942 HOR_21599 HOR_3081 HOR_3365 HOR_7552 HOR_8148 HOR_9043 Hockett Igri OUN333 Planet|

Several steps are required to run it, depending if you want to do [mapping] and/or [haplotype} analysis:

1. Create local folder for GMAP indices, outside the container, as these are bulky. In Linux you can do:

    mkdir -m 777 -p /path/to/local_gmap_db

2. Launch container and build indices, takes hours and up to 223GB of disk [mapping]

    docker run --rm -v /path/to/local_gmap_db/:/gmap_db/ -it eeadcsiccompbio/barleygraph:Med13-latest index

    # alternatively
    docker run --rm -v /path/to/local_gmap_db/:/gmap_db/ -it eeadcsiccompbio/barleygraph:Pan20-latest index
    
3. Map sequences inf FASTA files [mapping]

    docker run --rm -v ~/local_gmap_db/:/gmap_db/ -it eeadcsiccompbio/barleygraph:Med13-latest align2graph sequences.fa

