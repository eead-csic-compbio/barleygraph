#!/usr/bin/env bash

# expects i) path to graph folder ($1)

# precomputed for each graph
# agc listset ${1}/vcf_dbs/assemblies.agc > /${1}/genomes.txt

while read g; do
	if [ ! -f "/gmap_db/${g}/${g}.chromosome" ]; then
		echo "formatting ${g}"
		agc getset -o "${1}/data/${g}.fa" "${1}/vcf_dbs/assemblies.agc" $g;
		/barleygraph/gmap-2013-08-31/local/bin/gmap_build -D /gmap_db/ -d "${g}" "${1}/data/${g}.fa" &> "/gmap_db/${g}.log";
		rm "${1}/data/${g}.fa";		
	fi
done < ${1}/genomes.txt
