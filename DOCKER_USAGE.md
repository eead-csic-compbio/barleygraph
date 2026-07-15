# Docker Image Build and Usage Guide

## Building the Docker Image

### Prerequisites
- Ensure you have the Pan20 database files ready:
  - `Pan20.tgz` (compressed pangenome database)
  - `index_Pan20` (index file)
  
- These files should be in the repository root directory before building.

### Build Commands

**For Pan20:**
```bash
cd /scratch/PHG_barleymap/barleygraph
docker build \
  --build-arg graph=Pan20 \
  -t ghcr.io/eead-csic-compbio/barleygraph:pan20-$(date +%Y%m%d) \
  -f docker/Dockerfile .
```

**For Med13:**
```bash
cd /scratch/PHG_barleymap/barleygraph
docker build \
  --build-arg graph=Med13 \
  -t ghcr.io/eead-csic-compbio/barleygraph:med13-$(date +%Y%m%d) \
  -f docker/Dockerfile .
```

## Running the Docker Container

### Interactive Shell
```bash
docker run --rm \
  -v /path/to/local_gmap_db/:/gmap_db/ \
  -it ghcr.io/eead-csic-compbio/barleygraph:pan20-20251008 \
  /bin/bash
```

### Building Imputation Index
```bash
docker run --rm \
  -v /path/to/local_gmap_db/:/gmap_db/ \
  -v /path/to/working/dir:/workdir \
  -w /workdir \
  ghcr.io/eead-csic-compbio/barleygraph:pan20-20251008 \
  build_imputation_index -t 8
```

### Running Imputation
```bash
docker run --rm \
  -v /path/to/local_gmap_db/:/gmap_db/ \
  -v /path/to/working/dir:/workdir \
  -w /workdir \
  ghcr.io/eead-csic-compbio/barleygraph:pan20-20251008 \
  imputation reads/sample_r1.fq.gz reads/sample_r2.fq.gz -t 8
```

### Using the Index Script
```bash
docker run --rm \
  -v /path/to/local_gmap_db/:/gmap_db/ \
  -v /path/to/working/dir:/workdir \
  -w /workdir \
  ghcr.io/eead-csic-compbio/barleygraph:pan20-20251008 \
  index
```

## Environment Variables

- `CONDA_ENV_PREFIX`: Override the conda environment path (default: `/opt/conda/envs/phgv2.4/`)
- `_JAVA_OPTIONS`: Set to `-Xmx256g` by default for large memory allocation

## Volume Mounts

- `/gmap_db/`: Writable directory for GMAP database files (bound at runtime)
- `/workdir`: Your working directory with input reads and output location

## Notes

- The image includes:
  - PHG v2.4.74.229
  - GMAP v2013-08-31
  - Python conda environment (phgv2.4) with pandas, matplotlib, tqdm, pyyaml
  - align2graph script
  - Both Pan20 variants (gmap-geno, mmap-pro) if Pan20 was selected at build time

- When running `build_imputation_index` or `imputation` with Pan20, you will be prompted to choose between:
  1. gmap-geno
  2. mmap-pro
