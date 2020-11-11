# GBS-Typer-sanger-nf
An updated NextFlow version of Ben Metcalf's GBS Typer pipeline.

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/GBS-Typer-sanger-nf/blob/master/LICENSE)   
![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/sangerpathogens/gbs-typer-sanger-nf)   

### Installation locally
GBS Typer relies on Nextflow and Docker.
Download:
1. [Docker](https://www.docker.com/).
2. [Nextflow](https://www.nextflow.io/).
3. Clone repository:
```
git clone https://github.com/sanger-pathogens/GBS-Typer-sanger-nf.git
```
4. Pull Docker image
```
docker pull sangerpathogens/gbs-typer-sanger-nf:latest
```

### Pipeline test
```
cd GBS-Typer-sanger-nf
nextflow run main.nf --reads 'data/example_{1,2}.fastq.gz' --contigs 'data/example_contigs.fa'
```
Two fastq files should be found a new results directory.
