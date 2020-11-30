# GBS-Typer-sanger-nf
An updated NextFlow version of Ben Metcalf's GBS Typer pipeline.

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/GBS-Typer-sanger-nf/blob/main/LICENSE)   
![build](https://github.com/sanger-pathogens/GBS-Typer-sanger-nf/workflows/build/badge.svg)  
![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/sangerpathogens/gbs-typer-sanger-nf)   
[![codecov](https://codecov.io/gh/sanger-pathogens/GBS-Typer-sanger-nf/branch/main/graph/badge.svg)](https://codecov.io/gh/sanger-pathogens/GBS-Typer-sanger-nf)   

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
docker pull sangerpathogens/gbs-typer-sanger-nf:0.0.2
```

### Pipeline test (serotyping only)
```
cd GBS-Typer-sanger-nf
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'results'
```

To resume pipeline if incomplete:
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'results' -resume
```

### Run unit tests
```
python3 -m pytest
```
