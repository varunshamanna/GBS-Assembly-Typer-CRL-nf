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
docker pull sangerpathogens/gbs-typer-sanger-nf:0.0.1
```

### Pipeline test (serotyping only)
```
cd GBS-Typer-sanger-nf
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'Isolate_Sero_Res_Typing_results.txt'
```

To resume pipeline if incomplete:
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'Isolate_Sero_Res_Typing_results.txt' -resume
```

The reference resistance database for the resistance typing can be specified with --custom_res, --argannot and/or --resfinder, e.g. to run with ARGANNOT only:
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'Isolate_Sero_Res_Typing_results.txt' --argannot
```

To use all databases specify --all, e.g.
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'Isolate_Sero_Res_Typing_results.txt' --all
```

### Run unit tests
```
nosetests
```
