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

### Pipeline test
```
cd GBS-Typer-sanger-nf
run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'Isolate_Sero_Res_Typing_results.txt'
```
With samples 26237_7#5, 26077_6#118 and 25292_2#85, should get isolate results:
```
III	0	0	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	1
II	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0
III	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0
```

To resume pipeline if incomplete:
```
run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'Isolate_Sero_Res_Typing_results.txt' -resume
```
