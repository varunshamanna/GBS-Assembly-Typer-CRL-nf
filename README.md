# GBS-Typer-sanger-nf
An updated NextFlow version of [Ben Metcalf's GBS Typer pipeline](https://github.com/BenJamesMetcalf/GBS_Scripts_Reference).

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/GBS-Typer-sanger-nf/blob/main/LICENSE)   
![build](https://github.com/sanger-pathogens/GBS-Typer-sanger-nf/workflows/build/badge.svg)  
![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/sangerpathogens/gbs-typer-sanger-nf)   
[![codecov](https://codecov.io/gh/sanger-pathogens/GBS-Typer-sanger-nf/branch/main/graph/badge.svg)](https://codecov.io/gh/sanger-pathogens/GBS-Typer-sanger-nf)   

### Installation on local machine
GBS Typer relies on Nextflow and Docker.
Download:
1. [Docker](https://www.docker.com/).
2. [Nextflow](https://www.nextflow.io/).
3. Clone repository:
```
git clone https://github.com/sanger-pathogens/GBS-Typer-sanger-nf.git
```

### Running the pipeline
#### To run on one sample:
```
cd GBS-Typer-sanger-nf
nextflow run main.nf --reads 'data/sampleID_{1,2}.fastq.gz' --output 'sampleID'
```
This will create three tab-delimited files in a 'results' directory within the current directory:
1. **sampleID_serotype_res_incidence.txt** - Gives the serotype and presence/absence (i.e. +/-) of antibiotic resistance genes (GBS-specific alleles and ResFinder/ARG-ANNOT genes)
e.g. Isolate Strep B sample 25292_2#105 has serotype II and have genes: 23S1, 23S3, GYRA, LSAC and TETM

ID | Serotype | 23S1 | 23S3 | CAT | ERMB | ERMT | FOSA | GYRA | LNUB | LSAC | MEFA | MPHC | MSRA | MSRD | PARC | RPOBGBS-1 | RPOBGBS-2 | RPOBGBS-3 | RPOBGBS-4 | SUL2 | TETB | TETL | TETM | TETO | TETS
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
25292_2#105 | II | + | + | - | - | - | - | + | - | + | - | - | - | - | - | - | - | - | - | - | - | - | + | - | -

2. **sampleID_gbs_res_variants.txt** - Gives the SNP variants for GBS-specific resistance genes
e.g. Isolate Strep B sample 25292_2#105 have common variants 23S1, 23S3 and GYRA, but replacement of amino acid S by Q in position 17 of the PARC protein sequence

ID | 23S1 | 23S3 | GYRA | PARC | RPOBGBS-1 | RPOBGBS-2 | RPOBGBS-3 | RPOBGBS-4
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
25292_2#105 | 23S1 | 23S3 | GYRA | PARC-Q17S | | | | |

3. **sampleID_drug_cat_alleles.txt** - Gives the alleles for drug categories: EC (macrolides, lincosamides, streptogramins or oxazolidinones), FQ (fluoroquinolones), OTHER (other antibiotics) and TET (tetracyclines)
e.g. Isolate Step B sample 25292_2#105 have erythromycin-resistant 23S1 and 23S3, fluoroquinolone-resistant PARC and GYRA and tetracycline-resistant TETM-1 (allele of TETM)

ID | EC | FQ | OTHER | TET
:---: | :---: | :---: | :---: | :---:
25292_2#105 | 23S1:23S3 | PARC:GYRA | neg | TETM(TETM-1)

#### To run on multiple samples in a directory:
```
cd GBS-Typer-sanger-nf
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'all_samples'
```
This will produce combined tables of all_samples_serotype_res_incidence.txt, all_samples_gbs_res_variants.txt and all_samples_drug_cat_alleles.txt in the 'results' directory that can be identified by sample ID (i.e. the name of the file before _1.fastq.gz or _2.fastq.gz).

#### Additional useful options
    --restyper_min_read_depth   Minimum read depth where mappings to antibiotic resistance genes with fewer reads are excluded. (Default: 30)
    --serotyper_min_read_depth  Minimum read depth where mappings to serotyping genes with fewer reads are excluded. (Default: 30)

    -resume                   Resume if pipeline does not complete.

### Running on farm5
1. Load modules
```
module load ISG/singularity/3.6.4
module load nextflow/20.10.0-5430
```

2. When running (within bsub) add '-profile sanger' as an option, e.g.
```
nextflow run main.nf --reads 'data/sampleID_{1,2}.fastq.gz' --output 'sampleID' -profile sanger
```


### Run unit tests
```
python3 -m pytest
```

### Making a release
Once all changes have been pushed to the main branch, confirm that the CI has run successfully and then execute:
```
release.sh <version number (without v)>
```
This will tag the main branch.
