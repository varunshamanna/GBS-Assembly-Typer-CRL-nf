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

3. **sampleID_drug_cat_alleles_variants.txt** - Gives the GBS-specific variants and other resistance genes and alleles for drug categories: EC (macrolides, lincosamides, streptogramins or oxazolidinones), FQ (fluoroquinolones), OTHER (other antibiotics) and TET (tetracyclines)
e.g. Isolate Step B sample 25292_2#105 have GBS-specific variants: erythromycin-resistant 23S1 and 23S3, fluoroquinolone-resistant PARC and GYRA, and other resistance allele tetracycline-resistant tet(M)_1 of gene tet(M) (as specified by gene[allele])

ID | EC | FQ | OTHER | TET
:---: | :---: | :---: | :---: | :---:
25292_2#105 | 23S1:23S3 | PARC-Q17S:GYRA | neg | tet(M)[tet(M)_1]

#### To run on multiple samples in a directory:
```
cd GBS-Typer-sanger-nf
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'output_file_prefix'
```
This will produce combined tables of output_file_prefix_serotype_res_incidence.txt, output_file_prefix_gbs_res_variants.txt and output_file_prefix_drug_cat_alleles_variants.txt in the 'results' directory that can be identified by sample ID (i.e. the name of the file before _1.fastq.gz or _2.fastq.gz).

#### Additional useful options
##### Inputs
    --db_version                Database version. (Default: 0.0.2)
    --other_res_dbs             Paths to other resistance reference database(s). Must be FASTA format. Specify 'none' to omit using other resistance databases. (Default: 'db/0.0.2/ARGannot-DB/ARG-ANNOT.fasta db/0.0.2/ResFinder-DB/ResFinder.fasta')

##### Outputs
    --results_dir               Results directory for output files. (Default: './results')
    --tmp_dir                   Temporary directory for intermediate files. (Default: './tmp')

##### Parameters
    --gbs_res_min_coverage      Minimum coverage for mapping to the GBS resistance database. (Default: 99.9)
    --gbs_res_max_divergence    Maximum divergence for mapping to the GBS resistance database. (Default: 5)
    --other_res_min_coverage    Minimum coverage for mapping to other resistance reference database(s). Number of values must equal the number of resistance reference database files and must correspond to the order specified in --other_res_dbs. (Default: '70 70')
    --other_res_max_divergence  Maximum divergence for mapping to other resistance reference database(s). Number of values must equal the number of resistance reference database files and must correspond to the order specified in --other_res_dbs. (Default: '30 30')
    --restyper_min_read_depth   Minimum read depth where mappings to antibiotic resistance genes with fewer reads are excluded. (Default: 30)
    --serotyper_min_read_depth  Minimum read depth where mappings to serotyping genes with fewer reads are excluded. (Default: 30)

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
