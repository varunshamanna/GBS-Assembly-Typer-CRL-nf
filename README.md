[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/GBS-Typer-sanger-nf/blob/main/LICENSE)

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/sanger-bentley-group/GBS-Typer-sanger-nf)](https://github.com/sanger-bentley-group/GBS-Typer-sanger-nf/releases)

![GitHub Workflow Status](https://github.com/sanger-bentley-group/GBS-Typer-sanger-nf/actions/workflows/test.yml/badge.svg)

# GBS-Typer-sanger-nf
The GBS Typer is for characterising Group B Strep by serotyping, resistance typing, MLST, surface protein typing and penicillin-binding protein typing. It has been adapted from [Ben Metcalf's GBS Typer pipeline](https://github.com/BenJamesMetcalf/GBS_Scripts_Reference) in Nextflow for portability and reproducibility.

## Contents
- [ Running the pipeline ](#local)
    - [ Installation ](#localinstall)
    - [ Usage ](#localusage)
    - [ Running the pipeline on Sanger farm ](#sanger)
- [ Outputs ](#outputs)
    - [ Main Report ](#main)
    - [ Other Reports ](#other)
- [ Additional options ](#additional)
- [ Advanced ](#advanced)
    - [ Pencillin-binding protein Typing Workflow ](#pbp)
    - [ Other examples of running pipelines ](#examples)
    - [ Troubleshooting for errors](#errors)
    - [ Clean Up ](#cleanup)
    - [ Software dependencies ](#dependencies)
- [ Acknowledge us ](#acknowledge)
- [ Reporting issues ](#issues)

<a name="local"></a>
## Running the pipeline
- Running the pipeline requires an internet connection
- Currently it supports only paired-end reads

<a name="localinstall"></a>
### Installation
GBS Typer relies on Nextflow and Docker.
Download:
1. [Docker](https://www.docker.com/).
2. [Nextflow](https://www.nextflow.io/).
3. Move `nextflow` to your PATH
```
mv nextflow /usr/local/bin/
```
4. Clone repository:
```
git clone https://github.com/sanger-pathogens/GBS-Typer-sanger-nf.git
cd GBS-Typer-sanger-nf
```

<a name="localusage"></a>
### Usage

- To run on one sample (called `sampleID` in directory `data`)
```
nextflow run main.nf --reads data/sampleID_{1,2}.fastq.gz --results_dir my_results
```

- To run on multiple samples in a directory (called `data` and a results directory called `my_results`)
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --results_dir my_results
```

This will create a results directory (here called `my_results`) containing output files listed below.

<a name="sanger"></a>
### Running the pipeline on Sanger farm
Follow [these](README_sanger.md) instructions.

<a name="outputs"></a>
## Outputs

<a name="main"></a>
### Main Report
Running the command with will generate the main report `gbs_typer_report.txt`. This will include the serotype, MLST type, allelic frequencies from MLST, resistance gene incidence, surface protein types and GBS-specific resistance variants. You can find the description for each of the columns in the report [here](https://docs.google.com/spreadsheets/d/1R5FFvACC3a6KCKkTiluhTj492-4cCe74HcCoklqX-X0/edit?usp=sharing) where the `category` column is `in_silico_analysis`.

<a name="other"></a>
### Other Reports
1. **serotype_res_incidence.txt**

Gives the serotype and presence/absence (i.e. pos/neg) of antibiotic resistance genes, e.g. Isolate Strep B sample 25292_2#105 has serotype II and have genes: 23S1, 23S3, GYRA, lsaC and tetM

Sample_id | Serotype | 23S1 | 23S3 | CAT | ermB | ermT | FOSA | GYRA | lnuB | lsaC | mefA | MPHC | MSRA | msrD | PARC | RPOBGBS-1 | RPOBGBS-2 | RPOBGBS-3 | RPOBGBS-4 | SUL2 | tetB | tetL | tetM | tetO | tetS
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
25292_2#105 | II | pos | pos | neg | neg | neg | neg | pos | neg | pos | neg | neg | neg | neg | neg | neg | neg | neg | neg | neg | neg | neg | pos | neg | neg

2. **gbs_res_variants.txt**

Gives the SNP variants for GBS-specific resistance genes
e.g. Isolate Strep B sample 25292_2#105 have common variants 23S1, 23S3 and GYRA (shown with *), but replacement of amino acid S by Q in position 17 of the PARC protein sequence

Sample_id | 23S1 | 23S3 | GYRA | PARC | RPOBGBS-1 | RPOBGBS-2 | RPOBGBS-3 | RPOBGBS-4
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
25292_2#105 | * | * | * | Q17S | | | | |

3. **drug_cat_alleles_variants.txt**

Gives the GBS-specific variants and other resistance genes and alleles for drug categories: EC (macrolides, lincosamides, streptogramins or oxazolidinones), FQ (fluoroquinolones), OTHER (other antibiotics) and TET (tetracyclines)
e.g. Isolate Step B sample 25292_2#105 have GBS-specific variants: erythromycin-resistant 23S1 and 23S3, fluoroquinolone-resistant PARC and GYRA, and other resistance allele tetracycline-resistant tet(M)_1 of gene tet(M) (as specified by gene[allele])

Sample_id | EC | FQ | OTHER | TET
:---: | :---: | :---: | :---: | :---:
25292_2#105 | 23S1:23S3 | PARC-Q17S:GYRA | neg | tet(M)[tet(M)_1]

4. **new_mlst_alleles.log**

Indicates whether new MLST alleles have been found for each sample (where there are mismatches with sufficient read depth at least the value specified --mlst_min_read_depth [Default: 30]).

5. FASTA file **new_mlst_alleles.fasta** and a pileup file **new_mlst_pileup.txt**

If `new_mlst_alleles.log` includes "New MLST alleles found."

6. **existing_sequence_types.txt**

For other samples that have no new MLST alleles and only have existing sequence types. If "None found" for a sample then no sequence types were found (with sufficient read depth).

7. **surface_protein_incidence.txt**

This shows the incidence of different surface protein alleles in the Strep B sample(s), e.g.

Sample_id | alp1 | alp2/3 | alpha | hvgA | PI1 | PI2A1 | PI2A2 | PI2B | rib | srr1 | srr2
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
26189_8#338 | neg | pos | neg | neg | pos | pos | neg | neg | neg | pos | neg

8. **surface_protein_variants.txt**

This shows all the surface proteins in the Strep B sample(s), e.g.

Sample_id | ALPH | hvgA | PILI | SRR
:---: | :---: | :---: | :---: | :---:
26189_8#338 | alp2/3 | neg | PI1:PI2A1 | srr1

<a name="additional"></a>
## Additional options
### Inputs
    --contigs                       Path of file containing FASTA contigs. Only use when --run_pbptyper is specified. (Use wildcard '*' to specify multiple files, e.g. 'data/*.fa')
    --db_version                    Database version. (Default: latest in `db` directory)
    --other_res_dbs                 Path to other resistance reference database. Must be FASTA format. (Default: ResFinder)

### Parameters
    --gbs_res_min_coverage          Minimum coverage for mapping to the GBS resistance database. (Default: 99.9)
    --gbs_res_max_divergence        Maximum divergence for mapping to the GBS resistance database. (Default: 5, i.e. report only hits with <5% divergence)
    --other_res_min_coverage        Minimum coverage for mapping to other resistance reference database(s). (Default: 70)
    --other_res_max_divergence      Maximum divergence for mapping to other resistance reference database. (Default: 30, i.e. report only hits with <30% divergence)
    --restyper_min_read_depth       Minimum read depth where mappings to antibiotic resistance genes with fewer reads are excluded. (Default: 30)
    --serotyper_min_read_depth      Minimum read depth where mappings to serotyping genes with fewer reads are excluded. (Default: 30)

### Other Workflow Options
    --run_sero_res                  Run the main serotyping and resistance workflows. (Default: true)
                                    Use '--run_sero_res false' to override the default.
    --run_surfacetyper              Run the surface protein typing workflow. (Default: false)
    --run_mlst                      Run the MLST workflow to query existing sequence types and new MLST alleles. (Default: true)
    --run_pbptyper                  Run the PBP (Penicillin binding protein) allele typer workflow. Must also specify --contigs input. (Default: true)

### Other Parameters
    --mlst_min_read_depth           Minimum read depth where mappings to alleles in MLST with fewer reads are excluded. Only operational with --run_mlst. (Default: 30)
    --pbp_frac_align_threshold      Minimum fraction of sequence alignment length of PBP gene. Only operational with --run_pbptyper. (Default: 0.5)
    --pbp_frac_identity_threshold   Minimum fraction of alignment identity between PBP genes and assemblies. Only operational with --run_pbptyper. (Default: 0.5)
    --surfacetyper_min_coverage     Minimum coverage for mapping to the GBS surface protein database. Only operational with --run_surfacetyper. (Default: 70)
    --surfacetyper_max_divergence   Maximum divergence for mapping to the GBS surface protein database. Only operational with --run_surfacetyper. (Default: 8, i.e. report only hits with <8% divergence)
    --surfacetyper_min_read_depth   Minimum read depth for surface protein typing workflow. Only operational with --run_surfacetyper. (Default: 30)

<a name="advanced"></a>
## Advanced

<a name="pbp"></a>
### PBP (Penicillin-binding protein) Typing Workflow
To enable the PBP typing workflow provide the **--run_pbptyper** command line argument and specify the contig FASTA files using and **--contigs**:
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --results_dir my_results --run_pbptyper --contigs 'data/*.fa'
```

If existing PBP alleles are found, a tab-delimited file is created in the results directory. The file contains the sample IDs (that are determined from the contig FASTA file names e.g. 25292_2#85 from data/25292_2#85.fa), the contig identifiers with start, end and forward(+)/reverse(-) positions, and the PBP allele identifier.

Sample_id | Contig | PBP_allele
:---: | :---: | :---:
25292_2#85 | .25292_2_85.9:39495-40455(+) | 2\|\|GBS_1A
26077_6#118 | .26077_6_118.11:39458-40418(+) | 1\|\|GBS_1A

If a new PBP allele is found in a sample, a FASTA file of amino acids is created in the 'results' directory. For example, if contig .25292_2_85.9 from sample 25292_2#85 contained a new PBP-1A allele, then .25292_2_85.9_GBS1A-1_PBP_new_allele.faa is generated with contents:

    >.26077_6_118.11:39458-40418(+)
    DIYNSDTYIAYPNNELQIASTIMDATNGKVIAQLGGRHQNENISFGTNQSVLTDRDWGST
    MKPISAYAPAIDSGVYNSTGQSLNDSVYYWPGTSTQLYDWDRQYMGWMSMQTAIQQSRNV
    PAVRALEAAGLDEAKSFLEKLGIYYPEMNYSNAISSNNSSSDAKYGASSEKMAAAYSAFA
    NGGTYYKPQYVNKIEFSDGTNDTYAASGSRAMKETTAYMMTDMLKTVLTFGTGTKAAIPG
    VAQAGKTGTSNYTEDELAKIEATTGIYNSAVGTMAPDENFVGYTSKYTMAIWTGYKNRLT
    PLYGSQLDIATEVYRAMMSY


<a name="examples"></a>
### Other examples of running the pipeline
It is recommended you use the default parameters for specifying other resistance databases. However, to use different or multiple resistance databases with the GBS-specific resistance database, e.g. ARG-ANNOT and ResFinder in the `db/0.2.1` directory, both with a minimum coverage of 70 and maximum divergence of 30:
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --results_dir my_results --other_res_dbs 'db/0.2.1/ARGannot-DB/ARG-ANNOT.fasta db/0.2.1/ResFinder-DB/ResFinder.fasta' --other_res_min_coverage '70 70' --other_res_max_divergence '30 30'
```
To run **only** the surface protein typing workflow, use:
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --results_dir my_results --run_sero_res false --run_mlst false
```
To run **only** the MLST workflow, use:
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --results_dir my_results --run_sero_res false --run_surfacetyper false
```
To run **only** the serotyping and resistance typing workflows, use:
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --results_dir my_results --run_surfacetyper false --run_mlst false
```
To run **only** the PBP typing workflow, use:
```
nextflow run main.nf --results_dir my_results --run_sero_res false --run_surfacetyper false --run_mlst false --run_pbptyper --contigs 'data/*.fa'
```
Note: The **--reads** parameter is not needed for the PBP typing workflow.

<a name="errors"></a>
### Troubleshooting for errors
It is possible that the pipeline may not complete successfully due to issues with input files and/or individual steps of the pipeline. To troubleshoot potential issues, you can use the `-with-trace` parameter e.g. `nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --results_dir my_results -with-trace` to create a trace file in the current directory while the pipeline is running. The trace file will provide the status of the current directory that provides the status and the location of the log files for each sample and step. For example, if the step failed in `[00/8803ea]`, the command and the error from this step can be viewed by `cat work/00/8803ea01f8bc0fb43f68335d82831f/.command.log` (Hint: while typing `work/00/8803ea`, the complete file path can be completed with the TAB button.)

<a name="cleanup"></a>
### Clean up

The `work` directory keeps the intermediate files of the pipeline. You can use `nextflow clean [run_name|session_id]` to clean up or `rm -rf work` only when no other pipelines are running (more dangerous). The run names and session ids can be found by using `nextflow log -q`

<a name="dependencies"></a>
### Software dependencies
All pipeline dependencies are built into the [quay.io dependencies image](https://quay.io/repository/sangerpathogens/gbs-typer-sanger-nf), used by the pipeline.
The current program versions in this image are as follows:

Program | Version
:---: | :---:
bedtools | 2.29.2
biopython | 1.78
bowtie | 2.2.9
freebayes | 1.3.3+
prodigal | 1:2.6.3
python 2 | 2.7
python 3 | 3.8
samtools | 0.1.18
srst2 | 0.2.0

<a name="acknowledge"></a>
## Acknowledge us
If you use this software, please acknowledge the contributors as appropriate.

<a name="issues"></a>
## Reporting issues
If any questions or problems, please post them under [Issues](https://github.com/sanger-pathogens/GBS-Typer-sanger-nf/issues).
