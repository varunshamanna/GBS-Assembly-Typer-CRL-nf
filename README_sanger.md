## Running this pipeline on Sanger farm

**1. Clone repository**
```
git clone https://github.com/sanger-pathogens/GBS-Typer-sanger-nf.git
cd GBS-Typer-sanger-nf
```

**2. Load modules**
```
module load ISG/singularity
module load nextflow
```

**3. If running on farm5, you will need to set the http/https proxy**
```
export http_proxy=http://wwwcache.sanger.ac.uk:3128
export https_proxy=http://wwwcache.sanger.ac.uk:3128
```

**4. Run using bsub**

-  For a single sample
```
bsub -G <your_team> -J <job_name> -o %J.out -e %J.err -R "select[mem>1000] rusage[mem=1000]" -M1000 "nextflow run main.nf --reads data/sampleID_{1,2}.fastq.gz --results_dir my_results -profile sanger,lsf"
```

- For multiple samples
```
bsub -G <your_team> -J <job_name> -o %J.out -e %J.err -R "select[mem>1000] rusage[mem=1000]" -M1000 "nextflow run main.nf --reads data/*_{1,2}.fastq.gz --results_dir my_results -profile sanger,lsf"
```
Specifying `-profile sanger,lsf` will instruct Nextflow to run tasks as separate LSF jobs in parallel and will instruct the pipeline to build a local Singularity image from the [quay.io Docker image](https://quay.io/repository/sangerpathogens/gbs-typer-sanger-nf).

If you are processing many samples and would like to speed up the pipeline, you can increase the `queue_size` (default 100 i.e. a maximum 100 jobs running concurrently). For example, to increase the number of jobs running concurrently to 500, set `--queue_size 500`. (Note, the number of jobs may already be limited by your HPC's available resources.)
