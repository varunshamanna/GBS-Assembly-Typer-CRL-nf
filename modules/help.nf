
def printHelp() {
  log.info"""
  Usage:
    nextflow run main.nf --reads [paired reads] --output [results directory]

  Description:
    Characterisation of Group B Strep.

    ## Additional options
    ### Inputs
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
  """.stripIndent()
}
