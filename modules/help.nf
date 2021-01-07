
def printHelp() {
  log.info"""
  Usage:
    nextflow run main.nf --reads [paired reads] --output [output prefix]

  Description:
    Serotyping and resistance typing of Group B Strep.

  Additional options:

    Inputs:
      --db_version                  Database version. (Default: 0.0.2)
      --other_res_dbs               Paths to other resistance reference database(s). Must be FASTA format.
                                    Specify 'none' to omit using other resistance databases.
                                    (Default: 'db/0.0.2/ARGannot-DB/ARG-ANNOT.fasta db/0.0.2/ResFinder-DB/ResFinder.fasta')

    Outputs:
      --results_dir                 Results directory for output files. (Default: './results')
      --tmp_dir                     Temporary directory for intermediate files. (Default: './tmp')

    Parameters:
      --gbs_res_min_coverage        Minimum coverage for mapping to the GBS resistance database. (Default: 99.9)
      --gbs_res_max_divergence      Maximum divergence for mapping to the GBS resistance database. (Default: 5)
      --other_res_min_coverage      Minimum coverage for mapping to other resistance reference database(s).
                                    Number of values must equal the number of resistance reference database files
                                    and must correspond to the order specified in --other_res_dbs. (Default: '70 70')
      --other_res_max_divergence    Maximum divergence for mapping to other resistance reference database(s).
                                    Number of values must equal the number of resistance reference database files
                                    and must correspond to the order specified in --other_res_dbs. (Default: '30 30')
      --restyper_min_read_depth     Minimum read depth where mappings to antibiotic resistance genes with fewer
                                    reads are excluded. (Default: 30)
      --serotyper_min_read_depth    Minimum read depth where mappings to serotyping genes with fewer reads are excluded.
                                    (Default: 30)
  """.stripIndent()
}
