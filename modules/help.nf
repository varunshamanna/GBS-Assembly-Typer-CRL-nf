
def printHelp() {
  log.info"""
  Usage:
    nextflow run main.nf --reads [paired reads] --output [output prefix]
  Description:
    Serotyping and resistance typing of Group Strep B.
  Additional options:
    --restyper_min_read_depth   Minimum read depth where mappings to antibiotic resistance genes with fewer reads are excluded. (Default: 30)
    --serotyper_min_read_depth  Minimum read depth where mappings to serotyping genes with fewer reads are excluded. (Default: 30)
  """.stripIndent()
}
