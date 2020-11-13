
def printHelp() {
  log.info"""
  Usage:
    nextflow run main.nf --reads [paired reads] --contigs [contig fasta]
  Description: TODO
  """.stripIndent()
}
