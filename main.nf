/*
 * Nextflow pipeline for sero and resistance typing Group B Strep
 *
 */

// Enable DSL 2
nextflow.enable.dsl=2

// Import modules
include {printHelp} from './modules/help.nf'
include {serotyping} from './modules/serotyping.nf'

// Help message
if (params.help){
    printHelp()
    exit 0
}

// Check parameters
params.reads = ""
params.contigs = ""

if (params.reads == ""){
    println("Please specify reads with --reads")
    println("Print help with --help")
    System.exit(1)
}

if (params.contigs == ""){
    println("Please specify reads with --contigs")
    println("Print help with --help")
    System.exit(1)
}

// Import Serotyping DB
params.db = "./db/$params.dbversion"
params.db_serotyping = "$params.db/GBS_seroT_Gene-DB/GBS_seroT_Gene-DB_Final.fasta"

// Main workflow
workflow {

    // Create read pairs channel
    Channel.fromFilePairs( params.reads )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .set { read_pairs_ch }

    // Create contigs channel
    Channel.fromPath( params.contigs )
        .ifEmpty { error "Cannot find any contigs matching: ${params.contigs}" }
        .set { contigs_ch }

    // Serotyping
    sero_gene_db = file(params.db_serotyping)
    serotyping(read_pairs_ch, sero_gene_db)

}
