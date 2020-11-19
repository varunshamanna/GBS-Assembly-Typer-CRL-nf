/*
 * Nextflow pipeline for sero and resistance typing Group B Strep
 *
 */

// Enable DSL 2
nextflow.enable.dsl=2

// Import modules
include {printHelp} from './modules/help.nf'
include {serotyping} from './modules/serotyping.nf'
include {res_typer} from './modules/res_typer.nf'
include {combine_results} from './modules/combine.nf'

// Help message
if (params.help){
    printHelp()
    exit 0
}

// Check parameters
params.reads = ""
params.output = ""

if (params.reads == ""){
    println("Please specify reads with --reads")
    println("Print help with --help")
    System.exit(1)
}

if (params.output == ""){
    println("Please specify and output file with --output")
    println("Print help with --help")
    System.exit(1)
}


// Main workflow
workflow {

    // Create read pairs channel
    Channel.fromFilePairs( params.reads, checkIfExists: true )
        .set { read_pairs_ch }

    // Create new genotypes file
    Channel.fromPath( params.output )
        .set { output_ch }

    println("Using database version ${params.dbversion}")

    // Serotyping
    sero_gene_db = file(params.db_serotyping)
    serotyping(read_pairs_ch, sero_gene_db)

    // Resistance Typer
    res_typer(
        read_pairs_ch,
        file(params.db_restyper),
        file(params.db_argannot),
        file(params.db_resfinder))

    sero_res_ch = serotyping.out.join(res_typer.out)

    // Combine results
    combine_results(sero_res_ch)

    // Combine samples
    combine_results.out
        .collectFile() { item ->
            [ file(params.output), item + '\n']
        }
}
