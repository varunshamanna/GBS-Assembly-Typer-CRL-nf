/*
 * Nextflow pipeline for sero and resistance typing Group B Strep
 *
 */

// Enable DSL 2
nextflow.enable.dsl=2

// Import modules
include {printHelp} from './modules/help.nf'
include {serotyping} from './modules/serotyping.nf'
include {res_typer; srst2_for_res_typing} from './modules/res_typer.nf'
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

// Import databases
params.db = "./db/$params.dbversion"
params.db_serotyping = "$params.db/GBS_seroT_Gene-DB/GBS_seroT_Gene-DB_Final.fasta"
params.db_gbs_res_typer = "$params.db/GBS_resTyper_Gene-DB/GBS_Res_Gene-DB_Final.fasta"
params.db_argannot = "$params.db/ARGannot-DB/ARGannot_r1.fasta"
params.db_resfinder = "$params.db/ResFinder-DB/ResFinder.fasta"


// Main workflow
workflow {

    // Create read pairs channel
    Channel.fromFilePairs( params.reads, checkIfExists: true )
        .set { read_pairs_ch }

    // Create new genotypes file
    Channel.fromPath( params.output )
        .set { output_ch }

    // Serotyping
    // sero_gene_db = file(params.db_serotyping)
    // serotyping(read_pairs_ch, sero_gene_db)

    // Resistance Typer
    res_typer_gene_db = file(params.db_gbs_res_typer)
    srst2_for_res_typing(read_pairs_ch, res_typer_gene_db, 'RES', 99.9, 5)

    argannot_db = file(params.db_argannot)
    srst2_for_res_typing(read_pairs_ch, argannot_db, 'ARG', 70, 30)

    resfinder_db = file(params.db_resfinder)
    srst2_for_res_typing(read_pairs_ch, resfinder_db, 'ARG', 70, 30)
    
    //res_typer(read_pairs_ch, res_typer_gene_db, argannot_db, resfinder_db)

    //sero_res_ch = serotyping.out.join(res_typer.out)

    // Combine results
    //combine_results(sero_res_ch)

    // Combine samples
    //combine_results.out
    //    .collectFile() { item ->
    //        [ file(params.output), item + '\n']
    //    }
}
