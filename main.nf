/*
 * Nextflow pipeline for sero and resistance typing Group B Strep
 *
 */

// Enable DSL 2
nextflow.enable.dsl=2

// Import modules
include {printHelp} from './modules/help.nf'
include {serotyping} from './modules/serotyping.nf'
include {res_typer; srst2_for_res_typing; split_target_RES_seq_from_sam_file; split_target_RES_sequences; freebayes} from './modules/res_typer.nf'
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
params.db_res_targets = "$params.db/GBS_resTyper_Gene-DB/seqs_of_interest.txt"
params.db_argannot = "$params.db/ARGannot-DB/ARGannot_r1.fasta"
params.db_resfinder = "$params.db/ResFinder-DB/ResFinder.fasta"

// Resistance typing with the RES database
workflow RES {

    take:
        reads
        rt_RES

    main:
        res_typer_gene_db = file(params.db_gbs_res_typer)
        res_targets_file = file(params.db_res_targets)

        split_target_RES_sequences(res_typer_gene_db, res_targets_file)

        srst2_for_res_typing(reads, res_typer_gene_db, 'RES', 99.9, 5)
        split_target_RES_seq_from_sam_file(srst2_for_res_typing.out.bam_files, res_targets_file)

        freebayes(split_target_RES_seq_from_sam_file.out, split_target_RES_sequences.out)

    emit:
        rt_RES = srst2_for_res_typing.out.genes_files.join(freebayes.out)
}

// Resistance typing with the ARGANNOT database
workflow ARGANNOT {

    take:
        reads
        rt_ARG

    main:
        argannot_db = file(params.db_argannot)
        srst2_for_res_typing(reads, argannot_db, 'ARG', 70, 30)

    emit:
        rt_ARG = srst2_for_res_typing.out.genes_files
}

// Resistance typing with the ResFinder database
workflow ResFinder {

    take:
        reads
        rt_ResFinder

    main:
        resfinder_db = file(params.db_resfinder)
        srst2_for_res_typing(reads, resfinder_db, 'ARG', 70, 30)

    emit:
        rt_ResFinder = srst2_for_res_typing.out.genes_files
}

workflow {

    // Create new genotypes file
    Channel.fromPath( params.output )
        .set { output_ch }

    // Create read pairs channel
    Channel.fromFilePairs( params.reads, checkIfExists: true )
        .set { read_pairs_ch }

    // Create empty resistance typing channels
    rt_RES = channel.empty()
    rt_ARG = channel.empty()
    rt_ResFinder = channel.empty()

    main:
        serotyping(read_pairs_ch, file(params.db_serotyping))

        if( params.all | params.custom_res)
            RES(read_pairs_ch, rt_RES)

        if( params.all | params.argannot)
            ARGANNOT(read_pairs_ch, rt_ARG)

        if( params.all | params.resfinder)
            ResFinder(read_pairs_ch, rt_ResFinder)

        res_typer_ch = rt_ARG.join(rt_ResFinder.join(rt_RES))
        res_typer(res_typer_ch)




}
