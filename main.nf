/*
 * Nextflow pipeline for sero and resistance typing Group B Strep
 *
 */

// Enable DSL 2
nextflow.enable.dsl=2

// Import modules
include {printHelp} from './modules/help.nf'
include {serotyping} from './modules/serotyping.nf'
include {srst2_for_res_typing; split_target_RES_seq_from_sam_file; split_target_RES_sequences; freebayes} from './modules/res_alignments.nf'
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
    println("Please specify and output prefix with --output")
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

// Create tmp directory
tmp_dir = file('./tmp')
tmp_dir.mkdir()

// Output files
params.sero_res_incidence_out = "${params.output}_serotype_res_incidence.txt"
params.variants_out =  "${params.output}_gbs_res_variants.txt"
params.alleles_out = "${params.output}_drug_cat_alleles.txt"

// Resistance typing with the RES database
workflow RES {

    take:
        reads

    main:
        res_targets_file = file(params.db_res_targets)

        // Copy to tmp directory
        gbs_db = file(params.db_gbs_res_typer)
        gbs_db.copyTo(tmp_dir)

        split_target_RES_sequences(file(params.db_gbs_res_typer), res_targets_file)

        srst2_for_res_typing(reads, params.db_gbs_res_typer, tmp_dir, 'RES', 99.9, 5)
        split_target_RES_seq_from_sam_file(srst2_for_res_typing.out.bam_files, res_targets_file)

        freebayes(split_target_RES_seq_from_sam_file.out, split_target_RES_sequences.out, tmp_dir)

    emit:
        freebayes.out.id
}

// Resistance typing with the other resistance databases
workflow OTHER_RES {

    take:
        reads

    main:
        // Copy to tmp directory
        db_list = params.other_res_dbs.tokenize(' ')
        for (db in db_list) {
            other_db = file(db)
            other_db.copyTo(tmp_dir)
        }

        srst2_for_res_typing(reads, params.other_res_dbs, tmp_dir, 'OTHER_RES', 70, 30)

    emit:
        srst2_for_res_typing.out.id
}

workflow {

    // Create new genotypes file
    Channel.fromPath( params.output )
        .set { output_ch }

    // Create read pairs channel
    Channel.fromFilePairs( params.reads, checkIfExists: true )
        .set { read_pairs_ch }

    main:
        serotyping(read_pairs_ch, file(params.db_serotyping))

        RES(read_pairs_ch)
        OTHER_RES(read_pairs_ch)
        id_ch = RES.out.join(OTHER_RES.out)
        res_typer(id_ch, tmp_dir)

        // Combine results
        sero_res_ch = serotyping.out.join(res_typer.out)
        combine_results(sero_res_ch)

        // Combine samples
        combine_results.out.sero_res_incidence
            .collectFile(name: file(params.sero_res_incidence_out), keepHeader: true)

        combine_results.out.res_alleles
            .collectFile(name: file(params.alleles_out), keepHeader: true)

        combine_results.out.res_variants
            .collectFile(name: file(params.variants_out), keepHeader: true)
}
