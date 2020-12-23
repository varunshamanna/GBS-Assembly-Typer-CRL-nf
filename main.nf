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

if (params.other_res_dbs.toString().tokenize(' ').size() != params.other_res_min_coverage.toString().tokenize(' ').size()){
    println("Number of --other_res_min_coverage values is not equal to the number of --other_res_dbs files")
    println("Please specify an equal number of values.")
    System.exit(1)
}

if (params.other_res_dbs.toString().tokenize(' ').size() != params.other_res_max_divergence.toString().tokenize(' ').size()){
    println("Number of --other_res_max_divergence values is not equal to the number of --other_res_dbs files")
    println("Please specify an equal number of values.")
    System.exit(1)
}

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
        res_targets_file = file(params.gbs_res_targets_db)

        // Copy to tmp directory
        gbs_db = file(params.gbs_res_typer_db)
        gbs_db.copyTo(tmp_dir)

        split_target_RES_sequences(file(params.gbs_res_typer_db), res_targets_file)

        srst2_for_res_typing(reads, params.gbs_res_typer_db, tmp_dir, 'RES', params.gbs_res_min_coverage, params.gbs_res_max_divergence)
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
        db_list = params.other_res_dbs.toString().tokenize(' ')
        for (db in db_list) {
            other_db = file(db)
            other_db.copyTo(tmp_dir)
        }

        srst2_for_res_typing(reads, params.other_res_dbs, tmp_dir, 'OTHER_RES', params.other_res_min_coverage, params.other_res_max_divergence)

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
        serotyping(read_pairs_ch, file(params.serotyping_db))

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
