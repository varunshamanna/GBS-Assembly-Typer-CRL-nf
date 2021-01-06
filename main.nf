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

// Check if reads specified
if (params.reads == ""){
    println("Please specify reads with --reads")
    println("Print help with --help")
    System.exit(1)
}

// Check if output specified
if (params.output == ""){
    println("Please specify and output prefix with --output")
    println("Print help with --help")
    System.exit(1)
}

if (params.other_res_dbs.toString() != 'none'){

    // Check other resistance databases exist
    other_res_db_list = params.other_res_dbs.toString().tokenize(' ')
    for (db in other_res_db_list) {
        other_db = file(db, checkIfExists: true)
    }

    // Check number of resistance databases matches number of minimum coverage and max divergence parameters
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
}

// Create tmp directory if it doesn't already exist
tmp_dir = file('./tmp')
tmp_dir.mkdir()

// Results directory
results_dir = file('./results')

// Output files
params.sero_res_incidence_out = "${params.output}_serotype_res_incidence.txt"
params.variants_out =  "${params.output}_gbs_res_variants.txt"
params.alleles_out = "${params.output}_drug_cat_alleles.txt"

// Resistance mapping with the GBS resistance database
workflow GBS_RES {

    take:
        reads

    main:
        // Split GBS target sequences from GBS resistance database into separate FASTA files per sequence
        split_target_RES_sequences(file(params.gbs_res_typer_db), file(params.gbs_res_targets_db))

        // Get path and name of GBS resistance database
        db_path = file(params.gbs_res_typer_db).getParent()
        db_name = file(params.gbs_res_typer_db).getName()

        // Map genomes to GBS resistance database using SRST2
        srst2_for_res_typing(reads, db_name, db_path, tmp_dir, 'GBS_RES', params.gbs_res_min_coverage, params.gbs_res_max_divergence)

        // Split sam file for each GBS target sequence
        split_target_RES_seq_from_sam_file(srst2_for_res_typing.out.bam_files, file(params.gbs_res_targets_db))

        // Get consensus sequence using freebayes
        freebayes(split_target_RES_seq_from_sam_file.out, split_target_RES_sequences.out, tmp_dir)

    emit:
        freebayes.out.id
}

// Resistance mapping with the other resistance databases
workflow OTHER_RES {

    take:
        reads

    main:
        // Map genomes to resistance database using SRST2
        srst2_for_res_typing(reads, params.other_res_dbs, file('.'), tmp_dir, 'OTHER_RES', params.other_res_min_coverage, params.other_res_max_divergence)

    emit:
        srst2_for_res_typing.out.id
}

// Main Workflow
workflow {

    // Create new genotypes file
    Channel.fromPath( params.output )
        .set { output_ch }

    // Create read pairs channel
    Channel.fromFilePairs( params.reads, checkIfExists: true )
        .set { read_pairs_ch }

    main:

        // Serotyping Process
        serotyping(read_pairs_ch, file(params.serotyping_db))

        // Resistance Mapping Workflows
        GBS_RES(read_pairs_ch)

        if (params.other_res_dbs != 'none'){
            OTHER_RES(read_pairs_ch)
            id_ch = GBS_RES.out.join(OTHER_RES.out)
        } else {
            id_ch = GBS_RES.out
        }

        // Once GBS or both resistance workflows are complete, trigger resistance typing
        res_typer(id_ch, tmp_dir)

        // Combine serotype and resistance type results for each sample
        sero_res_ch = serotyping.out.join(res_typer.out)
        combine_results(sero_res_ch)

        // Combine samples and output results files
        combine_results.out.sero_res_incidence
            .collectFile(name: file("${results_dir}/${params.sero_res_incidence_out}"), keepHeader: true)

        combine_results.out.res_alleles
            .collectFile(name: file("${results_dir}/${params.alleles_out}"), keepHeader: true)

        combine_results.out.res_variants
            .collectFile(name: file("${results_dir}/${params.variants_out}"), keepHeader: true)
}
