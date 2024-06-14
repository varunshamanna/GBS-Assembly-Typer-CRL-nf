/*
 * Nextflow pipeline for sero and resistance typing Group B Strep
 *
 */

// Enable DSL 2
nextflow.enable.dsl=2


//assembly and qc modules
include { FILE_VALIDATION; PREPROCESS; READ_QC } from "$projectDir/modules/preprocess"
include { ASSEMBLY_UNICYCLER; ASSEMBLY_SHOVILL; ASSEMBLY_ASSESS; ASSEMBLY_QC } from "$projectDir/modules/assembly"
include { GET_REF_GENOME_BWA_DB; MAPPING; SAM_TO_SORTED_BAM; SNP_CALL; HET_SNP_COUNT; MAPPING_QC } from "$projectDir/modules/mapping"
include { GET_KRAKEN2_DB; TAXONOMY; TAXONOMY_QC } from "$projectDir/modules/taxonomy"
include { OVERALL_QC } from "$projectDir/modules/overall_qc"
include { GENERATE_SAMPLE_REPORT; GENERATE_OVERALL_REPORT } from "$projectDir/modules/output"
// Import modules
include {printHelp} from "$projectDir/modules/help.nf"
include {serotyping} from "$projectDir/modules/serotyping.nf"
include {srst2_for_res_typing; split_target_RES_seq_from_sam_file; split_target_RES_sequences; freebayes} from "$projectDir/modules/res_alignments.nf"
include {res_typer} from "$projectDir/modules/res_typer.nf"
include {surface_typer} from "$projectDir/modules/surface_typer.nf"
include {srst2_for_mlst; get_mlst_allele_and_pileup} from "$projectDir/modules/mlst.nf"
include {get_pbp_genes; get_pbp_alleles} from "$projectDir/modules/pbp_typer.nf"
include {finalise_sero_res_results; finalise_surface_typer_results; finalise_pbp_existing_allele_results; combine_results} from "$projectDir/modules/combine.nf"
include {get_version} from "$projectDir/modules/version.nf"



// Help message
if (params.help){
    printHelp()
    exit 0
}

if (params.reads == ""){
        println("Please specify reads with --reads.")
        println("Print help with nextflow main.nf --help")
        System.exit(1)
}

// Check if results_dir specified
if (params.output == ""){
    println("Please specify the results directory with --output.")
    println("Print help with nextflow main.nf --help")
    System.exit(1)
}

// Check parameters are within range
if (params.gbs_res_min_coverage < 0 | params.gbs_res_min_coverage > 100){
    println("--gbs_res_min_coverage value not in range. Please specify a value between 0 and 100.")
    System.exit(1)
}

if (params.gbs_res_max_divergence < 0 | params.gbs_res_max_divergence > 100){
    println("--gbs_res_max_divergence value not in range. Please specify a value between 0 and 100.")
    System.exit(1)
}

other_res_min_coverage_list = params.other_res_min_coverage.toString().tokenize(' ')
for (other_res_min_coverage in other_res_min_coverage_list){
    if (other_res_min_coverage.toDouble() < 0 | other_res_min_coverage.toDouble() > 100){
        println("--other_res_min_coverage value(s) not in range. Please specify a value between 0 and 100.")
        System.exit(1)
    }
}

other_res_max_divergence_list = params.other_res_max_divergence.toString().tokenize(' ')
for (other_res_max_divergence in other_res_max_divergence_list){
    if (other_res_max_divergence.toDouble() < 0 | other_res_max_divergence.toDouble() > 100){
        println("--other_res_max_divergence value(s) not in range. Please specify a value between 0 and 100.")
        System.exit(1)
    }
}

if (params.restyper_min_read_depth < 0){
    println("--restyper_min_read_depth value not in range. Please specify a value of 0 or above.")
    System.exit(1)
}

if (params.serotyper_min_read_depth < 0){
    println("--serotyper_min_read_depth value not in range. Please specify a value of 0 or above.")
    System.exit(1)
}

if (params.mlst_min_coverage < 0 | params.mlst_min_coverage > 100){
    println("--mlst_min_coverage value not in range. Please specify a value between 0 and 100.")
    System.exit(1)
}

if (params.mlst_min_read_depth < 0){
    println("--mlst_min_read_depth value not in range. Please specify a value of 0 or above.")
    System.exit(1)
}

if (params.surfacetyper_min_coverage < 0 | params.surfacetyper_min_coverage > 100){
    println("--surfacetyper_min_coverage value not in range. Please specify a value between 0 and 100.")
    System.exit(1)
}

if (params.surfacetyper_max_divergence < 0 | params.surfacetyper_max_divergence > 100){
    println("--surfacetyper_max_divergence value not in range. Please specify a value between 0 and 100.")
    System.exit(1)
}

if (params.surfacetyper_min_read_depth < 0){
    println("--surfacetyper_min_read_depth value not in range. Please specify a value of 0 or above.")
    System.exit(1)
}

// Create results directory if it doesn't already exist
output_dir = file(params.output)

if (output_dir.exists()){
    println("The output directory already exists. The results will be written in the existing one.")
} else {
    output_dir.mkdir()
}

// Resistance mapping with the GBS resistance database
workflow GBS_RES {

    take:
        reads

    main:

        gbs_res_typer_db = file(params.gbs_res_typer_db, checkIfExists: true)
        gbs_res_targets_db = file(params.gbs_res_targets_db, checkIfExists: true)

        // Split GBS target sequences from GBS resistance database into separate FASTA files per sequence
        split_target_RES_sequences(gbs_res_typer_db, gbs_res_targets_db)

        // Map genomes to GBS resistance database using SRST2
        srst2_for_res_typing(reads, gbs_res_typer_db, params.gbs_res_min_coverage, params.gbs_res_max_divergence)
        fullgenes = srst2_for_res_typing.out.fullgenes

        // Split sam file for each GBS target sequence
        split_target_RES_seq_from_sam_file(srst2_for_res_typing.out.bam_files, gbs_res_targets_db)

        // Get consensus sequence using freebayes
        freebayes(split_target_RES_seq_from_sam_file.out, split_target_RES_sequences.out)
        consensus = freebayes.out.consensus

    emit:
        fullgenes
        consensus
}

// Resistance mapping with the other resistance databases
workflow OTHER_RES {

    take:
        reads

    main:
        other_res_db = file(params.other_res_db, checkIfExists: true)
        // Map genomes to resistance database using SRST2
        srst2_for_res_typing(reads, other_res_db, params.other_res_min_coverage, params.other_res_max_divergence)
        fullgenes = srst2_for_res_typing.out.fullgenes

    emit:
        fullgenes
}

// MLST pipeline
workflow MLST {

    take:
        reads

    main:
        // Run SRST2 MLST
        srst2_for_mlst(reads, params.mlst_min_coverage)

        // Get new consensus allele and pileup data
        get_mlst_allele_and_pileup(srst2_for_mlst.out.bam_and_srst2_results, params.mlst_min_read_depth)

        // Collect outputs
        new_alleles = get_mlst_allele_and_pileup.out.new_alleles
        pileup = get_mlst_allele_and_pileup.out.pileup
        existing_alleles = get_mlst_allele_and_pileup.out.existing_alleles
        status = get_mlst_allele_and_pileup.out.new_alleles_status
        srst2_results = srst2_for_mlst.out.srst2_results

    emit:
        new_alleles
        pileup
        existing_alleles
        status
        srst2_results
}

// PBP-1A allele typing pipeline
workflow PBP1A {

    take:
        pbp_typer_output

    main:
        // Run
        get_pbp_alleles(pbp_typer_output, 'GBS1A-1', file(params.gbs_blactam_1A_db, checkIfExists: true))

        // Output new PBP alleles to results directory
        get_pbp_alleles.out.new_pbp.subscribe { it ->
            it.copyTo(file("${output_dir}"))
        }

        // Combine existing PBP alleles results in one file
        finalise_pbp_existing_allele_results(get_pbp_alleles.out.existing_pbp, file(params.config, checkIfExists: true))

    emit:
        // Emit existing PBP alleles for collection
        finalise_pbp_existing_allele_results.out
}

// PBP-2B allele typing pipeline
workflow PBP2B {

    take:
        pbp_typer_output

    main:
        // Run
        get_pbp_alleles(pbp_typer_output, 'GBS2B-1', file(params.gbs_blactam_2B_db, checkIfExists: true))

        // Output new PBP alleles to results directory
        get_pbp_alleles.out.new_pbp.subscribe { it ->
            it.copyTo(file("${output_dir}"))
        }

        // Combine existing PBP alleles results in one file
        finalise_pbp_existing_allele_results(get_pbp_alleles.out.existing_pbp, file(params.config, checkIfExists: true))

    emit:
        // Emit existing PBP alleles for collection
        finalise_pbp_existing_allele_results.out
}

// PBP-2X allele typing pipeline
workflow PBP2X {

    take:
        pbp_typer_output

    main:
        // Run
        get_pbp_alleles(pbp_typer_output, 'GBS2X-1', file(params.gbs_blactam_2X_db, checkIfExists: true))

        // Output new PBP alleles to results directory
        get_pbp_alleles.out.new_pbp.subscribe { it ->
            it.copyTo(file("${output_dir}"))
        }

        // Combine existing PBP alleles results in one file
        finalise_pbp_existing_allele_results(get_pbp_alleles.out.existing_pbp, file(params.config, checkIfExists: true))

    emit:
        // Emit existing PBP alleles for collection
        finalise_pbp_existing_allele_results.out
}

// Main Workflow
workflow {

    main:
         // Get path and prefix of Reference Genome BWA Database, generate from assembly if necessary
        GET_REF_GENOME_BWA_DB(params.ref_genome, params.db)

         // Get path to Kraken2 Database, download if necessary
        GET_KRAKEN2_DB(params.kraken2_db_remote, params.db)

        // Get read pairs into Channel raw_read_pairs_ch
        raw_read_pairs_ch = Channel.fromFilePairs("$params.reads/*_{,R}{1,2}{,_001}.{fq,fastq}{,.gz}", checkIfExists: true)

         // Basic input files validation
         // Output into Channel FILE_VALIDATION.out.result
        FILE_VALIDATION(raw_read_pairs_ch)

        // From Channel raw_read_pairs_ch, only output valid reads of samples based on Channel FILE_VALIDATION.out.result
        VALID_READS_ch = FILE_VALIDATION.out.result.join(raw_read_pairs_ch, failOnDuplicate: true, failOnMismatch: true)
                        .filter { it[1] == 'PASS' }
                        .map { it[0, 2..-1] }

        
        // Preprocess valid read pairs
        // Output into Channels PREPROCESS.out.processed_reads & PREPROCESS.out.json
        PREPROCESS(VALID_READS_ch)

        // From Channel PREPROCESS.out.json, provide Read QC status
        // Output into Channels READ_QC.out.bases, READ_QC.out.result, READ_QC.out.report
        READ_QC(PREPROCESS.out.json, params.length_low, params.depth)

        // From Channel PREPROCESS.out.processed_reads, only output reads of samples passed Read QC based on Channel READ_QC.out.result
        READ_QC_PASSED_READS_ch = READ_QC.out.result.join(PREPROCESS.out.processed_reads, failOnDuplicate: true, failOnMismatch: true)
                        .filter { it[1] == 'PASS' }
                        .map { it[0, 2..-1] }

        // From Channel READ_QC_PASSED_READS_ch, assemble the preprocess read pairs
        // Output into Channel ASSEMBLY_ch, and hardlink (default) the assemblies to $params.output directory
        switch (params.assembler) {
            case 'shovill':
                ASSEMBLY_ch = ASSEMBLY_SHOVILL(READ_QC_PASSED_READS_ch, params.min_contig_length, params.assembler_thread)
            break

            case 'unicycler':
                ASSEMBLY_ch = ASSEMBLY_UNICYCLER(READ_QC_PASSED_READS_ch, params.min_contig_length, params.assembler_thread)
            break
        }

        // From Channel ASSEMBLY_ch, assess assembly quality
        // Output into Channel ASSEMBLY_ASSESS.out.report
        ASSEMBLY_ASSESS(ASSEMBLY_ch)

        // From Channel ASSEMBLY_ASSESS.out.report and Channel READ_QC.out.bases, provide Assembly QC status
        // Output into Channels ASSEMBLY_QC.out.result & ASSEMBLY_QC.out.report
        ASSEMBLY_QC(
            ASSEMBLY_ASSESS.out.report
            .join(READ_QC.out.bases, failOnDuplicate: true),
            params.contigs,
            params.length_low,
            params.length_high,
            params.depth
        )

        // From Channel READ_QC_PASSED_READS_ch map reads to reference
        // Output into Channel MAPPING.out.sam
        MAPPING(GET_REF_GENOME_BWA_DB.out.path, GET_REF_GENOME_BWA_DB.out.prefix, READ_QC_PASSED_READS_ch)

        // From Channel MAPPING.out.sam, Convert SAM into sorted BAM and calculate reference coverage
        // Output into Channels SAM_TO_SORTED_BAM.out.sorted_bam and SAM_TO_SORTED_BAM.out.ref_coverage
        SAM_TO_SORTED_BAM(MAPPING.out.sam, params.lite)

        // From Channel SAM_TO_SORTED_BAM.out.sorted_bam calculates non-cluster Het-SNP site count
        // Output into Channel HET_SNP_COUNT.out.result
        SNP_CALL(params.ref_genome, SAM_TO_SORTED_BAM.out.sorted_bam, params.lite)
        HET_SNP_COUNT(SNP_CALL.out.vcf)

        // Merge Channels SAM_TO_SORTED_BAM.out.ref_coverage & HET_SNP_COUNT.out.result to provide Mapping QC Status
        // Output into Channels MAPPING_QC.out.result & MAPPING_QC.out.report
        MAPPING_QC(
        SAM_TO_SORTED_BAM.out.ref_coverage
        .join(HET_SNP_COUNT.out.result, failOnDuplicate: true, failOnMismatch: true),
        params.ref_coverage,
        params.het_snp_site
        )

        // From Channel READ_QC_PASSED_READS_ch assess GBS percentage in reads
        // Output into Channel TAXONOMY.out.report
        TAXONOMY(GET_KRAKEN2_DB.out.path, params.kraken2_memory_mapping, READ_QC_PASSED_READS_ch)

        // From Channel TAXONOMY.out.report, provide taxonomy QC status
         // Output into Channels TAXONOMY_QC.out.result & TAXONOMY_QC.out.report
        TAXONOMY_QC(TAXONOMY.out.report, params.gbs_percentage, params.non_gbs_percentage)

        // Merge Channels FILE_VALIDATION.out.result & READ_QC.out.result & ASSEMBLY_QC.out.result & MAPPING_QC.out.result & TAXONOMY_QC.out.result to provide Overall QC Status
        // Output into Channel OVERALL_QC.out.result & OVERALL_QC.out.report
        OVERALL_QC(
            raw_read_pairs_ch.map{ it[0] }
                .join(FILE_VALIDATION.out.result, failOnDuplicate: true, failOnMismatch: true)
                .join(READ_QC.out.result, failOnDuplicate: true, remainder: true)
                .join(ASSEMBLY_QC.out.result, failOnDuplicate: true, remainder: true)
                .join(MAPPING_QC.out.result, failOnDuplicate: true, remainder: true)
                .join(TAXONOMY_QC.out.result, failOnDuplicate: true, remainder: true)
        )
        // Generate sample reports by merging outputs from all result-generating modules
        GENERATE_SAMPLE_REPORT(
            raw_read_pairs_ch.map{ it[0] }
                .join(READ_QC.out.report, failOnDuplicate: true, remainder: true)
                .join(ASSEMBLY_QC.out.report, failOnDuplicate: true, remainder: true)
                .join(MAPPING_QC.out.report, failOnDuplicate: true, remainder: true)
                .join(TAXONOMY_QC.out.report, failOnDuplicate: true, remainder: true)
                .join(OVERALL_QC.out.report, failOnDuplicate: true, failOnMismatch: true)
            .map { [it[0], it[1..-1].minus(null)] } // Map Sample_ID to index 0 and all reports (with null entries removed) as a list to index 1
        )

        // Generate overall report based on sample reports, ARIBA metadata, resistance to MIC lookup table
        GENERATE_OVERALL_REPORT(GENERATE_SAMPLE_REPORT.out.report.collect())
        
        // From Channel READ_QC_PASSED_READS_ch, only output reads of samples passed overall QC based on Channel OVERALL_QC.out.result
        OVERALL_QC_PASSED_READS_ch = OVERALL_QC.out.result.join(READ_QC_PASSED_READS_ch, failOnDuplicate: true)
                        .filter { it[1] == 'PASS' }
                        .map { it[0, 2..-1] }

        // From Channel ASSEMBLY_ch, only output assemblies of samples passed overall QC based on Channel OVERALL_QC.out.result
        OVERALL_QC_PASSED_ASSEMBLIES_ch = OVERALL_QC.out.result.join(ASSEMBLY_ch, failOnDuplicate: true)
                            .filter { it[1] == 'PASS' }
                            .map { it[0, 2..-1] }

        // Serotyping Process
        serotyping(OVERALL_QC_PASSED_READS_ch, params.serotyper_min_read_depth)

        passed_reads_ch = OVERALL_QC_PASSED_READS_ch.map { sample_id, processed_one, processed_two, processed_unpaired ->  
        [processed_one, processed_two]
        }

        // Resistance Mapping Workflows
        GBS_RES(OVERALL_QC_PASSED_READS_ch)
        OTHER_RES(OVERALL_QC_PASSED_READS_ch)

        // Once GBS or both resistance workflows are complete, trigger resistance typing
        GBS_RES.out.fullgenes
        .join(GBS_RES.out.consensus)
        .join(OTHER_RES.out.fullgenes)
        .set { res_files_ch }

        res_typer(res_files_ch, params.restyper_min_read_depth, file(params.config, checkIfExists: true))

        // Combine serotype and resistance type results for each sample
        sero_res_ch = serotyping.out.join(res_typer.out.res_out)

        finalise_sero_res_results(sero_res_ch, file(params.config, checkIfExists: true))

        // Combine samples and output results files
        finalise_sero_res_results.out.sero_res_incidence
            .collectFile(name: file("${output_dir}/resistance_out/${params.sero_res_incidence_out}"), keepHeader: true)

        finalise_sero_res_results.out.res_alleles_variants
            .collectFile(name: file("${output_dir}/resistance_out/${params.alleles_variants_out}"), keepHeader: true)

        finalise_sero_res_results.out.res_variants
            .collectFile(name: file("${output_dir}/resistance_out/${params.variants_out}"), keepHeader: true)

        res_typer.out.res_accessions
                .collectFile(name: file("${output_dir}/resistance_out/${params.res_accessions_out}"))

        // MLST
        MLST(OVERALL_QC_PASSED_READS_ch)
        
        // Surface Typing Process
        surface_typer(OVERALL_QC_PASSED_READS_ch, file(params.gbs_surface_typer_db, checkIfExists: true),
            params.surfacetyper_min_read_depth, params.surfacetyper_min_coverage,
            params.surfacetyper_max_divergence)

        finalise_surface_typer_results(surface_typer.out, file(params.config, checkIfExists: true))

        // Combine results for surface typing
        finalise_surface_typer_results.out.surface_protein_incidence
            .collectFile(name: file("${output_dir}/surface_protein_out/${params.surface_protein_incidence_out}"), keepHeader: true)
        finalise_surface_typer_results.out.surface_protein_variants
            .collectFile(name: file("${output_dir}/surface_protein_out/${params.surface_protein_variants_out}"), keepHeader: true)

        // PBP Typer
        contig_paths = OVERALL_QC_PASSED_ASSEMBLIES_ch

        get_pbp_genes(contig_paths, file(params.gbs_blactam_db, checkIfExists: true), params.pbp_frac_align_threshold, params.pbp_frac_identity_threshold)

        // Get PBP existing and new alleles
        PBP1A(get_pbp_genes.out)
        PBP2B(get_pbp_genes.out)
        PBP2X(get_pbp_genes.out)

        PBP1A.out
        .concat(PBP2B.out, PBP2X.out)
        .set { PBP_all }

        PBP_all
        .collectFile(name: file("${output_dir}/pbp_alleles/${params.existing_pbp_alleles_out}"), keepHeader: true, sort: true)

        // Combine serotype, resistance, allelic profile, surface typer and GBS resistance variant
        // Get version of pipeline
        get_version()
        version_ch = get_version.out

        // Combine serotype and resistance type results for each sample
        combined_ch = serotyping.out
            .join(res_typer.out.res_out)
            .join(surface_typer.out)
            .join(MLST.out.srst2_results)

        combine_results(combined_ch, file(params.config, checkIfExists: true), version_ch)

        combine_results.out
            .collectFile(name: file("${output_dir}/${params.gbs_typer_report}"), keepHeader: true, sort: true)
}
