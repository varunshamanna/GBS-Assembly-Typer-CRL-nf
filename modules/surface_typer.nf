process surface_typer {

    input:
    tuple val(sample_id), path(read1), path(read2), path(unpaired)
    file(surface_protein_db)
    val(min_read_depth) // Minimum read depth threshold
    val(min_coverage) // Minimum coverage threshold
    val(max_divergence) // Maximum allowed divergence threshold

    output:
    tuple val(sample_id), file(inc_output_file), file(variants_output_file)

    script:
    inc_output_file="${sample_id}_surface_protein_incidence_sample.txt"
    variants_output_file="${sample_id}_surface_protein_variants_sample.txt"
    """
    set +e

    srst2 --samtools_args '\\-A' --input_pe ${read1} ${read2} --output ${sample_id}_SURFACE --log --save_scores --min_coverage ${min_coverage} --max_divergence ${max_divergence} --gene_db ${surface_protein_db}
    process_surface_typer_results.py --srst2_gbs_fullgenes ${sample_id}_SURFACE --surface_db ${surface_protein_db} --output_prefix ${sample_id} --min_read_depth ${min_read_depth}

    touch ${inc_output_file}
    touch ${variants_output_file}

    # Clean directory
    mkdir output
    mv ${inc_output_file} output
    mv ${variants_output_file} output
    find . -maxdepth 1 -type f -delete
    unlink ${surface_protein_db}
    mv output/${inc_output_file} .
    mv output/${variants_output_file} .
    rm -d output
    """
}
