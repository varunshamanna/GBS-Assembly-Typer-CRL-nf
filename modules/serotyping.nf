process serotyping {

    // Temporarily publish serotyping results - later collate
    publishDir './results'

    input:
    tuple val(pair_id), file(reads)
    file(sero_gene_db)

    output:
    file("${pair_id}_SeroType_Results.txt")

    """
    # Run SRST2 TODO: Make min_coverage and max_divergence options
    srst2 --samtools_args '\\-A' --input_pe ${reads[0]} ${reads[1]} --output SERO_${pair_id} --log --save_scores --min_coverage 99.0 --max_divergence 7 --gene_db ${sero_gene_db}

    # Process fullgenes results file: TODO: catch non-existent files
    sero_gene_db_name=\$(echo ${sero_gene_db} | sed 's/\.[^.]*//')
    process_serotyper_results.py --input SERO_${pair_id}__fullgenes__\${sero_gene_db_name}__results.txt --output ${pair_id}_SeroType_Results.txt --depth-threshold 10
    """
}
