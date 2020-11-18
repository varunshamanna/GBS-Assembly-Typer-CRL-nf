process serotyping {

    // Temporarily publish serotyping results - later collate
    publishDir './results'

    input:
    tuple val(pair_id), file(reads)
    file(sero_gene_db)

    output:
    tuple val(pair_id), file("${pair_id}_SeroType_Results.txt")

    """
    # Changed TEMP_SeroType_Results.txt to ${pair_id}_SeroType_Results.txt
    # in GBS_Serotyper.pl so output channel can pick it up
    srst2 --samtools_args '\\-A' --input_pe ${reads[0]} ${reads[1]} --output SERO_${pair_id} --log --save_scores --min_coverage 99.0 --max_divergence 7 --gene_db ${sero_gene_db}
    process_serotyper_results.py --srst2_output SERO_${pair_id} --sero_db ${sero_gene_db} --output ${pair_id}_SeroType_Results.txt --depth-threshold 10
    """
}
