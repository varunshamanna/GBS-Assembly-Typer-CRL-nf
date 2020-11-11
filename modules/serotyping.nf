process serotyping {

    // Temporarily publish serotyping results - later collate
    publishDir './results'

    input:
    tuple val(pair_id), file(reads)
    file(sero_gene_db)

    output:
    file("${pair_id}_SeroType_Results.txt")

    """
    # Changed TEMP_SeroType_Results.txt to ${pair_id}_SeroType_Results.txt
    # in GBS_Serotyper.pl so output channel can pick it up
    GBS_Serotyper.pl -1 ${reads[0]} -2 ${reads[1]} -r ${sero_gene_db} -n ${pair_id}
    """
}
