process res_typer {

    // Temporarily publish results - later collate
    publishDir './results'

    input:
    tuple val(pair_id), file(reads)
    file(res_typer_gene_db)
    file(argannot_db)
    file(resfinder_db)

    output:
    file("${pair_id}_Res_Results.txt")

    """
    GBS_Res_Typer.pl -1 ${reads[0]} -2 ${reads[1]} -a ${argannot_db} -b ${resfinder_db} -r ${res_typer_gene_db} -n ${pair_id}
    """
}