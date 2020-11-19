process res_typer {

    // Temporarily publish results - later collate
    publishDir './results'

    input:
    tuple val(pair_id), file(reads)
    file(res_typer_gene_db)
    file(argannot_db)
    file(resfinder_db)

    output:
    tuple val(pair_id), file("${pair_id}_BIN_Res_Results.txt")

    """

    GBS_Res_Typer.pl -1 ${reads[0]} -2 ${reads[1]} -a ${argannot_db} -b ${resfinder_db} -r ${res_typer_gene_db} -n ${pair_id}

    #srst2 --samtools_args '\\-A' --input_pe ${reads[0]} ${reads[1]} --output RES_${pair_id} --log --save_scores \
    #    --min_coverage 99.9 --max_divergence 5 --gene_db ${res_typer_gene_db}
    #srst2 --samtools_args '\\-A' --input_pe ${reads[0]} ${reads[1]} --output ARG_${pair_id} --log --save_scores \
    #    --min_coverage 70 --max_divergence 30 --gene_db ${argannot_db}
    #srst2 --samtools_args '\\-A' --input_pe ${reads[0]} ${reads[1]} --output RESFI_${pair_id} --log --save_scores \
    #    --min_coverage 70 --max_divergence 30 --gene_db ${resfinder_db}

    #process_res_typer_results.py \
    #    --srst2_gbs RES_${pair_id} \
    #    --srst2_argannot ARG_${pair_id} \
    #    --srst2_resfinder RESFI_${pair_id} \
    #    --output ${pair_id}_Res_Results.txt \
    #    --output_bin ${pair_id}_BIN_Res_Results.txt

    """

}
