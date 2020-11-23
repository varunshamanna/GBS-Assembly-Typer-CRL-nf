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
    """
}

process srst2_for_res_typing {

    input:
    tuple val(pair_id), file(reads)
    file(db)
    val(db_name)
    val(min_coverage)
    val(max_divergence)

    output:
    val(pair_id)
    file("${db_name}_${pair_id}*.bam")

    """
    srst2 --samtools_args '\\-A' --input_pe ${reads[0]} ${reads[1]} --output ${db_name}_${pair_id} --log --save_scores --min_coverage ${min_coverage} --max_divergence ${max_divergence} --gene_db ${db}
    """
}

process extract_target_seq {

    input:
    val(pair_id)
    file(bam_file)

    """
    samtools view -h ${bam_file} > ${bam_file}.sam
    #get_target_from_samfile -s ${bam_file}.sam -t <target> -o CHECK_target_seq.sam
    #samtools view -bS CHECK_target_seq.sam > CHECK_target_seq.bam
    #samtools index CHECK_target_seq.bam CHECK_target_seq.bai
    """
}
