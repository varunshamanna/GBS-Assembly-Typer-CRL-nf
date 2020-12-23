process res_typer {

    input:
    val(pair_id)
    path(res_dir)

    output:
    tuple val(pair_id), file("${pair_id}_res_incidence.txt"), file("${pair_id}_res_alleles.txt"), file("${pair_id}_res_gbs_variants.txt")

    """
    process_res_typer_results.py --srst2_gbs_fullgenes ${res_dir}/${pair_id}_RES_*__fullgenes__*__results.txt --srst2_gbs_consensus ${res_dir}/${pair_id}_consensus_seq.fna --srst2_other_fullgenes ${res_dir}/${pair_id}_OTHER_RES_*__fullgenes__*__results.txt --output_prefix ${pair_id}
    """
}
