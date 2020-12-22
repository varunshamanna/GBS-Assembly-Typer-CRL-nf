process res_typer {

    input:
    tuple val(pair_id), file(arg_fullgenes), file(resfinder_fullgenes), file(RES_fullgenes), file(RES_consensus)

    output:
    tuple val(pair_id), file("${pair_id}_res_incidence.txt"), file("${pair_id}_res_alleles.txt"), file("${pair_id}_res_gbs_variants.txt")

    """
    process_res_typer_results.py --srst2_gbs_fullgenes ${RES_fullgenes} --srst2_gbs_consensus ${RES_consensus} --srst2_other_fullgenes ${arg_fullgenes} ${resfinder_fullgenes} --min_read_depth ${params.restyper_min_read_depth} --output_prefix ${pair_id}
    """
}
