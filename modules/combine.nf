process combine_results {

    input:
    // ID, serotyping results, resistance incidence, resistance alleles and resistance variants results
    tuple val(pair_id), file(sero_results), file(res_incidence), file(res_alleles), file(res_variants)

    output:
    path("${pair_id}_sero_res_incidence.txt"), emit: sero_res_incidence
    path("${pair_id}_id_alleles_variants.txt"), emit: res_alleles_variants
    path("${pair_id}_id_variants.txt"), emit: res_variants

    """
    combine_results.py -i ${pair_id} -s ${sero_results} -r ${res_incidence} -a ${res_alleles} -v ${res_variants} -o ${pair_id}
    """
}
