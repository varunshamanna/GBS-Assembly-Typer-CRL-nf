process combine_results {

    input:
    tuple val(pair_id), file(sero_results), file(res_incidence), file(res_alleles), file(res_variants)

    output:
    path("${pair_id}_sero_res_incidence.txt"), emit: sero_res_incidence
    path(res_alleles), emit: res_alleles
    path(res_variants), emit: res_variants

    """
    combine_results.py -i ${pair_id} -s ${sero_results} -r ${res_incidence} -a ${res_alleles} -v ${res_variants} -o ${pair_id}
    """
}
