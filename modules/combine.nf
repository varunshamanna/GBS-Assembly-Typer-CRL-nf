process combine_results {

    input:
    tuple val(pair_id), file(sero_results), file(res_results)

    output:
    stdout

    """
    combine_results.py -i ${pair_id} -s ${sero_results} -r ${res_results}
    """
}
