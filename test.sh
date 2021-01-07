#!/usr/bin/env bash

in_dir=test_data/input_data
out_dir=test_data/output_data

# Run pipeline on test input data
nextflow run main.nf --reads "$in_dir/test_{1,2}.fastq.gz" --results_dir "$out_dir" --output 'test'

function file_diff {
    local test=${1}
    local reference=${2}

    diff_status=$(diff -s ${test} ${reference})
    # If files different then warning
    if [[ $diff_status != "Files ${test} and ${reference} are identical" ]]; then
        echo "${test} is missing or its contents not expected."
        return 1
    fi
}

error_status=0
# Check for test_drug_cat_alleles.txt output
file_diff ${out_dir}/test_drug_cat_alleles_variants.txt ${out_dir}/reference_drug_cat_alleles_variants.txt
out=$?
error_status=$(($error_status | $out))

# Check for test_gbs_res_variants.txt output
file_diff ${out_dir}/test_gbs_res_variants.txt ${out_dir}/reference_gbs_res_variants.txt
out=$?
error_status=$(($error_status | $out))

# Check for test_serotype_res_incidence.txt output
file_diff ${out_dir}/test_serotype_res_incidence.txt ${out_dir}/reference_serotype_res_incidence.txt
out=$?
error_status=$(($error_status | $out))

# Error if any output files missing or not expected
if [ $error_status = 1 ]; then
    echo ""
    echo "Test failed. Outputs listed above are missing or their contents not expected."
    exit 1
else
    echo "Test passed. All outputs expected."
fi
