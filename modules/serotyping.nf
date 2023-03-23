process serotyping {

    input:
    tuple val(pair_id), file(reads) // ID and paired read files
    val(min_read_depth) // Minimum read depth threshold

    output:
    tuple val(pair_id), file(output_file)

    script:
    output_file="${pair_id}_SeroType_Results.txt"
    sero_gene_db="GBS-SBG.fasta"

    """
    set +e
    # Must create an output file (empty if fails)

    # Get latest version of GBS Serotype Database
    git clone https://github.com/swainechen/GBS-SBG
    mv GBS-SBG/${sero_gene_db} .

    srst2 --samtools_args '\\-A' --input_pe ${reads[0]} ${reads[1]} --output SERO_${pair_id} --log --save_scores --gene_db ${sero_gene_db}
    process_serotyper_results.py --srst2_output SERO_${pair_id} --sero_db ${sero_gene_db} --output ${pair_id}_SeroType_Results.txt --min_read_depth ${min_read_depth}

    touch ${output_file}

    # Clean directory
    mkdir output
    mv ${output_file} output
    find . -maxdepth 1 -type f -delete
    unlink ${reads[0]}
    unlink ${reads[1]}
    mv output/${output_file} .
    rm -d output
    """
}
