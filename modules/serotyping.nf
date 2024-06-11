process serotyping {

    label 'gbs_typer_container'

    input:
    tuple val(sample_id), path(read1), path(read2), path(unpaired) // ID and paired read files
    val(min_read_depth) // Minimum read depth threshold

    output:
    tuple val(sample_id), file(output_file)

    script:
    output_file="${sample_id}_SeroType_Results.txt"
    sero_gene_db="GBS-SBG.fasta"

    """
    set +e
    # Must create an output file (empty if fails)

    # Get latest version of GBS Serotype Database
    git clone https://github.com/swainechen/GBS-SBG
    mv GBS-SBG/${sero_gene_db} .

    srst2 --samtools_args '\\-A' --input_pe "$read1" "$read2" --output SERO_${sample_id} --log --save_scores --gene_db ${sero_gene_db}
    process_serotyper_results.py --srst2_output SERO_${sample_id} --sero_db ${sero_gene_db} --output ${sample_id}_SeroType_Results.txt --min_read_depth ${min_read_depth}

    touch ${output_file}
    """
}
