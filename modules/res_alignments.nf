process srst2_for_res_typing {

    input:
    tuple val(pair_id), file(reads) // ID and paired read files
    val(dbs) // String of resistance database file(s)
    path(db_dir) // Path to database file(s)
    path(tmp_dir) // Path to temporary directory
    val(db_name) // Type of resistance database i.e. GBS or OTHER
    val(min_coverage) // String of minimum coverage parameter(s) for SRST2
    val(max_divergence) // String of maximum coverage parameter(s) for SRST2

    publishDir "./${tmp_dir}/${pair_id}", mode: 'move', overwrite: true, pattern: "${pair_id}_${db_name}_*__fullgenes__*__results.txt"

    output:
    tuple val(pair_id), file("${pair_id}*.bam"), emit: bam_files
    val(pair_id), emit: id
    file("${pair_id}_${db_name}_*__fullgenes__*__results.txt")

    """
    db_list='${dbs}'
    db_array=(\$db_list)

    min_cov_list='${min_coverage}'
    min_cov_array=(\$min_cov_list)

    max_div_list='${max_divergence}'
    max_div_array=(\$max_div_list)

    for ((i=0;i<\${#db_array[@]};i++));
    do
        mkdir -p ${tmp_dir}/${pair_id}
        cp ${db_dir}/\${db_array[i]} ${tmp_dir}/${pair_id}/
        db_file=\$(basename \${db_array[i]})
        srst2 --samtools_args '\\-A' --input_pe ${reads[0]} ${reads[1]} --output ${pair_id}_${db_name}_\${db_file} --log --save_scores --min_coverage \${min_cov_array[i]} --max_divergence \${max_div_array[i]} --gene_db ${tmp_dir}/${pair_id}/\${db_file}
        rm ${tmp_dir}/${pair_id}/\${db_file}*
    done
    """
}

process split_target_RES_sequences {

    input:
    file(fasta_file) // FASTA file of GBS target sequences
    file(targets_file) // Text file of GBS targets of interest

    output:
    file("CHECK_*")

    """
    get_targets_from_res_db.py -f ${fasta_file} -t ${targets_file} -o CHECK_
    """
}

process split_target_RES_seq_from_sam_file {

    input:
    tuple val(pair_id), file(bam_file) // ID and corresponding BAM file from mapping
    file(targets_file) // Text file of GBS targets of interest

    output:
    val(pair_id)
    file("CHECK_*_${pair_id}*.bam")
    file("CHECK_*_${pair_id}*.bai")

    """
    samtools view -h ${bam_file} > \$(basename ${bam_file} .bam).sam
    get_targets_from_samfile.py -s \$(basename ${bam_file} .bam).sam -t ${targets_file} -i ${pair_id} -o CHECK_
    for check_sam_file in CHECK_*_${pair_id}*.sam; do
        samtools view -bS \${check_sam_file} > \$(basename \${check_sam_file} .sam).bam
        samtools index \$(basename \${check_sam_file} .sam).bam \$(basename \${check_sam_file} .sam).bai
    done
    """
}

process freebayes {

    input:
    val(pair_id) // ID
    file(target_bam) // BAM file from a mapped target sequence of interest
    file(target_bai) // Corresponding BAM index file
    file(target_ref) // FASTA file of target sequence
    path(tmp_dir)

    publishDir "./${tmp_dir}/${pair_id}", mode: 'move', overwrite: true

    output:
    val(pair_id), emit: id
    file("${pair_id}_consensus_seq.fna")

    """
    for check_bam_file in CHECK_*_${pair_id}*.bam; do
        target=\$(echo \${check_bam_file} | sed 's/CHECK_//g' | sed 's/_${pair_id}.*//g')
        freebayes -q 20 -p 1 -f CHECK_\${target}_ref.fna \${check_bam_file} -v CHECK_\${target}_${pair_id}_seq.vcf
        bgzip CHECK_\${target}_${pair_id}_seq.vcf
        tabix -p vcf CHECK_\${target}_${pair_id}_seq.vcf.gz
        cat CHECK_\${target}_ref.fna | vcf-consensus CHECK_\${target}_${pair_id}_seq.vcf.gz >> ${pair_id}_consensus_seq.fna
    done
    """
}
