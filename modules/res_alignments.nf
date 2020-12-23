process srst2_for_res_typing {

    input:
    tuple val(pair_id), file(reads)
    val(dbs)
    path(db_dir)
    val(db_name)
    val(min_coverage)
    val(max_divergence)

    publishDir './tmp', mode: 'copy', overwrite: true

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
        db_file=\$(basename \${db_array[i]})
        srst2 --samtools_args '\\-A' --input_pe ${reads[0]} ${reads[1]} --output ${pair_id}_${db_name}_\${db_file} --log --save_scores --min_coverage \${min_cov_array[i]} --max_divergence \${max_div_array[i]} --gene_db ${db_dir}/\${db_file}
    done
    """
}

process split_target_RES_sequences {

    input:
    file(fasta_file)
    file(targets_file)

    output:
    file("CHECK_*")

    """
    get_targets_from_res_db.py -f ${fasta_file} -t ${targets_file} -o CHECK_
    """
}

process split_target_RES_seq_from_sam_file {

    input:
    tuple val(pair_id), file(bam_file)
    file(targets_file)

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
    val(pair_id)
    file(target_bam)
    file(target_bai)
    file(target_ref)
    path(out_dir)

    publishDir './tmp', mode: 'copy', overwrite: true

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
