process srst2_for_res_typing {

    input:
    tuple val(pair_id), file(reads)
    file(db)
    val(db_name)
    val(min_coverage)
    val(max_divergence)

    output:
    tuple val(pair_id), file("${db_name}_${pair_id}*.bam"), emit: bam_files
    tuple val(pair_id), file("*${pair_id}__fullgenes__*__results.txt"), emit: genes_files

    """
    srst2 --samtools_args '\\-A' --input_pe ${reads[0]} ${reads[1]} --output ${db_name}_${pair_id} --log --save_scores --min_coverage ${min_coverage} --max_divergence ${max_divergence} --gene_db ${db}
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

    output:
    tuple val(pair_id), file("${pair_id}_consensus_seq.fna")

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
