process srst2_for_mlst {

    input:
    tuple val(pair_id), file(reads) // ID and paired read files
    file(mlst_alleles) // Path to MLST alleles file
    file(mlst_definitions) // Path to MLST definition file
    val(min_coverage) // String of minimum coverage parameter(s) for SRST2

    //publishDir "./${tmp_dir}/${pair_id}", mode: 'move', overwrite: true, pattern: "${pair_id}_${db_name}_*__fullgenes__*__results.txt"

    output:
    tuple val(pair_id), file("${pair_id}*.bam"), file("${pair_id}__mlst__*__results.txt")

    """
    srst2 --samtools_args '\\-A' --mlst_delimiter '_' --input_pe ${reads[0]} ${reads[1]} --output ${pair_id} --save_scores --mlst_db ${mlst_alleles} --mlst_definitions ${mlst_definitions} --min_coverage ${min_coverage}
    """
}

process get_mlst_allele_and_pileup {

    input:
    tuple val(pair_id), file(bam_file), file(results_file)
    val(min_read_depth)
    file(mlst_alleles)

    output:
    file("${pair_id}_new_mlst_alleles.txt")

    """
    # Get alleles from mismatches in SRST2 MLST results file
    samtools index ${bam_file}
    get_alleles_from_srst2_mlst.py --mlst_results_file ${results_file} --min_read_depth ${min_read_depth} --output_file ${pair_id}_mlst_alleles.txt
    sed -n '1p' ${pair_id}_mlst_alleles.txt > ${pair_id}_mlst_alleles_header.txt
    sed '1d' ${pair_id}_mlst_alleles.txt > ${pair_id}_tmp.txt; mv ${pair_id}_tmp.txt ${pair_id}_mlst_alleles.txt
    num_alleles=\$(cat ${pair_id}_mlst_alleles.txt | wc -l)

    # Get consensus allele and variant pileup for each allele
    if [[ \${num_alleles} -gt 0 ]]; then
        for ((i=0;i<\${num_alleles};i++)); # Skip first line of alleles file
        do
            target_allele=\$(sed -n \"\${i+1}p\" ${pair_id}_mlst_alleles.txt)

            mlst_bam=\"\${target_allele}_${bam_file}\"
            mlst_vcf=\"\$(basename \${mlst_bam} .bam).vcf\"
            samtools view -b ${bam_file} \${target_allele} > \${mlst_bam}
            samtools index \${mlst_bam} \$(basename \${mlst_bam} .bam).bai
            echo \${target_allele} > \${target_allele}.txt
            get_targets_from_db.py -f ${mlst_alleles} -t \${target_allele}.txt -o CHECK_MLST_
            freebayes -q 20 -p 1 -f CHECK_MLST_\${target_allele}_ref.fna \${mlst_bam} -v \${mlst_vcf}
            bgzip \${mlst_vcf}
            tabix -p vcf \${mlst_vcf}.gz

            echo 'New MLST Allele Consensus:' >> ${pair_id}_new_mlst_alleles.txt
            cat CHECK_MLST_\${target_allele}_ref.fna | vcf-consensus \${mlst_vcf}.gz >> ${pair_id}_new_mlst_alleles.txt
            echo '\nNew MLST Allele Pileup:' >> ${pair_id}_new_mlst_alleles.txt
            samtools mpileup -f ${mlst_alleles} ${bam_file} -r \${target_allele} >> ${pair_id}_new_mlst_alleles.txt
            echo '\n' >> ${pair_id}_new_mlst_alleles.txt
        done
    else
        cat ${pair_id}_mlst_alleles_header.txt > ${pair_id}_new_mlst_alleles.txt
    fi

    """

}
