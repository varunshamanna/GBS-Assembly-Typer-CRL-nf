process srst2_for_mlst {

    publishDir "${params.output}/mlst_out/${sample_id}", mode: "copy", pattern: '*__mlst__*__results.txt'

    input:
    tuple val(sample_id), path(read1), path(read2), path(unpaired) // ID and paired read files
    val(min_coverage) // String of minimum coverage parameter(s) for SRST2

    //publishDir "./${tmp_dir}/${sample_id}", mode: 'move', overwrite: true, pattern: "${sample_id}_${db_name}_*__fullgenes__*__results.txt"

    output:
    tuple val(sample_id), file("${sample_id}*.bam"), file("${sample_id}__mlst__${mlst_name}__results.txt"), file(mlst_db), emit: bam_and_srst2_results
    tuple val(sample_id), file("${sample_id}__mlst__${mlst_name}__results.txt"), emit: srst2_results

    script:
    mlst_db="Streptococcus_agalactiae.fasta"
    mlst_name="Streptococcus_agalactiae"

    """
    set +e

    getmlst.py --species 'Streptococcus agalactiae'
    srst2 --samtools_args '\\-A' --input_pe ${read1} ${read2} --output ${sample_id} --save_scores --mlst_db ${mlst_db} --mlst_definitions profiles_csv --mlst_delimiter '_' --min_coverage ${min_coverage}

    touch ${sample_id}__mlst__${mlst_name}__results.txt

    # Clean directory
    mkdir output
    mv ${sample_id}*.bam output
    mv ${sample_id}__mlst__${mlst_name}__results.txt output
    mv ${mlst_db} output
    find . -maxdepth 1 -type f -delete
    unlink ${read1}
    unlink ${read2}
    mv output/${sample_id}*.bam .
    mv output/${sample_id}__mlst__${mlst_name}__results.txt .
    mv output/${mlst_db} .
    rm -d output
    """
}

process get_mlst_allele_and_pileup {

    tag "$sample_id"

    publishDir "${params.output}/mlst_out/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), file(bam_file), file(results_file), file(mlst_alleles)
    val(min_read_depth)

    output:
    path(output_new_mlst_alleles_fasta), emit: new_alleles, optional: true
    path(output_new_mlst_pileup), emit: pileup, optional: true
    path(output_existing_mlst_alleles), emit: existing_alleles, optional: true
    path(output_new_mlst_alleles_log), emit: new_alleles_status

    script:
    output_new_mlst_alleles_fasta="${sample_id}_new_mlst_alleles.fasta"
    output_new_mlst_pileup="${sample_id}_new_mlst_pileup.txt"
    output_existing_mlst_alleles="${sample_id}_existing_mlst_alleles.txt"
    output_new_mlst_alleles_log="${sample_id}_new_mlst_alleles.log"

    """

    # Get alleles from mismatches in SRST2 MLST results file
    samtools index ${bam_file}
    get_alleles_from_srst2_mlst.py --mlst_results_file ${results_file} --min_read_depth ${min_read_depth} --output_prefix ${sample_id}
    if [[ -f ${sample_id}_new_mlst_alleles.txt ]]
    then
        num_alleles=\$(cat ${sample_id}_new_mlst_alleles.txt | wc -l)
    else
        num_alleles=0
    fi

    # Get consensus allele and variant pileup for each allele
    if [[ \${num_alleles} -gt 1 ]]
    then
        echo "${sample_id}: New MLST alleles found." > tmp.log
        for ((i=2;i<=\${num_alleles};i++)); # Skip first line of alleles file
        do
            target_allele=\$(sed -n \"\${i}p\" ${sample_id}_new_mlst_alleles.txt)

            mlst_bam=\"\${target_allele}_${bam_file}\"
            mlst_vcf=\"\$(basename \${mlst_bam} .bam).vcf\"
            samtools view -b ${bam_file} \${target_allele} > \${mlst_bam}
            samtools index \${mlst_bam} \$(basename \${mlst_bam} .bam).bai
            echo \${target_allele} > \${target_allele}.txt
            get_targets_from_db.py -f ${mlst_alleles} -t \${target_allele}.txt -o CHECK_MLST_
            freebayes -q 20 -p 1 -f CHECK_MLST_\${target_allele}_ref.fna \${mlst_bam} -v \${mlst_vcf}
            bgzip \${mlst_vcf}
            tabix -p vcf \${mlst_vcf}.gz

            cat CHECK_MLST_\${target_allele}_ref.fna | vcf-consensus \${mlst_vcf}.gz >> tmp.fasta
            echo '\nNew MLST Allele Pileup:' >> tmp_pileup.txt
            samtools mpileup -f ${mlst_alleles} ${bam_file} -r \${target_allele} >> tmp_pileup.txt
        done
    elif [[ \${num_alleles} -eq 1 ]]
    then
        cat ${sample_id}_new_mlst_alleles.txt > tmp.log
    else
        echo "${sample_id}: No new MLST alleles found." > tmp.log
    fi

    if [ -f tmp.fasta ]
    then
        mv tmp.fasta ${output_new_mlst_alleles_fasta}
    fi

    if [ -f tmp_pileup.txt ]
    then
        mv tmp_pileup.txt ${output_new_mlst_pileup}
    fi

    mv tmp.log ${output_new_mlst_alleles_log}

    # Clean
    unlink ${bam_file}
    unlink ${results_file}
    unlink ${mlst_alleles}
    """

}
