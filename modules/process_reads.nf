process unzipReads {

    publishDir './results'

    input:
    tuple val(pair_id), file(reads)
    file(contigs)

    output:
    tuple file("${pair_id}_1.fastq"), file("${pair_id}_2.fastq")

    """
    gunzip -c ${reads[0]} > ${pair_id}_1.fastq
    gunzip -c ${reads[1]} > ${pair_id}_2.fastq
    """
}
