process GENERATE_SAMPLE_REPORT {
    label 'bash_container'
    label 'farm_low'

    tag "$sample_id"

    publishDir "${params.output}/sample_report", mode: "copy"

    input:
    tuple val(sample_id), path("${sample_id}_process_report_?.csv")

    output:
    path sample_report, emit: report

    script:
    sample_report="${sample_id}_report.csv"
    """
    SAMPLE_ID="$sample_id"
    SAMPLE_REPORT="$sample_report"

    source generate_sample_report.sh
    """
}

process GENERATE_OVERALL_REPORT {
    label 'csvtk_container'
    label 'farm_low'

    publishDir "${params.output}", mode: "copy"

    input:
    path '*'

    output:
    path "$overall_report", emit: report

    script:
    input_pattern='*_report.csv'
    overall_report='GBS_QC_results.tsv'
    """
    csvtk concat --out-delimiter \$\'\t\' *.csv > "$overall_report"
    """
}
