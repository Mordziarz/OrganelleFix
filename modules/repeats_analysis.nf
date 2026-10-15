process REPEATS_ANALYSIS {
    tag "${meta.sample_id}_${meta.organelle}"
    label 'process_medium'

    publishDir "${params.outdir}/repeats/per_sample", mode: 'copy'

    input:
    tuple val(meta), path(genome), path(cds), path(bed)

    output:
    path "*_repeats_results.csv", emit: sample_csv

    script:
    """
    Rscript ${workflow.projectDir}/bin/run_repeats.R \
        --fasta ${genome} \
        --gff ${bed} \
        --sample-id ${meta.sample_id}_${meta.organelle}
    """
}

process REPEATS_SUMMARY {
    label 'process_medium'

    publishDir "${params.outdir}/repeats/summary", mode: 'copy'

    input:
    path csv_files

    output:
    path "all_repeats_summary.csv", emit: summary_csv
    path "repeat*.pdf", emit: summary_plots

    script:
    """
    Rscript ${workflow.projectDir}/bin/summarize_repeats.R
    """
}
