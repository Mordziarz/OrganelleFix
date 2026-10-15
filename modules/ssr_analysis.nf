process SSR_ANALYSIS {
    tag "${meta.sample_id}_${meta.organelle}"
    label 'process_medium'

    publishDir "${params.outdir}/ssr/per_sample", mode: 'copy'

    input:
    tuple val(meta), path(genome), path(cds), path(bed)

    output:
    path "*_ssr_results.csv", emit: sample_csv

    script:
    """
    Rscript ${workflow.projectDir}/bin/run_ssr.R \
        --fasta ${genome} \
        --gff ${bed} \
        --sample-id ${meta.sample_id}_${meta.organelle}
    """
}

process SSR_SUMMARY {
    label 'process_medium'

    publishDir "${params.outdir}/ssr/summary", mode: 'copy'

    input:
    path csv_files

    output:
    path "all_ssr_summary.csv", emit: summary_csv
    path "ssr_*.pdf", emit: summary_plots

    script:
    """
    Rscript ${workflow.projectDir}/bin/summarize_ssr.R
    """
}
