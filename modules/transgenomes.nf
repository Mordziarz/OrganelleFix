process TRANSGENOMES {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/transgenomes_results/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(mt_fasta), path(mt_bed), path(pt_fasta), path(pt_bed)

    output:
    path "transfers_plot1.pdf", emit: plot1
    path "transfers_plot2.pdf", emit: plot2
    path "transfers_table.csv", emit: table

    script:
    """
    Rscript ${projectDir}/bin/run_transgenomes.R \
        --mt-fasta ${mt_fasta} \
        --mt-bed ${mt_bed} \
        --pt-fasta ${pt_fasta} \
        --pt-bed ${pt_bed} \
        --sample-id ${sample_id}
    """
}
