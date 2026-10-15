process RSCU_ANALYSIS {
    label 'process_medium'
    publishDir "${params.outdir}/rscu_results", mode: 'copy'

    input:
    tuple val(sample_ids), path(cds_files)

    output:
    path "combined_cds.fasta", emit: combined_fasta
    path "rscu_samples_table.csv", emit: table
    path "rscu_matrix_output.csv", emit: matrix_global
    path "rscu_matrix_output_plastome.csv", emit: matrix_plastome
    path "rscu_matrix_output_mitogenome.csv", emit: matrix_mitogenome
    path "rscu_heatmap_plastome.pdf", emit: heatmap_plastome
    path "rscu_heatmap_mitogenome.pdf", emit: heatmap_mitogenome
    path "rscu_histogram_*.pdf", emit: histograms

    script:
    def ids   = (sample_ids instanceof Collection) ? sample_ids : [sample_ids]
    def files = (cds_files instanceof Collection)   ? cds_files   : [cds_files]

    def args_list = [ids, files].transpose().collect { id, file -> "${id}=${file}" }.join(' ')
    """
    Rscript ${projectDir}/bin/run_rscu.R ${args_list}
    """
}
