process KAKS_ANALYSIS {
    tag "KaKs_${organelle}"
    label 'process_medium'

    publishDir "${params.outdir}/kaks/${organelle}", mode: 'copy'

    when:
    (organelle == 'plastome' && params.ks_ka_reference_plastome != null) ||
    (organelle == 'mitogenome' && params.ks_ka_reference_mitogenome != null)

    input:
    tuple val(organelle), val(sample_ids), path(cds_files)

    output:
    path "kaks_summary.csv", emit: summary_csv
    path "p*.pdf", optional: true, emit: plots

    script:
    def samples_str = (sample_ids instanceof Collection) ? sample_ids.join(',') : sample_ids
    def cds_str     = (cds_files instanceof Collection) ? cds_files.join(',') : cds_files
    
    def ref_org = (organelle == 'plastome') ? params.ks_ka_reference_plastome : params.ks_ka_reference_mitogenome
    
    """
    Rscript ${workflow.projectDir}/bin/calc_kaks.R \
        --organelle "${organelle}" \
        --samples "${samples_str}" \
        --cds_files "${cds_str}" \
        --reference "${ref_org}"

    Rscript ${workflow.projectDir}/bin/plot_kaks.R \
        --input kaks_summary.csv || true
    """
}
