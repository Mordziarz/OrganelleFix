process RNA_EDITING {
    tag "RNA_Edit_${meta.sample_id}_${meta.organelle}"
    label 'process_low'

    publishDir "${params.outdir}/rna_editing/${meta.organelle}", mode: 'copy'

    input:
    tuple val(meta), path(genome), path(cds), path(bed)

    output:
    path "*_rna_editing_predictions.csv", emit: editing_csv

    script:
    """
    python3 ${workflow.projectDir}/bin/predict_rna_editing.py \
        --sample "${meta.sample_id}" \
        --organelle "${meta.organelle}" \
        --cds "${cds}"
    """
}
