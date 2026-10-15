process PREPARE_GENOME {
    tag "$meta.sample_id ($meta.organelle)"
    label 'process_low'

    input:
    tuple val(meta), path(fasta), path(gbk), path(bed)

    output:
    tuple val(meta), path("${meta.sample_id}_${meta.organelle}_genome.fasta"), path("${meta.sample_id}_${meta.organelle}_CDS.fasta"), path("${meta.sample_id}_${meta.organelle}_annotations.bed"), emit: genome_data

    script:
    """
    python3 ${projectDir}/bin/prepare_genome.py \
        --sample-id "${meta.sample_id}" \
        --fasta "${fasta}" \
        --gbk "${gbk}" \
        --bed "${bed}"

    mv genome.fasta ${meta.sample_id}_${meta.organelle}_genome.fasta
    mv CDS.fasta ${meta.sample_id}_${meta.organelle}_CDS.fasta
    mv annotations.bed ${meta.sample_id}_${meta.organelle}_annotations.bed
    """
}
