process PHYLOGENY {
    tag "Phylogeny_${organelle}"
    label 'process_medium'

    publishDir "${params.outdir}/phylogeny/${organelle}", mode: 'copy'

    input:
    tuple val(organelle), val(sample_ids), path(cds_files)

    output:
    path "*.treefile", optional: true, emit: tree_file
    path "*.fasta", emit: alignment_file
    path "*.iqtree", optional: true, emit: report

    script:
    def samples_arg = sample_ids.join(' ')
    def cds_arg = cds_files.join(' ')
    
    def outgroup_val = (organelle == 'plastome') ? params.outgroup_plastome : params.outgroup_mitogenome
    def outgroup_param = outgroup_val ? "--outgroup \"${outgroup_val}\"" : ""

    """
    python3 ${workflow.projectDir}/bin/build_phylogeny.py \
        --organelle "${organelle}" \
        --samples ${samples_arg} \
        --cds_files ${cds_arg} \
        ${outgroup_param}
    """
}
