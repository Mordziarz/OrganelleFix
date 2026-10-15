nextflow.enable.dsl = 2

params.ks_ka_reference_plastome  = null
params.ks_ka_reference_mitogenome = null
params.outgroup_plastome  = null
params.outgroup_mitogenome = null

include { PREPARE_GENOME } from './modules/prepare_genome'
include { RSCU_ANALYSIS }  from './modules/rscu_analysis'
include { TRANSGENOMES }   from './modules/transgenomes'
include { SSR_ANALYSIS; SSR_SUMMARY } from './modules/ssr_analysis'
include { REPEATS_ANALYSIS; REPEATS_SUMMARY } from './modules/repeats_analysis'
include { KAKS_ANALYSIS } from './modules/kaks_analysis'
include { RNA_EDITING }   from './modules/rna_editing'
include { PHYLOGENY }     from './modules/phylogeny'

workflow {
    if (!params.samples_table) {
        error "Missing required parameter: --samples_table."
    }

    input_ch = Channel
        .fromPath(params.samples_table, checkIfExists: true)
        | splitCsv(header: true)
        | map { row -> 
            if (!row.sample_id || !row.organelle || !row.comparison_group) {
                error "Error in CSV file! Check if headers include: sample_id, organelle, comparison_group, fasta_path, gbk_path, bed_path"
            }
          
            def meta = [
                sample_id: row.sample_id,
                organelle: row.organelle,
                group: row.comparison_group
            ]
          
            def fasta = (row.fasta_path && row.fasta_path != "") ? file(row.fasta_path, checkIfExists: true) : file("NO_FASTA_${row.sample_id}")
            def gbk   = (row.gbk_path   && row.gbk_path   != "") ? file(row.gbk_path,   checkIfExists: true) : file("NO_GBK_${row.sample_id}")
            def bed   = (row.bed_path   && row.bed_path   != "") ? file(row.bed_path,   checkIfExists: true) : file("NO_BED_${row.sample_id}")
                    
            return [ meta, fasta, gbk, bed ]
        }

    PREPARE_GENOME(input_ch)

    ch_all_cds = PREPARE_GENOME.out.genome_data
        .map { meta, genome, cds, bed -> 
            def unique_id = "${meta.sample_id}_${meta.organelle}"
            return tuple(unique_id, cds) 
        }
        .toList()
        .map { list -> 
            def sample_ids = list.collect { it[0] }
            def cds_files  = list.collect { it[1] }
            return [ sample_ids, cds_files ]
        }

    RSCU_ANALYSIS(ch_all_cds)

    ch_plastome_data = PREPARE_GENOME.out.genome_data
        .filter { meta, genome, cds, bed -> meta.organelle == 'plastome' }

    ch_plastome_kaks = ch_plastome_data
        .map { meta, genome, cds, bed -> tuple(meta.sample_id, cds) }
        .toList()
        .map { list -> [ 'plastome', list.collect { it[0] }, list.collect { it[1] } ] }

    ch_mitogenome_data = PREPARE_GENOME.out.genome_data
        .filter { meta, genome, cds, bed -> meta.organelle == 'mitogenome' }

    ch_mitogenome_kaks = ch_mitogenome_data
        .map { meta, genome, cds, bed -> tuple(meta.sample_id, cds) }
        .toList()
        .map { list -> [ 'mitogenome', list.collect { it[0] }, list.collect { it[1] } ] }

    ch_transgenomes_input = PREPARE_GENOME.out.genome_data
        .map { meta, genome, cds, bed -> tuple(meta.sample_id, meta.organelle, genome, bed) }
        .groupTuple(by: 0)
        .filter { sample_id, organelles, genomes, beds -> 
            organelles.contains('plastome') && organelles.contains('mitogenome')
        }
        .map { sample_id, organelles, genomes, beds ->
            def pt_idx = organelles.indexOf('plastome')
            def mt_idx = organelles.indexOf('mitogenome')
            return tuple(sample_id, genomes[mt_idx], beds[mt_idx], genomes[pt_idx], beds[pt_idx])
        }

    TRANSGENOMES(ch_transgenomes_input)

    SSR_ANALYSIS(PREPARE_GENOME.out.genome_data)
    SSR_SUMMARY(SSR_ANALYSIS.out.sample_csv.collect())

    REPEATS_ANALYSIS(PREPARE_GENOME.out.genome_data)
    REPEATS_SUMMARY(REPEATS_ANALYSIS.out.sample_csv.collect())

    KAKS_ANALYSIS(ch_plastome_kaks.mix(ch_mitogenome_kaks))
    RNA_EDITING(PREPARE_GENOME.out.genome_data)
    
    PHYLOGENY(ch_plastome_kaks.mix(ch_mitogenome_kaks))
}
