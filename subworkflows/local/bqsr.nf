//
// Run the BQSR procedure given Deduplicated bam
//

include { GATK4_BASERECALIBRATOR } from '../../modules/nf-core/modules/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR } from '../../modules/nf-core/modules/gatk4/applybqsr/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/modules/samtools/index/main'

workflow BQSR {
    take:
        bam                     // channel: [val(meta), path(bam)]
        bai                     // channel: [val(meta), path(bai)]
        target_interval_list     // path(target_interval_list)
        ref_fasta               // path(fasta)
        ref_fai                 // path(fai)
        ref_dict                // path(dict)
        ref_dbsnp               // path(ref_dnsnp)
        ref_dbsnp_index         // path(ref_dbsnp_index)
        ref_knownindels         // path(ref_knownindels)
        ref_knownindels_index   // path(ref_knownindels_index)

    main:
        ch_versions = Channel.empty()

        GATK4_BASERECALIBRATOR (
            bam.join(bai).combine(target_interval_list),        // tuple val(meta), path(input), path(input_index), path(intervals)
            ref_fasta,                                          // path fasta
            ref_fai,                                            // path fai
            ref_dict,                                           // path dict
            ref_dbsnp.mix(ref_knownindels),                     // path knownSites
            ref_dbsnp_index.mix(ref_knownindels_index)          // path knownSites_tbi
        )
        ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions.first())

        GATK4_APPLYBQSR (
            bam
            .join(bai)
            .join(GATK4_BASERECALIBRATOR.out.table)
            .combine(target_interval_list),                      // tuple val(meta), path(input), path(input_index), path(bqsr_table), path(intervals)
            ref_fasta,                                           // path  fasta
            ref_fai,                                             // path  fai
            ref_dict                                             // path  dict
        )
        ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions.first())

        SAMTOOLS_INDEX ( GATK4_APPLYBQSR.out.bam )

    emit:

        bqsr_bam               = GATK4_APPLYBQSR.out.bam        // channel: [ val(meta), [ marked_bam ] ]
        bqsr_bai               = SAMTOOLS_INDEX.out.bai         // channel: [ val(meta), [ marked_bai ] ]

        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]

}
