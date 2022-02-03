//
// somatic CNV workflow, prepare PON
// https://gatk.broadinstitute.org/hc/en-us/articles/360035531132
// @vincenthhu DEC 2021

include { GATK4_COLLECTREADCOUNTS as GATK4_COLLECTREADCOUNTS_NORMAL       }                     from '../../modules/local/gatk4/collectreadcounts'
include { GATK4_CREATEREADCOUNTPANELOFNORMALS                             }                     from '../../modules/nf-core/modules/gatk4/gatk4_createreadcountpanelofnormals/main'
include { GATK4_COLLECTALLELICCOUNTS as GATK4_COLLECTALLELICCOUNTS_NORMAL }                     from '../../modules/local/gatk4/collectalleliccounts'



workflow SOMATIC_CNV_PON {
    take:
        normal_bam                      // channel: [val(meta), path(bam)]
        normal_bai                      // channel: [val(meta), path(bai)]
        processed_interval_list         // path(processed_interval_list, interval_list)
        annotated_interval_list_tsv     // path(annotated_interval_list_tsv, in tsv format)
        ref_fasta                       // path(fasta)
        ref_fai                         // path(fai)
        ref_dict                        // path(dict)


    main:
        ch_versions = Channel.empty()

        // Step 1: Build PON from normal bams
        GATK4_COLLECTREADCOUNTS_NORMAL (
            normal_bam.join(normal_bai).combine(processed_interval_list),           // tuple val(meta), path(input), path(input_index), path(intervals)
            ref_fasta,                                                              // path fasta
            ref_fai,                                                                // path fai
            ref_dict                                                                // path dict
        )
        versions = GATK4_COLLECTREADCOUNTS_NORMAL.out.versions.first())


        // consoliate all normal read counts as input for GATK4_CREATEREADCOUNTPANELOFNORMALS
        GATK4_COLLECTREADCOUNTS_NORMAL.out.hdf5.map {
            meta, hdf5 ->
            new_meta = [:]
            new_meta.id = 'CNV_SOMATIC_PON'
            [new_meta, hdf5]
        }.groupTuple(by:0)
        .set {ch_cnv_somatic_PON_input}

        GATK4_CREATEREADCOUNTPANELOFNORMALS (
            ch_cnv_somatic_PON_input,
            annotated_interval_list_tsv
        )
        versions = GATK4_CREATEREADCOUNTPANELOFNORMALS.out.versions.first())

        GATK4_COLLECTALLELICCOUNTS_NORMAL (
            normal_bam.join(normal_bai).combine(common_sites),
            ref_fasta,
            ref_fai,
            ref_dict
        )
        versions = GATK4_COLLECTALLELICCOUNTS_NORMAL.out.versions.first())

    emit:
        pon_hdf5              = GATK4_CREATEREADCOUNTPANELOFNORMALS.out.hdf5
        normal_allelic_tsv    = GATK4_COLLECTALLELICCOUNTS_NORMAL.out.allelic_tsv
        versions              = ch_versions.ifEmpty(null)                            // channel: [ versions.yml ]

}
