//
// somatic CNV workflow, prepare PON
// https://gatk.broadinstitute.org/hc/en-us/articles/360035531132
// @vincenthhu DEC 2021

include { GATK4_COLLECTREADCOUNTS as GATK4_COLLECTREADCOUNTS_TUMOR       }               from '../../modules/local/gatk4/collectreadcounts'
include { GATK4_COLLECTALLELICCOUNTS as GATK4_COLLECTALLELICCOUNTS_TUMOR }               from '../../modules/local/gatk4/collectalleliccounts'
include { GATK4_DENOISEREADCOUNTS as GATK4_DENOISEREADCOUNTS_TUMOR       }               from '../../modules/local/gatk4/denoisereadcounts'



include { GATK4_CREATEREADCOUNTPANELOFNORMALS }         from '../../modules/nf-core/modules/gatk4/gatk4_createreadcountpanelofnormals/main'




workflow SOMATIC_CNV_PAIR_GATK4 {
    take:
        tumor_bam
        tumor_bai
        processed_interval_list         // path(processed_interval_list, interval_list)
        annotated_interval_list_tsv
        common_sites
        cnv_pon
        ref_fasta
        ref_fai
        ref_dict

    main:
        GATK4_COLLECTREADCOUNTS_TUMOR (
            tumor_bam.join(tumor_bai).combine(processed_interval_list),
            ref_fasta,                                                              // path fasta
            ref_fai,                                                                // path fai
            ref_dict                                                                // path dict
        )
        versions = GATK4_COLLECTREADCOUNTS_TUMOR.out.versions.first()) // channel: [ versions.yml ]

        GATK4_COLLECTALLELICCOUNTS_TUMOR (
            tumor_bam.join(tumor_bai).combine(common_sites),
            ref_fasta,
            ref_fai,
            ref_dict
        )

        GATK4_DENOISEREADCOUNTS_TUMOR (
            GATK4_COLLECTREADCOUNTS_TUMOR.out.hdf5,
            cnv_pon
        )

        GATK4_COLLECTALLELICCOUNTS_TUMOR.out.allelicCounts_tsv.map {
            meta, allelicCounts_tsv ->
            normal_allelicCounts_tsv = allelicCounts_tsv.getParent() + '/' + meta.matched_sample + '.allelicCounts.tsv'
            [meta, allelicCounts_tsv, normal_allelicCounts_tsv]
        }.join(GATK4_COLLECTREADCOUNTS_TUMOR.out.denoisedCR_tsv)
        . set {ch_cnv_somatic_modelsegments_input}

        GATK4_MODELSEGMENTS (
            ch_cnv_somatic_modelsegments_input
        )

        GATK4_CALLCOPYRATIOSEGMENTS (
            GATK4_MODELSEGMENTS.out.cr_seg
        )

        GATK4_PLOTDENOISEDCOPYRATIOS (
            GATK4_DENOISEREADCOUNTS_TUMOR.out.standardizedCR_tsv.join(GATK4_DENOISEREADCOUNTS_TUMOR.out.denoisedCR_tsv)
        )

        GATK4_PLOTMODELEDSEGMENTS (
            GATK4_DENOISEREADCOUNTS_TUMOR.out.denoisedCR_tsv.join(GATK4_MODELSEGMENTS.out.hets_csv).join(GATK4_MODELSEGMENTS.out.modelFinal_seg),
            ref_dict
        )
    emit:




}
