//
// somatic SNV variant calling for matched T/N sample pairs using GATK4 Mutect2
// https://gatk.broadinstitute.org/hc/en-us/articles/360035531132
// @vincenthhu DEC 2021

include { GATK4_MUTECT2 }                           from '../../modules/local/gatk4/mutect2'
include { GATK4_LEARNREADORIENTATIONMODEL }         from '../../modules/nf-core/modules/gatk4/learnreadorientationmodel/main'
include { GATK4_GETPILEUPSUMMARIES }                from '../../modules/nf-core/modules/gatk4/getpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION}             from '../../modules/nf-core/modules/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS}                  from '../../modules/nf-core/modules/gatk4/filtermutectcalls/main'

workflow SOMATIC_VC_GATK4 {
    take:
        tumor_bam               // channel: [val(meta), path(bam)]
        tumor_bai               // channel: [val(meta), path(bai)]
        normal_bam
        normal_bai
        target_interval_list     // path(target_interval_list)
        ref_fasta               // path(fasta)
        ref_fai                 // path(fai)
        ref_dict                // path(dict)
        ref_germline_resource
        ref_germline_resource_index
        pon_vcf
        pon_vcf_index

    main:
        ch_versions = Channel.empty()

        // ===================================================
        // get channel for matched bams for each subject
        // ===================================================
        tumor_bam.join( tumor_bai ).map {
            meta, bam, bai ->
            subject = meta.subject
            sample = meta.sample
            [subject, sample, bam, bai]
        }.groupTuple(by:0)
        .set {ch_tumor}

        normal_bam.join( normal_bai ).map {
            meta, bam, bai ->
            subject = meta.subject
            sample = meta.sample
            [subject, sample, bam, bai]
        }.groupTuple(by:0)
        .set {ch_normal}

        ch_tumor.cross(ch_normal).map {
            tumor, normal ->
            meta = [:]
            meta.id = tumor[0]                                             // meta.id = subject
            [meta, tumor[2] + normal[2], tumor[3] + normal[3], normal[1]]  // tuple val(meta), path list [bam, bam, bam], path list [bai, bai, bai], name list [normal, normal, normal]
        }.set {ch_mutect2_input}

        // ====================================================
        // Step 1. Mutect2
        // ====================================================
        // // There are three steps to the filter. First, run Mutect2 with the --f1r2-tar-gz argument.
        // // This creates an output with raw data used to learn the orientation bias model.
        // // Previously this was done by CollectF1R2Counts.
        // // By absorbing it into Mutect2, we eliminated the cost of CollectF1R2Counts with almost no effect on the runtime of Mutect2.
        // // When multiple tumor samples are specified, you only need a single --f1r2-tar-gz output, which contains data for each tumor sample.
        // // gatk Mutect2 -R ref.fasta \
        //     -L intervals.interval_list \
        //     -I tumor.bam \
        //     -germline-resource af-only-gnomad.vcf \
        //     -pon panel_of_normals.vcf   \
        //     --f1r2-tar-gz f1r2.tar.gz \
        //     -O unfiltered.vcf

        // // Matched tumor/normal mode for each subject_id of multiple samples
        // gatk Mutect2 \
        //     -R reference.fa \
        //     -I tumor1.bam \
        //     -I tumor2.bam \
        //     -I normal1.bam \
        //     -I normal2.bam \
        //     -normal normal1_sample_name \
        //     -normal normal2_sample_name \
        //     --germline-resource af-only-gnomad.vcf.gz \
        //     --panel-of-normals pon.vcf.gz \
        //     -O somatic.vcf.gz

        GATK4_MUTECT2 (
            ch_mutect2_input,                                              // tuple val(meta) , path(input) , path(input_index) , val(which_norm)
            false,                                                         // val  run_single
            false,                                                         // val  run_pon
            false,                                                         // val  run_mito
            '',                                                            // val  interval_label
            target_interval_list,                                           // path target_interval_list
            ref_fasta,                                                     // path fasta
            ref_fai,                                                       // path fai
            ref_dict,                                                      // path dict
            ref_germline_resource,                                         // path germline_resource
            ref_germline_resource_index,                                   // path germline_resource_tbi
            pon_vcf,                                                       // path panel_of_normals
            pon_vcf_index                                                  // path panel_of_normals_tbi
        )
        ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions.first())

        // =========================================================
        // Step 2. Filtering
        // =========================================================
        // // 2.1 Next, pass this raw data to LearnReadOrientationModel:
        // gatk LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz

        GATK4_LEARNREADORIENTATIONMODEL ( GATK4_MUTECT2.out.f1r2 )
        ch_versions = ch_versions.mix(GATK4_LEARNREADORIENTATIONMODEL.out.versions.first())

        // // 2.2 Run GetPileupSummaries to summarize read support for a set number of known variant sites.
        // gatk GetPileupSummaries \
        //      -I tumor.bam \
        //      -V chr17_small_exac_common_3_grch38.vcf.gz \
        //      -L chr17_small_exac_common_3_grch38.vcf.gz \
        //      -O getpileupsummaries.table
        GATK4_GETPILEUPSUMMARIES (
            tumor.bam.join(tumor.bai),                                     // tuple val(meta), path(bam), path(bai)
            ref_germline_resource,                                         // path variants
            ref_germline_resource_index,                                   // path variants_tbi
            target_interval_list                                           // path sites
        )
        ch_versions = ch_versions.mix(GATK4_GETPILEUPSUMMARIES.out.versions.first())

        // // 2.3 Estimate contamination with CalculateContamination.
        // gatk CalculateContamination \
        //     -I getpileupsummaries.table \
        //     -tumor-segmentation segments.table \
        //     -O calculatecontamination.table
        GATK4_CALCULATECONTAMINATION(
            GATK4_GETPILEUPSUMMARIES.out.table.combine(Channel.of('/')),   // tuple val(meta), path(pileup), path(matched)
            true                                                           // val segmentout
        )
        ch_versions = ch_versions.mix(GATK4_CALCULATECONTAMINATION.out.versions.first())

        // // 2.4 Finally, pass the learned read orientation model to FilterMutectCallswith the -ob-priors argument:
        GATK4_CALCULATECONTAMINATION.out.contamination.map {
            meta, contamination ->
            new_meta = [:]
            new_meta.id = meta.subject
            [new_meta, contamination]
        }.groupTuple(by:0).set { ch_contamination_list_per_subject}

        GATK4_CALCULATECONTAMINATION.out.segmentation.map {
            meta, segmentation ->
            new_meta = [:]
            new_meta.id = meta.subject
            [new_meta, segmentation]
        }.groupTuple(by:0).set { ch_segmentation_list_per_subject}
        // gatk FilterMutectCalls -V unfiltered.vcf \
        //     [--tumor-segmentation segments.table] \
        //     [--contamination-table contamination.table] \
        //     --ob-priors read-orientation-model.tar.gz \
        //     -O filtered.vcf
        GATK4_FILTERMUTECTCALLS (
            GATK4_MUTECT2.out.vcf.join(GATK4_MUTECT2.out.tbi),                                  // tuple val(meta), path(vcf), path(tbi),
                                 .join(GATK4_MUTECT2.out.stats),                                //                  path(stats),
                                 .join(GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior),      //                  path(orientationbias),
                                 .join(ch_segmentation_list_per_subject),                       //                  path(segmentation),
                                 .join(ch_contamination_list_per_subject),                      //                  path(contaminationfile),
                                 .combine(Channel.of(false))                                    //                  val(contaminationest)
            ref_fasta,                                                                          // path fasta
            ref_fai,                                                                            // path fai
            ref_dict                                                                            // path dict
        )
        ch_versions = ch_versions.mix(GATK4_FILTERMUTECTCALLS.out.versions.first())


    emit:

        vcf                   = GATK4_FILTERMUTECTCALLS.out.vcf
        tbi                   = GATK4_FILTERMUTECTCALLS.out.tbi
        stats                 = GATK4_FILTERMUTECTCALLS.out.stats

        versions              = ch_versions.ifEmpty(null)                            // channel: [ versions.yml ]

}
