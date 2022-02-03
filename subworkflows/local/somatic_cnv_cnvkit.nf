include { CNVKIT_BATCH } from '../../modules/nf-core/modules/cnvkit/batch/main'

workflow SOMATIC_CNV_CNVKIT {
    take:
        tumor_bam,
        tumor_bai,
        normal_bam,
        normal_bai,
        fasta,
        fai,
        dict,
        targets

    main:
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
            meta.tumor_sample = tumor[1]                                   // name list [tumor, tumor, tumor, ...]
            meta.normal_sample = normal[1]                                 // name list [normal, normal, normal, ...]
            [meta, tumor[2], normal[2]]  // tuple val(meta), path list [tumor_bam, tumor_bam, tumor_bam, ...], path list [normal_bam, normal_bam, normal_bam, ...]
        }.set {ch_cnvkit_input}

        CNVKIT_BATCH (
            ch_cnvkit_input,
            fasta,
            targets,
            Channel.of('/')
        )



    emit:
}
