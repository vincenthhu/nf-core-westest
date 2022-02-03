//
// PREPARE PON for mutect2 somatic
// https://gatk.broadinstitute.org/hc/en-us/articles/360035531132
// DEC 2021

include { GATK4_MUTECT2 as GATK4_MUTECT2_PON                                         } from '../../modules/nf-core/modules/gatk4/mutect2/main'
include { GATK4_GENOMICSDBIMPORT as GATK4_GENOMICSDBIMPORT_PON                       } from '../../modules/local/gatk4/genomicsdbimport'
include { GATK4_CREATESOMATICPANELOFNORMALS as GATK4_CREATESOMATICPANELOFNORMALS_PON } from '../../modules/local/gatk4/createsomaticpanelofnormals'

workflow PREPARE_PON {
    take:
        bam                             // channel: [val(meta), path(bam)]
        bai                             // channel: [val(meta), path(bai)]
        target_interval_list             // path(target_interval_list)
        ref_fasta                       // path(fasta)
        ref_fai                         // path(fai)
        ref_dict                        // path(dict)
        ref_germline_resource
        ref_germline_resource_index


    main:
        ch_versions = Channel.empty()

        //
        // The three steps to create a panel of normals are:
        //
        // Step 1: Run Mutect2 in tumor-only mode for each normal sample:
        // gatk Mutect2 -R reference.fasta -I normal1.bam --max-mnp-distance 0 -O normal1.vcf.gz
        GATK4_MUTECT2_PON (
            bam.join(bai).combine(Channel.of(null)),            // tuple val(meta) , path(input) , path(input_index) , val(which_norm)
            false,                                              // val  run_single
            true,                                               // val  run_pon
            false,                                              // val  run_mito
            '',                                                 // val  interval_label
            ref_fasta,                                          // path fasta
            ref_fai,                                            // path fai
            ref_dict,                                           // path dict
            Channel.of('/'),                                    // path germline_resource
            Channel.of('/'),                                    // path germline_resource_tbi
            Channel.of('/'),                                    // path panel_of_normals
            Channel.of('/')                                     // path panel_of_normals_tbi

        )
        ch_versions = ch_versions.mix(GATK4_MUTECT2_PON.out.versions.first())

        // Step 2: Create a GenomicsDB from the normal Mutect2 calls:
        // consolidate vcf of normal samples into list
        GATK4_MUTECT2_PON.out.vcf.map {
            meta, vcf ->
            new_meta = [:]
            new_meta.id = 'PON_db'                           // access the .id attribute of meta, will be the name of gendb
            [new_meta, vcf]
        }                                                       // end the closure to return newly modified channel
        .groupTuple(by: 0)                                      // group vcf paths into list [ [id:PON], [vcf path, vcf path, ..] ]
        .set{ ch_pon_gendb_vcf_input }                          // create a new multi-channel named VCFs

        GATK4_MUTECT2_PON.out.tbi.map {
            meta, tbi ->
            new_meta = [:]
            new_meta.id = 'PON_db'                           // access the .id attribute of meta
            new_meta.phenotype = meta.phenotype
            [new_meta, tbi]
        }                                                       // end the closure to return newly modified channel
        .groupTuple(by: 0)                                      // group tbi paths into list [ [id:PON], [tbi path, tbi path, ..] ]
        .set{ ch_pon_gendb_tbi_input }                          // create a new multi-channel named tbi

        // gatk GenomicsDBImport -R reference.fasta -L intervals.interval_list \
        // --genomicsdb-workspace-path pon_db \
        // -V normal1.vcf.gz \
        // -V normal2.vcf.gz \
        // -V normal3.vcf.gz

        GATK4_GENOMICSDBIMPORT_PON (
            ch_pon_gendb_vcf_input
            .join(ch_pon_gendb_tbi_input)
            .combine(target_interval_list)                     // tuple val(meta), path(vcf), path(tbi), path(intervalfile)
            ref_fasta,
            ref_fai,
            ref_dict
        )
        ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT_PON.out.versions.first())

        //
        // Step 3: Combine the normal calls using CreateSomaticPanelOfNormals:

        //gatk CreateSomaticPanelOfNormals -R reference.fasta \
        //  --germline-resource af-only-gnomad.vcf.gz \
        //  -V gendb://pon_db \
        //  -O pon.vcf.gz

        GATK4_CREATESOMATICPANELOFNORMALS_PON (
            GATK4_GENOMICSDBIMPORT_PON.out.genomicsdb,
            ref_fasta,
            ref_fai,
            ref_dict,
            ref_germline_resource,
            ref_germline_resource_index
        )
        ch_versions = ch_versions.mix(GATK4_CREATESOMATICPANELOFNORMALS_PON.out.versions.first())

    emit:

        pon_vcf               = GATK4_CREATESOMATICPANELOFNORMALS_PON.out.vcf        // channel: [ val(meta), [ vcf ] ]
        pon_vcf_index         = GATK4_CREATESOMATICPANELOFNORMALS_PON.out.tbi        // channel: [ val(meta), [ tbi ] ]

        versions              = ch_versions.ifEmpty(null)                            // channel: [ versions.yml ]

}


