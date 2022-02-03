//
// A quality check subworkflow for processed bams.
//

include { GATK4_HAPLOTYPECALLER  } from '../../modules/nf-core/modules/gatk4/haplotypecaller/main'
include { GATK4_INDEXFEATUREFILE } from '../../modules/nf-core/modules/gatk4/indexfeaturefile/main'
include { GATK4_GENOMICSDBIMPORT } from '../../modules/nf-core/modules/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS } from '../../modules/nf-core/modules/gatk4/genotypegvcfs/main'


workflow GERMLINE_VC_GATK4 {

    take:
        bam                     // channel: [ val(meta), path(bam) ]
        bai                     // channel: [ val(meta), path(bai) ]
        fasta                   // path: genome.fasta
        fai                     // path: genome.fai
        dict                    // path: genome.dict
        dbsnp                   // path: genome.dbsnp
        dbsnp_index             // path: genome.dbsnp_index
        target_interval_list    // path: target_interval_list


    main:
        ch_versions = Channel.empty()

        // run gatk4 haplotypecaller for germline variant calling
        GATK4_HAPLOTYPECALLER (
            bam.join(bai).combine(target_interval_list),                         // tuple val(meta), path(input), path(input_index), path(intervals)
            fasta,                                                              // path fasta
            fai,                                                                // path fai
            dict,                                                               // path dict
            dbsnp,                                                              // path dbsnp
            dbsnp_tbi                                                           // path dbsnp_tbi
        )
        ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())

        // GATK4_INDEXFEATUREFILE (GATK4_HAPLOTYPECALLER.out.vcf)               // tuple val(meta), path(feature_file)
        // ch_versions = ch_versions.mix(GATK4_INDEXFEATUREFILE.out.versions.first())


        GATK4_HAPLOTYPECALLER.out.vcf.map {
            meta, vcf ->
            new_meta = [:]                                                      // clone to avoid overriding the global meta
            new_meta.id = 'GERMLINE'                                            // access the .id attribute of meta
            [new_meta, vcf]
        }                                                                       // end the closure to return newly modified channel
        .groupTuple(by: 0)                                                      // group VCF paths into [ [GERMLINE], [VCF path, VCF path, ..] ]
        .set{ ch_gendb_vcf_input }                                              // create a new multi-channel named VCFs

        GATK4_HAPLOTYPECALLER.out.tbi.map {
            meta, tbi ->
            new_meta = [:]                                                      // clone to avoid overriding the global meta
            new_meta.id = 'GERMLINE'                                            // access the .id attribute of meta
            [new_meta, tbi]
        }                                                                       // end the closure to return newly modified channel
        .groupTuple(by: 0)                                                      // group tbi paths into [ [GERMLINE], [tbi path, tbi path, ..] ]
        .set{ ch_gendb_tbi_input }                                              // create a new multi-channel named tbi


        GATK4_GENOMICSDBIMPORT (
            ch_gendb_vcf_input.join(ch_gendb_tbi_input)     // tuple val(meta), path(vcf), path(tbi), path(intervalfile), val(intervalval), path(wspace)
                                .combine(target_interval_list)
                                .combine(Channel.of(null))
                                .combine(Channel.of('/')),                  //NOTE: use Channel.of('/') for null path in tuple!!!!!!!
            false,                                                              // val run_intlist
            false,                                                              // val run_updatewspace
            false                                                               // val input_map
        )
        ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions.first())


        GATK4_GENOTYPEGVCFS (
            GATK4_GENOMICSDBIMPORT.out.genomicsdb                               // tuple val(meta), path(gvcf), path(gvcf_index), path(intervals)
                                        .combine(Channel.of('/'))
                                        .combine(target_interval_list)
            fasta                                                               // path  fasta
            fai                                                                 // path  fasta_index
            dict                                                                // path  fasta_dict
            dbsnp                                                               // path  dbsnp
            dbsnp_index                                                         // path  dbsnp_index
        )
        ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions.first())


    emit:
        vcf                     = GATK4_GENOTYPEGVCFS.out.vcf                   // channel: [ val(meta), path(*.vcf.gz) ]
        tbi                     = GATK4_GENOTYPEGVCFS.out.tbi                   // channel: [ val(meta), path(*.tbi) ]
        versions                = ch_versions.ifEmpty(null)                     // channel: [ versions.yml ]
}
