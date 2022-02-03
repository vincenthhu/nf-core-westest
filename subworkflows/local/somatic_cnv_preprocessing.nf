//
// somatic SNV variant calling for matched T/N sample pairs using GATK4 Mutect2
// https://gatk.broadinstitute.org/hc/en-us/articles/360035531132
// @vincenthhu DEC 2021

include { GATK4_PREPROCESSINTERVALS }               from '../../modules/local/gatk4/preprocessintervals'
include { GATK4_ANNOTATEINTERVALS   }               from '../../modules/local/gatk4/annotateintervals'
include { GATK4_INTERVALLISTTOBED   }               from '../../modules/local/gatk4/intervallisttobed'
include { GATK4_BEDTOINTERVALLIST   }               from '../../modules/nf-core/modules/gatk4/bedtointervallist/main'

include { BEDTOOLS_INTERSECT        }               from '../../modules/nf-core/modules/bedtools/intersect/main'
include { BEDTOOLS_SORT             }               from '../../modules/nf-core/modules/bedtools/sort/main'



workflow SOMATIC_CNV_PREPROCESSING {
    take:
        target_interval_list                 // path target interva_list of WES
        blacklist_intervals                 // path blacklist_intervals
        ref_fasta
        ref_fai
        ref_dict
        ref_dnsnp
        ref_dbsnp_index

    main:
        GATK4_PREPROCESSINTERVALS (
            [[id:'CNV'], target_interval_list, blacklist_intervals]
            250,            // padding 250
            100000,         // bin-length 100k
            ref_fasta,
            ref_fai
            ref_dict
        )
        versions = GATK4_PREPROCESSINTERVALS.out.versions.first())

        GATK4_ANNOTATEINTERVALS (
            GATK4_PREPROCESSINTERVALS.out.interval_list,
            ref_fasta,
            ref_fai
            ref_dict
        )
        versions = GATK4_ANNOTATEINTERVALS.out.versions.first())

        GATK4_INTERVALLISTTOBED (
            GATK4_PREPROCESSINTERVALS.out.interval_list
        )
        versions = GATK4_INTERVALLISTTOBED.out.versions.first())

        BEDTOOLS_INTERSECT (
            GATK4_INTERVALLISTTOBED.out.bed.combine(ref_dnsnp),
            'bed'
        )
        versions = BEDTOOLS_INTERSECT.out.versions.first())

        BEDTOOLS_SORT (
            BEDTOOLS_INTERSECT.out.intersect,
            'bed' )
        versions = BEDTOOLS_SORT.out.versions.first())

        GATK4_BEDTOINTERVALLIST (
            BEDTOOLS_SORT.out.sorted,
            ref_dict)
        versions = GATK4_BEDTOINTERVALLIST.out.versions.first())


    emit:
        processed_interval_list              = GATK4_PREPROCESSINTERVALS.out.interval_list  // tuple val(meta), path('*.interval_list')
        annotated_interval_list_tsv          = GATK4_ANNOTATEINTERVALS.out.tsv              // tuple val(meta), path('*.tsv')
        common_sites                         = GATK4_BEDTOINTERVALLIST.out.interval_list    // tuple val(meta), path('*.interval_list')

}
