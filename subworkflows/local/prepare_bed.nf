//
// Prepare reference bed files
//

include { TABIX_TABIX as TABIX_PT } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_BGZIPTABIX as TABIX_PBT } from '../../modules/nf-core/modules/tabix/bgziptabix/main'
include { GATK4_BEDTOINTERVALLIST } from '../../modules/nf-core/modules/gatk4/bedtointervallist/main'

workflow PREPARE_BED {
    take:
        bed_name        // file: bed file
        dict            // path: ref_dict
        refFlat_name

    main:
        tab_out = Channel.empty()
        if (bed_name) {
            bed_file = file(bed_name)
            id       = bed_name.split('/')[-1]
            ch_bed   = Channel.fromList([[['id':id], bed_file]])

            if ( bed_name.endsWith(".gz") && file(bed_name, checkIfExists:true) ) {
                tbi_out = TABIX_PT (ch_bed).tbi
                tab_out = ch_bed.join(tbi_out)
                ch_interval_list = GATK4_BEDTOINTERVALLIST(ch_bed, dict).out.interval_list
            } else if ( file(bed_name, checkIfExists:true) ) {
                tab_out = TABIX_PBT (ch_bed).gz_tbi
                ch_interval_list = GATK4_BEDTOINTERVALLIST(ch_bed, dict).out.interval_list

                )
            } else {
                exit 1, "ERROR: Please check target bed file does not exist!\n${bed_file}"
            }
        }

        ch_refFlat = Channel.empty()
        if (refFlat_name) {
            refFlat_file = file(refFlat_name)
            id = refFlat_name.split('/')[-1]
            ch_refFlat = Channel.fromList([[['id':id], refFlat_file]])
        }

    emit:
        idx             = tab_out
        interval_list   = ch_interval_list
        refFlat         = ch_refFlat
}
