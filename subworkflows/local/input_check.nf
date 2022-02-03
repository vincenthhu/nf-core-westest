//
// Check input samplesheet and get read channels
// Input format: 10 cols
// sample,lane,subject_id,gender,phenotype,paternal_id,maternal_id,single_end,fastq_1,fastq_2

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.tsv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .set(sheet)

    ch_case_info = sheet.first()
                        .map { create_case_channels(it) }
    reads        = sheet.map { create_fastq_channels(it) }
    samples      = sheet.map { create_sample_channels(it) }

    emit:
    ch_case_info                              // channel: [ subject ]
    reads                                     // channel: [ val(meta), [ reads ] ]
    samples                                   // channel: [ sample ]

    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq1, fastq2 ] ]
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id             = row.sample + ',' + row.lane
    meta.sample         = row.sample
    meta.lane           = row.lane
    meta.matched_sample = row.matched_sample
    meta.subject        = row.subject_id
    meta.gender         = row.gender
    meta.phenotype      = returnStatus(row.phenotype.toInteger())
    meta.paternal       = row.paternal_id
    meta.maternal       = row.maternal_id
    meta.single_end     = row.single_end.toBoolean()
    meta.read_group     = "\'@RG\\tID:" + row.sample + "\\tSM:${row.sample}\tPL:illimina\tLB:${row.sample}"

    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return array
}

// Function to get a list of metadata (e.g. pedigree, case id) from the sample; [ meta ]
def create_sample_channels(LinkedHashMap row) {
    def sample              = [:]
    sample.subject          = row.subject_id
    sample.gender           = row.gender
    sample.phenotype        = returnStatus(row.phenotype.toInteger())
    sample.paternal         = row.paternal_id
    sample.maternal         = row.maternal_id
    sample.sample           = row.sample
    sample.lane             = row.lane
    sample.matched_sample   = row.matched_sample

    return sample
}

// Function to get a list of metadata (e.g. case id) for the case [ meta ]
def create_case_channel(LinkedHashMap row) {
    def case_info   = [:]
    case_info.id    = row.subject_id

    return case_info
}
