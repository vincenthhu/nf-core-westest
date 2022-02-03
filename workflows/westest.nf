/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowWestest.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
checkPathParamList = [
    params.input, params.multiqc_config,
    params.fasta, params.fasta_fai, params.dict, params.chromsize,
    params.bwa,
    params.dbsnp, params.g1000snp, params.hapmap, params.mills, params.omni, params.axiom, params.knownindels, params.cosmic, params.gnomad
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Local modules
//



//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                 } from '../subworkflows/local/input_check'
include { PREPARE_GENOME              } from '../subworkflows/local/prepare_genome'
include { PREPARE_BED                 } from '../subworkflows/local/prepare_bed'
include { BQSR                        } from '../subworkflows/local/bqsr'
include { GERMLINE_VC_GATK4           } from '../subworkflows/local/germline_vc_gatk4'
include { PREPARE_PON                 } from '../subworkflows/local/prepare_pon'
include { SOMATIC_VC_GATK4            } from '../subworkflows/local/somatic_vc_gatk4'
include { SOMATIC_CNV_PREPROCESSING   } from '../subworkflows/local/somatic_cnv_preprossing'
include { SOMATIC_CNV_PON             } from '../subworkflows/local/somatic_cnv_pon'
include { SOMATIC_CNV_PAIR_GATK4      } from '../subworkflows/local/somatic_cnv_pair_gatk4'
include { SOMATIC_CNV_CNVKIT          } from '../subworkflows/local/somatic_cnv_cnvkit'


include { ALIGN_BWAMEM2               } from '../subworkflows/nf-core/align_bwamem2'
/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { SAMTOOLS_SORT               } from '../modules/nf-core/modules/samtools/sort/main'
include { PICARD_SORTSAM              } from '../modules/nf-core/modules/picard/sortsam/main'
include { PICARD_MARKDUPLICATES       } from '../modules/nf-core/modules/picard/markduplicates/main'
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { UMITOOLS_EXTRACT            } from '../modules/nf-core/modules/umitools/extract/main'
include { TRIMGALORE                  } from '../modules/nf-core/modules/trimgalore/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow WESTEST {

    ch_versions = Channel.empty()

    /*
    ================================================================================
                                    CHECKING REFERENCES
    ================================================================================
    */
    PREPARE_GENOME ()
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    ch_target_bed = Channel.empty()
    ch_interval_list = Channel.empty()
    if (params.target_bed) {
        PREPARE_BED (
            params.target_bed,
            params.dict
        )

        ch_target_bed = PREPARE_BED.out.idx
        ch_interval_list = PREPARE_BED.out.interval_list
    }
    /*
    ================================================================================
                                    Preprocess input
    ================================================================================
    */

    //
    // Read in csv meta sheet
    //
    INPUT_CHECK (ch_input)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // Run FastQC
    //
    FASTQC( INPUT_CHECK.out.reads )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Run trimming
    //
    params.save_trimmed = true      // save trimmed fastq files
    TRIMGALORE ( INPUT_CHECK.out.reads)
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())

    //
    // Mapping with BWA mem2
    //
    ch_marked_bam = Channel.empty()
    ch_marked_bai = Channel.empty()

    // fastqs from different lanes of same sample are merged in ALIGN_BWAMEM2()!!!!!!!!!!!
    if (params.aligner == 'bwamem2') {
        ALIGN_BWAMEM2 (
            TRIMGALORE.out.reads,
            PREPARE_GENOME.out.ref_bwamem2
        )

        ch_marked_bam = ALIGN_BWAMEM2.out.marked_bam
        ch_marked_bai = ALIGN_BWAMEM2.out.marked_bai

        ch_versions = ch_versions.mix(ALIGN_BWAMEM2.out.versions.first())
    }

    //
    // Recalibration
    //
    BQSR (
        ch_marked_bam,
        ch_marked_bai,
        ch_interval_list.map { it [1] },
        PREPARE_GENOME.out.fasta
        PREPARE_GENOME.out.fasta_fai
    )
    ch_versions = ch_versions.mix(BQSR.out.versions.ifEmpty(null))

    //
    // BAM quality check
    //
    QC_BAM (
        BQSR.out.bqsr_bam,
        PREPARE_GENOME.out.ref_fasta,
        PREPARE_GENOME.out.ref_chromsize
    )
    ch_versions = ch_versions.mix(QC_BAM.out.versions.ifEmpty(null))
    /*
    ================================================================================
                                    Germline SNVariant calling
    ================================================================================
    */
    GERMLINE_VC_GATK4 (
        BQSR.out.bam,
        BQSR.out.bai,
        PREPARE_GENOME.out.ref_fasta,
        PREPARE_GENOME.out.ref_fai,
        PREPARE_GENOME.out.ref_dict,
        PREPARE_GENOME.out.ref_dnsnp,
        PREPARE_GENOME.out.ref_dbsnp_index
        ch_interval_list.map { it [1] },
    )
    ch_versions = ch_versions.mix(GERMLINE_VC_GATK4.out.versions.first())

    /*
    ================================================================================
                                    Somatic SNVariant calling
    ================================================================================
    */

    ch_normal_bam = Channel.empty()
    ch_normal_bai = Channel.empty()
    ch_tumor_bam = Channel.empty()
    ch_tumor_bai = Channel.empty()

    def criteria = branchCriteria {
                    normal: it.meta.phenotype == 0
                    tumor:  it.meta.phenotype == 1
                    }

    BQSR.out.bqsr_bam.branch(criteria).set { ch_SV_bam }
    ch_normal_bam = ch_SV_bam.normal
    ch_tumor_bam = ch_SV_bam.tumor

    BQSR.out.bqsr_bai.branch(criteria).set { ch_SV_bai }
    ch_normal_bai = ch_SV_bai.normal
    ch_tumor_bai = ch_SV_bai.tumor

    //
    // Build PON
    //
    PREPARE_PON (
        ch_normal_bam,                      // channel: [val(meta), path(bam)]
        ch_normal_bai,                      // channel: [val(meta), path(bai)]
        ch_interval_list.map { it [1] },    // path(target_interval_list)
        PREPARE_GENOME.out.ref_fasta,
        PREPARE_GENOME.out.ref_fai,
        PREPARE_GENOME.out.ref_dict
    )
    ch_versions = ch_versions.mix(PREPARE_PON.out.versions.first())

    //
    // variant calling
    //
    SOMATIC_VC_GATK4 (
        ch_tumor_bam,
        ch_tumor_bai,
        ch_normal_bam,
        ch_normal_bai,
        ch_interval_list.map {it[1]},
        PREPARE_GENOME.out.ref_fasta,
        PREPARE_GENOME.out.ref_fai,
        PREPARE_GENOME.out.ref_dict,
        PREPARE_GENOME.out.ref_germline_resource,
        PREPARE_GENOME.out.ref_germline_resource_index,
        PREPARE_PON.out.pon_vcf,
        PREPARE_PON.out.pon_vcf_index
    )
    ch_versions = ch_versions.mix(SOMATIC_VC_GATK4.out.versions.first())


    /*
    ================================================================================
                                    Germline CNV calling
    ================================================================================
    */



    /*
    ================================================================================
                                    Somatic CNV calling
    ================================================================================
    */

    SOMATIC_CNV_PREPROCESSING (
        ch_interval_list.map { it[1] }
        PREPARE_GENOME.out.blacklist_intervals
    )

    SOMATIC_CNV_PON (
        ch_normal_bam,
        ch_normal_bai,
        SOMATIC_CNV_PREPROCESSING.out.processed_interval_list.map {it[1]},
        SOMATIC_CNV_PREPROCESSING.out.annotated_interval_list_tsv.map {it[1]},
        PREPARE_GENOME.out.ref_fasta,
        PREPARE_GENOME.out.ref_fai,
        PREPARE_GENOME.out.ref_dict
    )

    SOMATIC_CNV_PAIR_GATK4 (
        ch_tumor_bam,
        ch_tumor_bai,
        SOMATIC_CNV_PREPROCESSING.out.processed_interval_list.map {it[1]},
        SOMATIC_CNV_PREPROCESSING.out.annotated_interval_list_tsv.map {it[1]},
        SOMATIC_CNV_PREPROCESSING.out.common_sites.map {it[1]},
        SOMATIC_CNV_PON.out.pon_hdf5,
        PREPARE_GENOME.out.ref_fasta,
        PREPARE_GENOME.out.ref_fai,
        PREPARE_GENOME.out.dict
    )


    /*
    =======================================================k=========================
                                    CNV calling with CNVkit
    ================================================================================
    */

    SOMATIC_CNV_CNVKIT (
        ch_tumor_bam,
        ch_tumor_bai,
        ch_normal_bam,
        ch_normal_bai,
        PREPARE_GENOME.out.ref_fasta,
        PREPARE_GENOME.out.ref_fai,
        PREPARE_GENOME.out.ref_dict,
        PREPARE_BED.out.refFlat.map { it[1] }
    )



























    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )


    //
    //
    //










    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowWestest.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    Extra functions
========================================================================================
*/
// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
}

// Check parameter existence
def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
        log.warn "Unknown parameter: ${it}"
        return false
    }
    return true
}

// Compare each parameter with a list of parameters
def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
}

// Channeling the TSV file containing BAM.
// Format is: "subject gender status sample bam bai"
def extractBam(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            checkNumberOfItem(row, 6)
            def idPatient = row[0]
            def gender    = row[1]
            def status    = returnStatus(row[2].toInteger())
            def idSample  = row[3]
            def bamFile   = returnFile(row[4])
            def baiFile   = returnFile(row[5])

            if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
            if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"

            return [idPatient, gender, status, idSample, bamFile, baiFile]
        }
}

// Create a channel of germline FASTQs from a directory pattern: "my_samples/*/"
// All FASTQ files in subdirectories are collected and emitted;
// they must have _R1_ and _R2_ in their names.
def extractFastqFromDir(pattern) {
    def fastq = Channel.create()
    // a temporary channel does all the work
    Channel
        .fromPath(pattern, type: 'dir')
        .ifEmpty { error "No directories found matching pattern '${pattern}'" }
        .subscribe onNext: { sampleDir ->
            // the last name of the sampleDir is assumed to be a unique sample id
            sampleId = sampleDir.getFileName().toString()

            for (path1 in file("${sampleDir}/**_R1_*.fastq.gz")) {
                assert path1.getName().contains('_R1_')
                path2 = file(path1.toString().replace('_R1_', '_R2_'))
                if (!path2.exists()) error "Path '${path2}' not found"
                (flowcell, lane) = flowcellLaneFromFastq(path1)
                String random = org.apache.commons.lang.RandomStringUtils.random(8, true, true) // random string to avoid duplicate names
                patient = sampleId
                gender = 'ZZ'  // unused
                status = 0  // normal (not tumor)
                rgId = "${flowcell}.${sampleId}.${lane}.${random}"
                result = [patient, gender, status, sampleId, rgId, path1, path2]
                fastq.bind(result)
            }
    }, onComplete: { fastq.close() }
    fastq
}

// Extract gender and status from Channel
def extractInfos(channel) {
    def genderMap = [:]
    def statusMap = [:]
    channel = channel.map{ it ->
        def idPatient = it[0]
        def gender = it[1]
        def status = it[2]
        def idSample = it[3]
        genderMap[idPatient] = gender
        statusMap[idPatient, idSample] = status
        [idPatient] + it[3..-1]
    }
    [genderMap, statusMap, channel]
}

// Channeling the TSV file containing FASTQ or BAM
// Format is: "subject gender status sample lane fastq1 fastq2"
// or: "subject gender status sample lane bam"
def extractFastq(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            def idPatient  = row[0]
            def gender     = row[1]
            def status     = returnStatus(row[2].toInteger())
            def idSample   = row[3]
            def idRun      = row[4]
            def file1      = returnFile(row[5])
            def file2      = "null"
            if (hasExtension(file1, "fastq.gz") || hasExtension(file1, "fq.gz") || hasExtension(file1, "fastq") || hasExtension(file1, "fq")) {
                checkNumberOfItem(row, 7)
                file2 = returnFile(row[6])
                if (!hasExtension(file2, "fastq.gz") && !hasExtension(file2, "fq.gz")  && !hasExtension(file2, "fastq") && !hasExtension(file2, "fq"))
                    exit 1, "File: ${file2} has the wrong extension. See --help for more information"
                if (hasExtension(file1, "fastq") || hasExtension(file1, "fq") || hasExtension(file2, "fastq") || hasExtension(file2, "fq")) {
                    exit 1, "We do recommend to use gziped fastq file to help you reduce your data footprint."
                }
            }
            //else if (hasExtension(file1, "bam")) checkNumberOfItem(row, 6)
            else "No recognisable extention for input file: ${file1}"

            [idPatient, gender, status, idSample, idRun, file1, file2]
        }
}

// Channeling the TSV file containing mpileup
// Format is: "subject gender status sample pileup"
def extractPileup(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            checkNumberOfItem(row, 5)
            def idPatient = row[0]
            def gender    = row[1]
            def status    = returnStatus(row[2].toInteger())
            def idSample  = row[3]
            def mpileup   = returnFile(row[4])

            if (!hasExtension(mpileup, "pileup")) exit 1, "File: ${mpileup} has the wrong extension. See --help for more information"

            return [idPatient, gender, status, idSample, mpileup]
        }
}

// Channeling the TSV file containing Recalibration Tables.
// Format is: "subject gender status sample bam bai recalTable"
def extractRecal(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            checkNumberOfItem(row, 7)
            def idPatient  = row[0]
            def gender     = row[1]
            def status     = returnStatus(row[2].toInteger())
            def idSample   = row[3]
            def bamFile    = returnFile(row[4])
            def baiFile    = returnFile(row[5])
            def recalTable = returnFile(row[6])

            if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
            if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"
            if (!hasExtension(recalTable, "recal.table")) exit 1, "File: ${recalTable} has the wrong extension. See --help for more information"

            [idPatient, gender, status, idSample, bamFile, baiFile, recalTable]
    }
}

// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    InputStream fileStream = new FileInputStream(path.toFile())
    InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream)
    Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
    BufferedReader buffered = new BufferedReader(decoder)
    def line = buffered.readLine()
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(' ')[0].split(':')
    String fcid
    int lane
    if (fields.size() == 7) {
        // CASAVA 1.8+ format
        fcid = fields[2]
        lane = fields[3].toInteger()
    } else if (fields.size() == 5) {
        fcid = fields[0]
        lane = fields[1].toInteger()
    }
    [fcid, lane]
}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Return file if it exists
def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}

// Remove .ann .gz and .vcf extension from a VCF file
def reduceVCF(file) {
    return file.fileName.toString().minus(".ann").minus(".vcf").minus(".gz")
}

// Return status [0,1]
// 0 == Normal, 1 == Tumor
def returnStatus(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it
}
/*
========================================================================================
    THE END
========================================================================================
*/
