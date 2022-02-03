//
// prepare reference genome files
//

workflow PREPARE_GENOME {
    take:
    prepare_tool_indices // list  : tools to prepare indices for

    main:

    ch_versions = Channel.empty()

    ch_ref_fasta = Channel.empty()
    if (params.fasta) {
        if (params.fasta.endsWith('.fasta') || params.fasta.endsWith('.fa')) {
            if (file(params.fasta).exists()) {
                ch_ref_fasta = file( params.fasta )
            } else {
                exit 1, "Genome reference file: ${params.fasta} does not exist. See --help for more information"
            }
        } else {
            exit 1, "Genome reference: ${params.fasta} has the wrong extension. See --help for more information"
        }
    }

    ch_ref_fai = Channel.empty()
    if (params.fai) {
        if (params.fai.endsWith('.fasta.fai') || params.fasta.endsWith('.fa.fai')) {
            if (file(params.fai).exists()) {
                ch_ref_fai = file( params.fasta )
            } else {
                exit 1, "Genome reference file: ${params.fai} does not exist. See --help for more information"
            }
        } else {
            exit 1, "Genome reference: ${params.fai} has the wrong extension. See --help for more information"
        }
    }

    ch_ref_dict = Channel.empty()
    if (params.dict) {
        if (params.dict.endsWith('.dict')) {
            if (file(params.dict).exists()) {
                ch_ref_dict = file( params.dict )
            } else {
                exit 1, "Genome dictionary: ${params.dict} does not exist. See --help for more information"
            }
        } else {
            exit 1, "Genome dictionary: ${params.dict} has the wrong extension. See --help for more information"
        }
    }

    ch_ref_chromsize = Channel.empty()
    if (params.chromsize) {
        if (params.chromsize.endsWith('.chromsize')) {
            if (file(params.chromsize).exists()) {
                ch_ref_chromsize = file( params.chromsize )
            } else {
                exit 1, "Genome chromsize: ${params.chromsize} does not exist. See --help for more information"
            }
        } else {
            exit 1, "Genome chromsize: ${params.chromsize} has the wrong extension. See --help for more information"
        }
    }

    ch_ref_bwamem2 = Channel.empty()
    if (params.bwamem2) {
        if (params.bwamem2.endsWith('.fasta') || params.bwamem2.endsWith('.fa')) {
            if (!file(params.bwamem2 + '.amb').exists()) {
                exit 1, "BWAMEM2 index: ${params.bwamem2}.amb does not exist. See --help for more information"
            }
            if (!file(params.bwamem2 + '.ann').exists()) {
                exit 1, "BWAMEM2 index: ${params.bwamem2}.ann does not exist. See --help for more information"
            }
            if (!file(params.bwamem2 + '.bwt').exists()) {
                exit 1, "BWAMEM2 index: ${params.bwamem2}.bwt does not exist. See --help for more information"
            }
            if (!file(params.bwamem2 + '.pac').exists()) {
                exit 1, "BWAMEM2 index: ${params.bwamem2}.pac does not exist. See --help for more information"
            }
            if (!file(params.bwamem2 + '.sa').exists()) {
                exit 1, "BWAMEM2 index: ${params.bwamem2}.sa does not exist. See --help for more information"
            }

            if (file(params.bwamem2).exists()) {
                ch_ref_bwamem2 = file( params.bwamem2 )
            } else {
                exit 1, "BWAMEM2 index: ${params.bwa} does not exist. See --help for more information"
            }
        } else {
            exit 1, "BWAMEM2 index: ${params.bwamem2} has the wrong extension. See --help for more information"
        }
    }

    ch_ref_dbsnp = Channel.empty()
    if (params.dbsnp) {
        if (params.dbsnp.endsWith('.vcf.gz')) {
            if (file(params.dbsnp).exists()) {
                ch_ref_dbsnp = file( params.dbsnp )
            } else {
                exit 1, "dbSNP reference: ${params.dbsnp} does not exist. See --help for more information"
            }
        } else {
            exit 1, "dbSNP reference: ${params.dbsnp} has the wrong extension. See --help for more information"
        }
    }

    ch_ref_dbsnp_index = Channel.empty()
    if (params.dbsnp_index) {
        if (params.dbsnp_index.endsWith('.vcf.gz.tbi')) {
            if (file(params.dbsnp_index).exists()) {
                ch_ref_dbsnp_index = file( params.dbsnp_index )
            } else {
                exit 1, "dbSNP reference index: ${params.dbsnp_index} does not exist. See --help for more information"
            }
        } else {
            exit 1, "dbSNP reference index: ${params.dbsnp_index} has the wrong extension. See --help for more information"
        }
    }


    ch_ref_g1000snp = Channel.empty()
    if (params.g1000snp) {
        if (params.g1000snp.endsWith('.vcf.gz')) {
            if (file(params.g1000snp).exists()) {
                ch_ref_g1000snp = file( params.g1000snp )
            } else {
                exit 1, "1000 Genomes phase I SNPs high confidence reference: ${params.g1000snp} does not exist. See --help for more information"
            }
        } else {
            exit 1, "1000 Genomes phase I SNPs high confidence reference: ${params.g1000snp} has the wrong extension. See --help for more information"
        }
    }

    ch_ref_hapmap = Channel.empty()
    if (params.hapmap) {
        if (params.hapmap.endsWith('.vcf.gz')) {
            if (file(params.hapmap).exists()) {
                ch_ref_hapmap = file( params.hapmap )
            } else {
                exit 1, "Hapmap reference: ${params.hapmap} does not exist. See --help for more information"
            }
        } else {
            exit 1, "Hapmap reference: ${params.hapmap} has the wrong extension. See --help for more information"
        }
    }

    ch_ref_mills = Channel.empty()
    if (params.mills) {
        if (params.mills.endsWith('.vcf.gz')) {
            if (file(params.mills).exists()) {
                ch_ref_mills = file( params.mills )
            } else {
                exit 1, "Mills and 1000 Genomes gold standard indels: ${params.mills} does not exist. See --help for more information"
            }
        } else {
            exit 1, "Mills and 1000 Genomes gold standard indels: ${params.mills} has the wrong extension. See --help for more information"
        }
    }

    ch_ref_omni = Channel.empty()
    if (params.omni) {
        if (params.omni.endsWith('.vcf.gz')) {
            if (file(params.omni).exists()) {
                ch_ref_omni = file( params.omni )
            } else {
                exit 1, "1000 Genomes omni 2.5 reference: ${params.omni} does not exist. See --help for more information"
            }
        } else {
            exit 1, "1000 Genomes omni 2.5 reference: ${params.omni} has the wrong extension. See --help for more information"
        }
    }

    ch_ref_axiom = Channel.empty()
    if (params.axiom) {
        if (params.axiom.endsWith('.vcf.gz')) {
            if (file(params.axiom).exists()) {
                ch_ref_axiom = file( params.axiom )
            } else {
                exit 1, "Axiom Exome Plus genotypes all populations poly reference: ${params.axiom} does not exist. See --help for more information"
            }
        } else {
            exit 1, "Axiom Exome Plus genotypes all populations poly reference: ${params.axiom} has the wrong extension. See --help for more information"
        }
    }

    ch_ref_knownindels = Channel.empty()
    if (params.knownindels) {
        if (params.knownindels.endsWith('.vcf.gz')) {
            if (file(params.knownindels).exists()) {
                ch_ref_knownindels = file( params.knownindels )
            } else {
                exit 1, "Known indel reference: ${params.knownindels} does not exist. See --help for more information"
            }
        } else {
            exit 1, "Known indel reference: ${params.knownindels} has the wrong extension. See --help for more information"
        }
    }

    ch_ref_knownindels_index = Channel.empty()
    if (params.knownindels_index) {
        if (params.knownindels_index.endsWith('.vcf.gz.tbi')) {
            if (file(params.knownindels_index).exists()) {
                ch_ref_knownindels_index = file( params.knownindels_index )
            } else {
                exit 1, "Known indels reference index: ${params.knownindels_index} does not exist. See --help for more information"
            }
        } else {
            exit 1, "Known indels reference index: ${params.knownindels_index} has the wrong extension. See --help for more information"
        }
    }


    ch_ref_cosmic = Channel.empty()
    if (params.cosmic) {
        if (params.cosmic.endsWith('.vcf.gz')) {
            if (file(params.cosmic).exists()) {
                ch_ref_cosmic = file( params.cosmic )
            } else {
                exit 1, "COSMIC coding and noncoding chr M sorted reference: ${params.cosmic} does not exist. See --help for more information"
            }
        } else {
            exit 1, "COSMIC coding and noncoding chr M sorted reference: ${params.cosmic} has the wrong extension. See --help for more information"
        }
    }

    ch_ref_germline_resource = Channel.empty()
    if (params.germline_resource) {
        if (params.germline_resource.endsWith('.vcf.gz')) {
            if (file(params.germline_resource).exists()) {
                ch_ref_germline_resource = file( params.germline_resource )
            } else {
                exit 1, "germline_resource: ${params.germline_resource} does not exist. See --help for more information"
            }
        } else {
            exit 1, "germline_resource: ${params.germline_resource} has the wrong extension. See --help for more information"
        }
    }

    ch_ref_germline_resource_index = Channel.empty()
    if (params.germline_resource_index) {
        if (params.germline_resource_index.endsWith('.vcf.gz.tbi')) {
            if (file(params.germline_resource_index).exists()) {
                ch_ref_germline_resource_index = file( params.germline_resource_index )
            } else {
                exit 1, "germline_resource index: ${params.germline_resource_index} does not exist. See --help for more information"
            }
        } else {
            exit 1, "germline_resource index: ${params.germline_resource_index} has the wrong extension. See --help for more information"
        }
    }



    emit:
    ref_fasta               = ch_ref_fasta              // path: genome.fasta
    ref_fai                 = ch_ref_fai                // path: genome.fai
    ref_dict                = ch_ref_dict               // path: genome.dict
    ref_chromsize           = ch_ref_chromsize          // path: genome.chromsize
    ref_bwamem2             = ch_ref_bwamem2            // path: genome.bwamem2
    ref_dbsnp               = ch_ref_dbsnp              // path: genome.dbsnp
    ref_dbsnp_index         = ch_ref_dbsnp_index        // path: genome.dbsnp_index
    ref_g1000snp            = ch_ref_g1000snp           // path: genome.g1000snp
    ref_hapmap              = ch_ref_hapmap             // path: genome.hapmap
    ref_mills               = ch_ref_mills              // path: genome.mills
    ref_omni                = ch_ref_omni               // path: genome.omni
    ref_axiom               = ch_ref_axiom              // path: genome.axiom
    ref_knownindels         = ch_ref_knownindels        // path: genome.knownindels
    ref_knownindels_index   = ch_ref_knownindels_index  // path: genome.knownindels_index
    ref_cosmic              = ch_ref_cosmic             // path: genome.cosmic
    ref_gnomad              = ch_ref_gnomad             // path: genome.gnomad
    ref_gnomad_index        = ch_ref_gnomad_index       // path: genome.gnomad_index

    versions                = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
