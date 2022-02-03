process GATK4_GENOMICSDBIMPORT {
    tag "$meta.id"
    label 'process_low'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "bioconda::gatk4=4.2.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.4.1--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.4.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(intervalfile)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("${prefix}")      , emit: genomicsdb
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    inputs_command = "${'-V ' + vcf.join(' -V ')}"
    dir_command = "--genomicsdb-workspace-path ${prefix}"
    intervals_command = " -L ${intervalfile} "

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK GenomicsDBImport] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    gatk --java-options "-Xmx${avail_mem}g" GenomicsDBImport \\
        -R ${fasta} \\
        $dir_command \\
        $intervals_command \\
        $inputs_command \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
