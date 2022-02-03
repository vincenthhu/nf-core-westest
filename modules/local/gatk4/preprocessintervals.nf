process GATK4_PREPROCESSINTERVALS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.4.1--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.4.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(interval_list), path(blacklist_intervals)
    val padding
    val bin_length
    path ref_fasta
    path ref_fai
    path ref_dict

    output:
    tuple val(meta), path("*.interval_list"), emit: interval_list
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (padding) {
        padding_option = "--padding ${padding}"
    }
    if (bin_length) {
        bin_length_option = "--bin-length ${bin_length}"
    }
    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK PreprocessIntervals] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" PreprocessIntervals \\
        -L $interval_list \\
        -XL $blacklist_intervals \\
        -R $ref_fasta
        -O ${prefix}.interval_list \\
        $padding_option \\
        $bin_length_option \\
        --interval-merging-rule OVERLAPPING_ONLY
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
