process GATK4_ANNOTATEINTERVALS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.4.1--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.4.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(interval_list)
    path ref_fasta
    path ref_fai
    path ref_dict

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK AnnotateIntervals] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" AnnotateIntervals \\
        -L $interval_list \\
        -R $ref_fasta
        -O ${prefix}.tsv \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
