process GATK4_PLOTDENOISEDCOPYRATIOS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.4.1--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.4.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(standardizedCR_tsv), path(denoisedCR_tsv)
    path dict

    output:
    tuple val(meta), path("*.denoised.png"          ), emit: denoised_png
    tuple val(meta), path("*.standardizedMAD.txt"   ), emit: standardizedMAD_txt
    tuple val(meta), path("*.denoisedMAD.txt"       ), emit: denoisedMAD_txt
    tuple val(meta), path("*.deltaMAD.txt"          ), emit: deltaMAD_txt
    tuple val(meta), path("*.scaledDeltaMAD.txt"    ), emit: scaledDeltaMAD_txt

    path "versions.yml"                              , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK ModelSegments] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" PlotDenoisedCopyRatios \\
            --standardized-copy-ratios $standardizedCR_tsv \\
            --denoised-copy-ratios  $denoisedCR_tsv \\
            --sequence-dictionary ${dict} \\
            --output './' \\
            --output-prefix ${prefix} \\
            $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
