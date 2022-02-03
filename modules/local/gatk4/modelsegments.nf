
process GATK4_MODELSEGMENTS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.4.1--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.4.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_allelicCounts_tsv), path(normal_allelicCounts_tsv), path(tumor_denoisedCR_tsv)

    output:
    tuple val(meta), path("*.hets.csv"), emit: hets_csv
    tuple val(meta), path("*.hets.normal.csv"), emit: hets_normal_csv
    tuple val(meta), path("*.cr.seg"), emit: cr_seg
    tuple val(meta), path("*.cr.igv.seg"), emit: cr_igv_seg
    tuple val(meta), path("*.af.igv.seg"), emit: af_igv_seg
    tuple val(meta), path("*.modelBegin.seg"), emit: modelBegin_seg
    tuple val(meta), path("*.modelBegin.cr.param"), emit: modelBegin_cr_param
    tuple val(meta), path("*.modelBegin.af.param"), emit: modelBegin_af_param
    tuple val(meta), path("*.modelFinal.seg"), emit: modelFinal_seg
    tuple val(meta), path("*.modelFinal.cr.param"), emit: modelFinal_cr_param
    tuple val(meta), path("*.modelFinal.af.param"), emit: modelFinal_af_param

    path "versions.yml"           , emit: versions

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
    gatk --java-options "-Xmx${avail_mem}g" ModelSegments \\
            --denoised-copy-ratios $tumor_denoisedCR_tsv \\
            --allelic-counts $tumor_allelicCounts_tsv \\
            --normal-allelic-counts $normal_allelicCounts_tsv \\
            --maximum-number-of-segments-per-chromosome 0 \\
            --output './' \\
            --output-prefix ${prefix} \\
            $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
