process CLEANHEADERS_FULL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.3':
        'biocontainers/coreutils:9.3' }"

    input:
    tuple val(meta) , path(fasta)

    output:
    tuple val(meta), path('*.orthologs.full.fasta'), emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    sed_version = '4.0' // BusyBox sed --version outputs "This is not GNU sed version 4.0"
    if ("${fasta}" == "${prefix}.fasta") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    sed 's|:.*\$||g' ${fasta} > ${prefix}.orthologs.full.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: ${sed_version}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.orthologs.probe.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: ${sed_version}
    END_VERSIONS
    """
}
