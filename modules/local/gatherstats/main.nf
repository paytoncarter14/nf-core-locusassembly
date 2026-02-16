process GATHERSTATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.3':
        'biocontainers/coreutils:9.3' }"

    input:
    tuple val(meta), path(general)

    output:
    path('summary.csv'), emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def coreutils_version = '9.3'
    """
    echo 'sample,loci,kmer_coverage,full_length,probe_length' > summary.csv
    cat *.stats.summary.csv >> summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $coreutils_version
    END_VERSIONS
    """

    stub:
    """
    touch summary.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $coreutils_version
    END_VERSIONS
    """
}
