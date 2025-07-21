process GATHERSTATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'biocontainers/python:3.12'

    input:
    tuple val(meta), path(general)

    output:
    path('summary.csv'), emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo 'sample,kmer_coverage,full_length,probe_length' > summary.csv
    cat *.stats.summary.csv >> summary.csv
    """
}
