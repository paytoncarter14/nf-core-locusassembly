process CLEANHEADERS_PROBE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12':
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta) , path(fasta)

    output:
    tuple val(meta), path('*.orthologs.probe.fasta'), emit: fasta
    path('*.stats.kmer_coverage.csv'), emit: kmer_coverage, optional: true
    path('*.stats.length_full.csv'), emit: full_length, optional: true
    path('*.stats.length_probe.csv'), emit: probe_length, optional: true
    path('*.stats.summary.csv'), emit: general, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def python_version = '3.12.2'
    if ("${fasta}" == "${prefix}.fasta") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    template "cleanheaders.py"

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ( fasta.name.tokenize('.')[-2] == 'probe_ortho' )
    """
    touch ${prefix}.orthologs.probe.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${python_version}
    END_VERSIONS
    """
    else
    """
    touch ${prefix}.orthologs.probe.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: ${python_version}
    END_VERSIONS
    """
}
