process CLEANHEADERS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'biocontainers/python:3.12'

    input:
    tuple val(meta) , path(fasta)

    output:
    tuple val(meta), path('*.fasta'), emit: fasta
    path('*.stats.kmer_coverage.csv'), emit: kmer_coverage, optional: true
    path('*.stats.length_full.csv'), emit: full_length, optional: true
    path('*.stats.length_probe.csv'), emit: probe_length, optional: true
    path('*.stats.summary.csv'), emit: general, optional: true
    // path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${fasta}" == "${prefix}.fasta") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    if ( fasta.name.tokenize('.')[-2] == 'probe_ortho' ) template "cleanheaders.py"
    else
    """
    sed 's|:.*\$||g' ${fasta} > ${prefix}.orthologs.full.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version | head -1 | sed 's|sed (GNU sed) ||g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ( fasta.name.tokenize('.')[-2] == 'probe_ortho' )
    """
    touch ${prefix}.orthologs.probe.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version | head -1 | sed 's|sed (GNU sed) ||g')
    END_VERSIONS
    """
    else
    """
    touch ${prefix}.orthologs.full.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version | head -1 | sed 's|sed (GNU sed) ||g')
    END_VERSIONS
    """
}
