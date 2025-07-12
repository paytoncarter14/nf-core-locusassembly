process GETSPADESCOVERAGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'biocontainers/python:3.12'

    input:
    tuple val(meta) , path(fasta)

    output:
    tuple val(meta), path('*.txt'), emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk '
    BEGIN {
        FS = ":"
        print "locus,length,coverage"
    }
    /^>/ {
        name = \$1
        sub(/^>/, "", name)

        if (match(\$0, /length_([0-9]+)/)) {
            length_val = substr(\$0, RSTART + 7, RLENGTH - 7)
        } else {
            length_val = "N/A"
        }

        if (match(\$0, /cov_([0-9.]+)/)) {
            coverage_val = substr(\$0, RSTART + 4, RLENGTH - 4)
        } else {
            coverage_val = "N/A"
        }

        print name "," length_val "," coverage_val
    } ' ${fasta} > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version | head -1 | sed 's|sed (GNU sed) ||g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version | head -1 | sed 's|sed (GNU sed) ||g')
    END_VERSIONS
    """
}
