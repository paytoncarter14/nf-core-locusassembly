process GETLOCUSLENGTH {
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
        print "locus,length"
        current_id = ""
        current_seq_length = 0
    }

    /^>/ {
        if (current_id != "") {
            print current_id "," current_seq_length
        }
        
        current_id = substr(\$0, 2)
        sub(/ .*/, "", current_id)

        current_seq_length = 0
        next
    }

    !/^>/ {
        current_seq_length += length(gensub(/[^A-Za-z]/, "", "g", \$0))
    }

    END {
        if (current_id != "") {
            print current_id "," current_seq_length
        }
    }
    ' ${fasta} > ${prefix}.txt

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
