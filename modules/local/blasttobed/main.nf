process BLASTTOBED {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(txt)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.key"), emit: key
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #awk 'BEGIN{FS=OFS="\t"} {if(\$9 < \$10) {split(\$1, x, "_"); print x[1] "_" \$2, \$9, \$10} else {split(\$1, y, "_"); print y[1] "_" \$2, \$10, \$9}}' ${txt} > ${prefix}.bed
        
    cut -f 1-2 ${txt} > ${prefix}.key
        
    awk 'BEGIN{FS=OFS="\t"} {if(\$9 < \$10) {split(\$1, x, "_"); print \$2, \$9, \$10, \$1} else {split(\$1, y, "_"); print \$2, \$10, \$9, \$1}}' ${txt} > ${prefix}.bed


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | sed -e '1!d' -e 's/,.*\$//' -e 's/GNU Awk //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        awk: \$(awk --version | sed -e '1!d' -e 's/,.*\$//' -e 's/GNU Awk //')
    END_VERSIONS
    """
}
