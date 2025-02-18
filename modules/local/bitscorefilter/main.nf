process BITSCOREFILTER {
    tag "$meta.id"
    label 'process_single'

    // conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/python:3.12':
    //     'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(txt)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${txt}" == "${prefix}.txt") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    sort -k1,1 -k12,12nr ${txt} | awk 'BEGIN{FS=OFS="\t"}{if(\$1!=prev_seqname){prev_seqname=\$1;max_bit_score=\$12}if(\$12>=max_bit_score*0.8){print \$0}}' > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(sort --version | sed -n 's/sort (GNU coreutils) \\(.*\\)\$/\\1/p')
        awk: \$(awk --version | sed -e '1!d' -e 's/,.*\$//' -e 's/GNU Awk //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        sort: \$(sort --version | sed -n 's/sort (GNU coreutils) \\(.*\\)\$/\\1/p')
        awk: \$(awk --version | sed -e '1!d' -e 's/,.*\$//' -e 's/GNU Awk //')
    END_VERSIONS
    """
}
