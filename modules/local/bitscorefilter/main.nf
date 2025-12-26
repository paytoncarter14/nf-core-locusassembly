process BITSCOREFILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.3':
        'biocontainers/coreutils:9.3' }"

    input:
    tuple val(meta), path(txt)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.key"), emit: key
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${txt}" == "${prefix}.txt") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    sort -k1,1 -k12,12nr ${txt} | awk 'BEGIN{FS=OFS="\t"}{if(\$1!=prev_seqname){prev_seqname=\$1;max_bit_score=\$12}if(\$12>=max_bit_score*0.8){print \$0}}' > ${prefix}.txt
    cut -f 1-2 ${prefix}.txt > ${prefix}.key
    awk 'BEGIN{FS=OFS="\t"} {if(\$9 < \$10) {split(\$1, x, "_"); print \$2, \$9, \$10, \$1} else {split(\$1, y, "_"); print \$2, \$10, \$9, \$1}}' ${prefix}.txt > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(sort --version | sed -n 's/sort (GNU coreutils) \\(.*\\)\$/\\1/p')
        cut: \$(cut --version | sed -n 's/cut (GNU coreutils) \\(.*\\)\$/\\1/p')
        awk: \$(awk --version | sed -e '1!d' -e 's/,.*\$//' -e 's/GNU Awk //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    touch ${prefix}.bed
    touch ${prefix}.key

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        sort: \$(sort --version | sed -n 's/sort (GNU coreutils) \\(.*\\)\$/\\1/p')
        awk: \$(awk --version | sed -e '1!d' -e 's/,.*\$//' -e 's/GNU Awk //')
    END_VERSIONS
    """
}
