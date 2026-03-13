process SPADESFILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(scaffolds)

    output:
    tuple val(meta), path("*.scaffolds_filtered.fa"), emit: scaffolds

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    gunzip ${scaffolds}

    awk '
    /^>/ {
	keep=0;
	match(\$0, /cov_([0-9]+\\.[0-9]+)/, arr);
	match(\$0, /length_([0-9]+)/, len);
	if(arr[1]+0 >= ${params.spades_min_cov} && len[1]+0 >= ${params.spades_min_len}) keep=1
    }
    keep' ${meta.id}.scaffolds.fa > ${meta.id}.scaffolds_filtered.fa

    """

    stub:
    """
    ${create_cmd} ${prefix}.${suffix}
    """
}
