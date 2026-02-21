process PARALOGFILTER {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.3':
        'biocontainers/coreutils:9.3' }"

    input:
    tuple val(meta), path(input)
    val(ratio_threshold)

    output:
    tuple val(meta), file( "${output_file}" )   , emit: sorted
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args     ?: ''
    def args2       = task.ext.args2    ?: ''
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    suffix          = task.ext.suffix   ?: "${input.extension}"
    output_file     = "${prefix}.${suffix}"
    def VERSION     = "9.3" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    if ("$input" == "$output_file") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
# Sort by locus and bitscore (descending),
# then use awk to filter any loci where the second best hit
# is more than [ratio_threshold, e.g. 0.9] times as good as the best hit (i.e. it's a paralog)
sort -k1,1 -k12,12nr ${input} | awk '

BEGIN {
    count = 1
    locus = "init"
    prev_locus = ""
    best_score = ""
    is_bad = "true"
}

{
    locus = \$1
    score = \$12
    if (locus != prev_locus) {
        best_score = score
        prev_locus = locus
        if (is_bad == "false") {
            print best_line
        }
        is_bad = "false"
        best_line = \$0
        count = 1
    } else {
        count++
        if (count == 2) {
            if (score / best_score >= ${ratio_threshold}) {
                is_bad = "true"
            }
        }
    }
}

END {
    if (count == 1) {
        print \$0
    }
}

' > ${output_file}

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    coreutils: $VERSION
END_VERSIONS
    """

    stub:
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    suffix          = task.ext.suffix   ?: "${input.extension}"
    output_file     = "${prefix}.${suffix}"
    def VERSION     = "9.3"

    if ("$input" == "$output_file") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${output_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """
}
