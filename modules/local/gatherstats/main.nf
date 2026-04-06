process GATHERSTATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:9.3':
        'biocontainers/coreutils:9.3' }"

    input:
    tuple val(meta), path(general)

    output:
    path('summary_mqc.csv'), emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def coreutils_version = '9.3'
    """
    echo '# plot_type: "generalstats"' > summary_mqc.csv
    cat *.stats.summary.csv | sort -k1,1 > summary-sorted.csv
    if [ -e *.contam_index.txt ]; then
        echo 'sample,loci,kmer_coverage,full_length,probe_length,contam_index' >> summary_mqc.csv
        cat *.contam_index.txt | sort -k1,1 > contam-index-sorted.csv
        join -t ',' summary-sorted.csv contam-index-sorted.csv >> summary_mqc.csv
    else
        echo 'sample,loci,kmer_coverage,full_length,probe_length' >> summary_mqc.csv
        cat summary-sorted.csv >> summary_mqc.csv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $coreutils_version
    END_VERSIONS
    """

    stub:
    """
    touch summary_mqc.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $coreutils_version
    END_VERSIONS
    """
}
