process COLLECTSTATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'biocontainers/python:3.12'

    input:
    tuple val(meta) , path(probe_fasta, stageAs: 'probe_fasta/*')
    tuple val(meta2) , path(full_fasta, stageAs: 'full_fasta/*')
    tuple val(meta3), path(coverage_txt, stageAs: 'coverage_txt/*')
    tuple val(meta4), path(kmer_txt, stageAs: 'kmer_txt/*')
    path(probe_reference)

    output:
    tuple val(meta), path('locus_lengths.csv'), emit: locus_lengths
    tuple val(meta), path('mean_mapping_coverage.csv'), emit: mapping_coverage
    tuple val(meta), path('spades_kmer_coverage.csv'), emit: kmer_coverage

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    template 'collectstats.py'

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
