

process FASTQTOOLS_SORT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "biocontainers/fastq-tools:0.8.3--h1104d80_5"
    
    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    if ("$fastq" == "${prefix}.fastq.gz" ) error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    zcat ${fastq[0]} | fastq-sort -i | gzip > ${prefix}.R1.fastq.gz
    zcat ${fastq[1]} | fastq-sort -i | gzip > ${prefix}.R2.fastq.gz
    """

    stub:
    """
    touch ${prefix}.fastq.gz
    """
}