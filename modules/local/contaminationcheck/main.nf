process CONTAMINATIONCHECK {
    tag "$meta.id"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
    ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'
    : 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'}"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta)

    output:
    tuple val(meta), path("*.contam_summary.txt"), emit: contam_summary
    tuple val(meta), path("*.contam_index.txt"), emit: contam_index

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo -e "CHROM\\tPOS\\tAD\\tGQ\\tRATIO" > ${meta.id}.txt
    bcftools mpileup \
        --fasta-ref ${fasta} \
        --min-BQ 30 \
        --min-MQ 20 \
        --annotate AD \
        ${bam} | \
    bcftools call \
        --multiallelic-caller \
        --variants-only \
        --annotate GQ | \
    bcftools query \
        --include 'GT="het" && GQ>=20' \
        --format '[%CHROM\\t%POS\\t%AD\\t%GQ\\n]' | \
    awk '
    BEGIN {OFS="\\t"}
    {
        n = split(\$3, a, ",")
        sum = 0
        for (i=1; i<=n; i++) sum += a[i]
        if (sum >= 16) print \$0, a[1]/sum
    }' >> ${meta.id}.contam_summary.txt

    awk '
    BEGIN { OFS = "," }
    NR==1 { next }
    {
        if (\$5 < 0.3) below++
        total++
    }
    END { print "${meta.id}", (total > 0 ? below/total : 0) }' ${meta.id}.contam_summary.txt > ${meta.id}.contam_index.txt
    """
}