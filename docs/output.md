# nf-core/targetassembly: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

### fastp

Since SPAdes has its own error correction, fastp is configured by default in this pipeline to only perform adapter trimming. It will also generate QC metrics (e.g. number of reads, %>Q30, GC content).

<details markdown="1">

  <summary>Output files</summary>

- `fastp/`
  - `*_1.fastp.fastq.gz`: R1 of fastp-filtered input FASTQ.
  - `*_2.fastp.fastq.gz`: R2 of fastp-filtered input FASTQ.
  - `*.fastp.html`: fastp HTML report.
  - `*.fastp.json`: fastp JSON results.
  - `*.fastp.log`: fastp log.

</details>

### SPAdes

SPAdes is a versatile de novo genome assembler. Originally intended for small genomes,  it also works exceptionally well with target enrichment data.

<details markdown="1">

  <summary>Output files</summary>

- `spades/`
  - `*.assembly.gfa.gz`: SPAdes assembly graph and scaffolds paths.
  - `*.contigs.fa.gz`: SPAdes contigs.
  - `*.scaffolds.fa.gz`: SPAdes scaffolds.
  - `*.spades.log`: SPAdes log.

</details>

### SPAdes filter

This process removes SPAdes scaffolds that don't meet a certain length and kmer coverage requirement.

<details markdown="1">

  <summary>Output files</summary>

- `spadesfilter/`
  - `*.scaffolds_filtered.fa`: SPAdes scaffolds filtered to remove those shorter than the `spades_min_len` parameter (default = 100) and with lower kmer coverage than the `spades_min_cov` parameter (default = 1).

</details>

### BLAST

The pipeline uses NCBI BLAST to perform several searches, all with the dc-megablast task. These searches include:

- SPAdes scaffolds to probe target sequences: used to initially assign putative orthology.
- Probe target sequences to reference genome: used to ensure probe target sequences aren't paralogs, and used to confirm orthology with the scaffold to reference search.
- SPAdes scaffolds to reference genome: used with the probe to reference search to ensure putative orthologs and their respective probe sequences hit the same place on the reference genome.

<details markdown="1">

  <summary>Output files</summary>

- `blast/`
  - `*.probe2scaffold.txt`: BLAST results of probe target sequences to the SPAdes scaffolds.
  - `*.scaffold2reference.txt`: BLAST results of SPAdes scaffolds to the reference genome.
  - `*.probe2reference.txt`: BLAST results of the probe target sequences to the reference genome.
  - One folder per sample ID and one for the reference genome with the generated BLAST databases.

</details>

### paralogfilter

This process filters paralogs by removing any loci where the probe target sequences have ambiguous hits on the reference genome. If a probe target sequence has a secondary hit in a different location with a similar bitscore to the best hit (within a ratio set by the `paralogfilter_ratio_threshold` parameter, default = 0.9), that locus is discarded.

<details markdown="1">

  <summary>Output files</summary>

- `paralogfilter/`:
  - `*.probes_filtered.txt`: The probe to reference BLAST results with any paralogs removed.

</details>

### bitscorefilter

This process is used to filter probe to scaffold BLAST results. It only keep hits with bitscores at least 80% of the top bitscore per probe.

<details markdown="1">

  <summary>Output files</summary>

- `bitscorefilter/`
  - `*.filtered.bed`: Coordinates of BLAST probe hits on scaffolds.
  - `*.filtered.key`: Dictionary of locus names to associated scaffolds.
  - `*.filtered.txt`: BLAST results of probes to scaffolds, filtered to only include hits with bitscores at least 80% of top bitscore per probe.

</details>

### bedtools

This tool is used to pull regions from the SPAdes scaffolds that had probe to scaffold BLAST hits (putative orthologs). It is also used to pull the final orthologs from the scaffolds after the ortholog filter.

<details markdown="1">

  <summary>Output files</summary>

- `bedtools/`
  - `*.fa`: FASTA file with regions of scaffolds (putative orthologs) that had BLAST probe hits.
- `getorthologs/`
  - `*.full_ortho.fa`: The full scaffolds of all orthologs that passed the orthology filter.
  - `*.probe_ortho.fa`: The scaffolds of all orthologs that passed the orthology filter, trimmed to include only the probe to scaffold BLAST hit coordinates.

</details>

### GNU sort

This tool is used to keep only the top ortholog to reference hit by bitscore for each ortholog.

<details markdown="1">

  <summary>Output files</summary>

- `gnu/`
  - `*.txt`: The scaffold to reference BLAST results, filtered to only keep the top ortholog/reference hit by bitscore for each ortholog.

</details>

### orthologfilter

This process ensures that BLAST searches for the putative orthologs and their associated probes hit the same location on the reference genome.

<details markdown="1">

  <summary>Output files</summary>

- `orthologfilter/`
  - `*.ortho_full.bed`: BED file with all loci that pass the ortholog filter.
  - `*.ortho_probe.bed`: BED file with all loci that pass the ortholog filter and coordinates that only include the probe to scaffold BLAST hit coordinates.

</details>

### cleanheaders

This process removes SPAdes metrics from the final ortholog FASTA files so only the sample ID remains. The extracted metrics are saved in separate files.

<details markdown="1">

  <summary>Output files</summary>

- `cleanheaders/`
  - `*.orthologs.full.fasta`: The final orthologs FASTA file, providing the entire SPAdes scaffold, with SPAdes information removed from the FASTA headers.
  - `*.orthologs.probe.fasta`: The final orthologs FASTA file, with scaffolds trimmed to include only the probe to scaffold BLAST hit coordinates and with SPAdes information removed from the FASTA headers.
  - `*.stats.kmer_coverage.csv`: CSV with SPAdes kmer coverage values for each locus.
  - `*.stats.length_full.csv`: CSV with lengths of each full scaffold for each locus.
  - `*.stats.length_probe.csv`: CSV with lengths of the probe coordinate trimmed orthologs for each locus.
  - `*.stats.summary.csv`: CSV with number of loci recovered, mean SPAdes kmer coverage of all loci, mean full length of each locus, and mean probe length of each locus.

</details>

### gatherstats

This process concatenates all `*.stats.summary.csv` files from the cleanheaders process.

<details markdown="1">

  <summary>Output files</summary>

- `gatherstats/`
  - `summary_mqc.csv`: CSV of all `cleanheaders/*.stats.summary.csv` results concatenated.

</details>

### minimap2

The sequencing reads are mapped to the resulting full-length orthologs. This is done to calculate a contamination index, as well as provide quality metrics with samtools stats.

<details markdown="1">

  <summary>Output files</summary>

- `minimap2/`
  - `*.bam`: The BAM result of mapping the reads back to the full-length scaffolds.
  - `*.bam.bai`: Companion BAM index file.

</details>

### samtools

The pipeline uses samtools stats to generate stats relevant to target enrichment sequencing, such as percent mapped/on-target.

<details markdown="1">

  <summary>Output files</summary>

- `samtools/`
  - `*.stats`: Results of samtools stats.

</details>

### contaminationcheck

This process uses BCFtools to call variants on the mapped reads, which are then filtered to only include high-quality heterozygous variants. For diploid organisms, heterozygous variants should have allele frequencies of 0.5, and substantial deviations could reflect certain kinds of contamination. The contamination index is calculated as the proportion of heterozygous calls with alternative allele frequencies below the `contaminationcheck_af_threshold` parameter (default = 0.3).

<details markdown="1">

  <summary>Output files</summary>

- `contaminationcheck/`
  - `*.contam_index.txt`: CSV with the calculated contamination index, i.e.  proportion of heterozygous calls with alternative allele frequencies below the `contaminationcheck_af_threshold` parameter (default = 0.3).
  - `*.contam_summary.txt`: Tab-separated text file with heterozygous calls that pass all filtering parameters.

</details>

### MultiQC

MultiQC generates a single HTML report with interactive tables and figures. It compiles metrics from fastp, SPAdes, samtools stats, and the contamination check.

<details markdown="1">

  <summary>Output files</summary>

- `multiqc/`
  - `multiqc_data/`: Data for MultiQC report.
  - `multiqc_report.html`: MultiQC HTML report.

</details>

### Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

<details markdown="1">

  <summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>