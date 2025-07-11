/* TODO:
    file endings and prefixes
    check and lint all local modules
    biocontainers/python:3.12 is probably a good general container, already used for ORTHOLOGFILTER
    maybe biocontainers/coreutils:9.3 as well
    publish tblastx as nf-core module
*/

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_targetassembly_pipeline'

include { BLAST_MAKEBLASTDB                            } from '../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN                                 } from '../modules/nf-core/blast/blastn/main'
include { GNU_SORT                                     } from '../modules/local/gnu_sort/main'
include { GNU_SORT as GNU_SORT2                        } from '../modules/local/gnu_sort/main'
include { FASTQTOOLS_SORT                              } from '../modules/local/fastqtools_sort/main'
include { FASTP                                        } from '../modules/nf-core/fastp/main'
include { SPADES                                       } from '../modules/nf-core/spades/main'
include { VSEARCH_CLUSTER                              } from '../modules/nf-core/vsearch/cluster/main'
include { VSEARCH_SORT                                 } from '../modules/nf-core/vsearch/sort/main'
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB2      } from '../modules/nf-core/blast/makeblastdb/main'
include { BLAST_TBLASTX                                } from '../modules/local/blast/tblastx/main'
include { BITSCOREFILTER                               } from '../modules/local/bitscorefilter/main'
include { BEDTOOLS_GETFASTA                            } from '../modules/nf-core/bedtools/getfasta/main'
include { BLAST_TBLASTX as BLAST_TBLASTX2              } from '../modules/local/blast/tblastx/main'
include { ORTHOLOGFILTER                               } from '../modules/local/orthologfilter/main'
include { BEDTOOLS_GETFASTA as ORTHOLOGS_PROBEGETFASTA } from '../modules/nf-core/bedtools/getfasta/main'
include { BEDTOOLS_GETFASTA as ORTHOLOGS_FULLGETFASTA  } from '../modules/nf-core/bedtools/getfasta/main'
include { CLEANHEADERS as CLEANHEADERS_PROBE           } from '../modules/local/cleanheaders/main'
include { CLEANHEADERS as CLEANHEADERS_FULL            } from '../modules/local/cleanheaders/main'
include { QUAST                                        } from '../modules/nf-core/quast/main'

workflow TARGETASSEMBLY {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_reference = file(params.reference)
    ch_probes = file(params.probes)
    ch_versions = Channel.empty()

    /* -------------------------
    Prepare reference and probes
    ------------------------- */

    // make reference genome blast db
    BLAST_MAKEBLASTDB ([[id: ch_reference.baseName], ch_reference])
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    // blastn probes to reference genome
    BLAST_BLASTN ([[id: ch_probes.baseName], ch_probes], BLAST_MAKEBLASTDB.out.db)
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions)

    // TODO: can these two sorts be combined?
    // sort probe/reference hits by bitscore and keep only best probe/reference hit by bitscore
    GNU_SORT (BLAST_BLASTN.out.txt)
    ch_versions = ch_versions.mix(GNU_SORT.out.versions)

    /* -------------
    Sample preparation
    ------------- */

    // TODO: option to skip this?
    // filter adapters, gather sequencing qc with fastp
    FASTP (ch_samplesheet, [], false, false, false)
    ch_versions = ch_versions.mix(FASTP.out.versions)

    // TODO: is this necessary?
    // sort fastqs to avoid errors with SPAdes
    FASTQTOOLS_SORT ( FASTP.out.reads )

    /* -----
    Assembly
    ----- */

    // assemble scaffolds with SPAdes
    SPADES (FASTQTOOLS_SORT.out.fastq.map{[it[0], it[1], [], []]}, [], [])
    ch_versions = ch_versions.mix(SPADES.out.versions)

    // collapse similar scaffolds
    VSEARCH_CLUSTER (SPADES.out.scaffolds)
    ch_versions = ch_versions.mix(VSEARCH_CLUSTER.out.versions)

    // sort scaffolds
    VSEARCH_SORT (VSEARCH_CLUSTER.out.centroids, '--sortbylength')
    ch_versions = ch_versions.mix(VSEARCH_SORT.out.versions)

    /* -------------
    Orthology filter
    ------------- */

    // make blast db from scaffolds
    BLAST_MAKEBLASTDB2 (VSEARCH_SORT.out.fasta)
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB2.out.versions)

    // tblastx probes (query) to scaffolds (db)
    BLAST_TBLASTX (BLAST_MAKEBLASTDB2.out.db.map{[it[0], ch_probes]}, BLAST_MAKEBLASTDB2.out.db)
    ch_versions = ch_versions.mix(BLAST_TBLASTX.out.versions)

    // keep probe/scaffold hits with bit scores at least 80% of top bit score per probe
    // and transform probe/scaffold blast output to bed for bedtools
    BITSCOREFILTER (BLAST_TBLASTX.out.txt)
    ch_versions = ch_versions.mix(BITSCOREFILTER.out.versions)

    // pull regions of scaffold sequences (putative orthologs) that had tblastx probe hits
    together = BITSCOREFILTER.out.bed.join(VSEARCH_SORT.out.fasta)
    BEDTOOLS_GETFASTA ( together.map{it[0..1]}, together.map{it[2]} )
    ch_versions = ch_versions.mix(BEDTOOLS_GETFASTA.out.versions)

    // tblastx putative orthologs (query) to reference genome (db)
    BLAST_TBLASTX2 ( BEDTOOLS_GETFASTA.out.fasta, BLAST_MAKEBLASTDB.out.db )
    ch_versions = ch_versions.mix(BLAST_TBLASTX2.out.versions)

    // keep only top ortholog/reference hit by bitscore for each ortholog
    GNU_SORT2 ( BLAST_TBLASTX2.out.txt )
    ch_versions = ch_versions.mix(GNU_SORT2.out.versions)

    // ortholog filter: make sure putative orthologs intersect the same coordinates as the probe/reference blast
    ORTHOLOGFILTER ( GNU_SORT2.out.sorted, GNU_SORT.out.sorted )
    ch_versions = ch_versions.mix(ORTHOLOGFILTER.out.versions)

    // pull full and probe orthologs
    probe_together = ORTHOLOGFILTER.out.probe.join(VSEARCH_SORT.out.fasta)
    full_together = ORTHOLOGFILTER.out.full.join(VSEARCH_SORT.out.fasta)
    ORTHOLOGS_PROBEGETFASTA ( probe_together.map{it[0..1]}, probe_together.map{it[2]} )
    ORTHOLOGS_FULLGETFASTA ( full_together.map{it[0..1]}, full_together.map{it[2]} )
    ch_versions = ch_versions.mix(ORTHOLOGS_PROBEGETFASTA.out.versions)
    ch_versions = ch_versions.mix(ORTHOLOGS_FULLGETFASTA.out.versions)

    // remove spades information and keep only locus name in fasta headers
    CLEANHEADERS_PROBE ( ORTHOLOGS_PROBEGETFASTA.out.fasta )
    CLEANHEADERS_FULL ( ORTHOLOGS_FULLGETFASTA.out.fasta )
    ch_versions = ch_versions.mix(CLEANHEADERS_PROBE.out.versions)
    ch_versions = ch_versions.mix(CLEANHEADERS_FULL.out.versions)

    // run QUAST on probe orthologs
    QUAST ( CLEANHEADERS_PROBE.out.fasta, [[id: ch_reference.baseName], ch_reference], [[id: ''], []] )

    // Collate and save software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'targetassembly_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}
