/* TODO:
    BITSCOREFILTER, BLASTTOBED doesn't have container
    probably check and lint all local modules
    biocontainers/python:3.12 is probably a good general container, already used for ORTHOLOGFILTER
    add fastp '-g' to modules.conf (didn't yet because I didn't want pipeline to restart totally over)
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_targetassembly_pipeline'

include { BLAST_MAKEBLASTDB                            } from '../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN                                 } from '../modules/nf-core/blast/blastn/main'
include { GNU_SORT                                     } from '../modules/nf-core/gnu/sort/main'
include { GNU_SORT as GNU_SORT2                        } from '../modules/nf-core/gnu/sort/main'
include { FASTP                                        } from '../modules/nf-core/fastp/main'
include { SPADES                                       } from '../modules/nf-core/spades/main'
include { VSEARCH_CLUSTER                              } from '../modules/nf-core/vsearch/cluster/main'
include { VSEARCH_SORT                                 } from '../modules/nf-core/vsearch/sort/main'
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB2      } from '../modules/nf-core/blast/makeblastdb/main'
include { BLAST_TBLASTX                                } from '../modules/local/blast/tblastx/main'
include { BITSCOREFILTER                               } from '../modules/local/bitscorefilter/main'
include { BLASTTOBED                                   } from '../modules/local/blasttobed/main'
include { BEDTOOLS_GETFASTA                            } from '../modules/nf-core/bedtools/getfasta/main'
include { BLAST_TBLASTX as BLAST_TBLASTX2              } from '../modules/local/blast/tblastx/main'
include { GNU_SORT as GNU_SORT3                        } from '../modules/nf-core/gnu/sort/main'
include { GNU_SORT as GNU_SORT4                        } from '../modules/nf-core/gnu/sort/main'
include { ORTHOLOGFILTER                               } from '../modules/local/orthologfilter/main'
include { BEDTOOLS_GETFASTA as ORTHOLOGS_PROBEGETFASTA } from '../modules/nf-core/bedtools/getfasta/main'
include { BEDTOOLS_GETFASTA as ORTHOLOGS_FULLGETFASTA  } from '../modules/nf-core/bedtools/getfasta/main'
include { CLEANHEADERS as CLEANHEADERS_PROBE           } from '../modules/local/cleanheaders/main'
include { CLEANHEADERS as CLEANHEADERS_FULL            } from '../modules/local/cleanheaders/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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

    // sort probe/reference hits by bitscore
    GNU_SORT (BLAST_BLASTN.out.txt)
    ch_versions = ch_versions.mix(GNU_SORT.out.versions)

    // keep only best probe/reference hit by bitscore
    GNU_SORT2 (GNU_SORT.out.sorted)
    ch_versions = ch_versions.mix(GNU_SORT2.out.versions)

    /* -------------
    Sample preparation
    ------------- */

    // filter adapters, gather sequencing qc with fastp
    FASTP (ch_samplesheet, [], false, false, false)
    ch_versions = ch_versions.mix(FASTP.out.versions)

    /* -----
    Assembly
    ----- */

    // assemble scaffolds with SPAdes
    SPADES (FASTP.out.reads.map{[it[0], it[1], [], []]}, [], [])
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
    BITSCOREFILTER (BLAST_TBLASTX.out.txt)
    ch_versions = ch_versions.mix(BITSCOREFILTER.out.versions)

    // transform probe/scaffold blast output to bed for bedtools
    BLASTTOBED ( BITSCOREFILTER.out.txt )
    ch_versions = ch_versions.mix(BLASTTOBED.out.versions)

    // pull regions of scaffold sequences (putative orthologs) that had tblastx probe hits
    together = BLASTTOBED.out.bed.join(VSEARCH_SORT.out.fasta)
    BEDTOOLS_GETFASTA ( together.map{it[0..1]}, together.map{it[2]} )
    ch_versions = ch_versions.mix(BEDTOOLS_GETFASTA.out.versions)

    // tblastx putative orthologs (query) to reference genome (db)
    BLAST_TBLASTX2 ( BEDTOOLS_GETFASTA.out.fasta, BLAST_MAKEBLASTDB.out.db )
    ch_versions = ch_versions.mix(BLAST_TBLASTX2.out.versions)

    // keep only top ortholog/reference hit by bitscore for each ortholog
    GNU_SORT3 ( BLAST_TBLASTX2.out.txt )
    GNU_SORT4 ( GNU_SORT3.out.sorted )
    ch_versions = ch_versions.mix(GNU_SORT3.out.versions)
    ch_versions = ch_versions.mix(GNU_SORT4.out.versions)

    // ortholog filter: make sure putative orthologs intersect the same coordinates as the probe/reference blast
    ORTHOLOGFILTER ( GNU_SORT4.out.sorted, GNU_SORT2.out.sorted )
    ch_versions = ch_versions.mix(ORTHOLOGFILTER.out.versions)

    // pull full and probe orthologs
    probe_together = ORTHOLOGFILTER.out.probe.join(VSEARCH_SORT.out.fasta)
    full_together = ORTHOLOGFILTER.out.full.join(VSEARCH_SORT.out.fasta)
    ORTHOLOGS_PROBEGETFASTA ( probe_together.map{it[0..1]}, probe_together.map{it[2]} )
    ORTHOLOGS_FULLGETFASTA ( full_together.map{it[0..1]}, full_together.map{it[2]} )
    ch_versions = ch_versions.mix(ORTHOLOGS_PROBEGETFASTA.out.versions)
    ch_versions = ch_versions.mix(ORTHOLOGS_FULLGETFASTA.out.versions)

    // get rid of everything but locus name in fasta headers
    CLEANHEADERS_PROBE ( ORTHOLOGS_PROBEGETFASTA.out.fasta )
    CLEANHEADERS_FULL ( ORTHOLOGS_FULLGETFASTA.out.fasta )
    ch_versions = ch_versions.mix(CLEANHEADERS_PROBE.out.versions)
    ch_versions = ch_versions.mix(CLEANHEADERS_FULL.out.versions)

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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
