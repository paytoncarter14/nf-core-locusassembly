include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_targetassembly_pipeline'

include { BLAST_MAKEBLASTDB                            } from '../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN as BLAST_BLASTNPROBE2REF        } from '../modules/nf-core/blast/blastn/main'
include { BLAST_BLASTN as BLAST_BLASTNPROBE2SCAFFOLD   } from '../modules/nf-core/blast/blastn/main'
include { BLAST_BLASTN as BLAST_BLASTNSCAFFOLD2REF     } from '../modules/nf-core/blast/blastn/main'
include { BLAST_TBLASTX as BLAST_TBLASTXPROBE2REF      } from '../modules/local/blast/tblastx/main'
include { PARALOGFILTER                                } from '../modules/local/paralogfilter/main'
include { GNU_SORT as GNU_SORT2                        } from '../modules/local/gnu_sort/main'
include { FASTP                                        } from '../modules/nf-core/fastp/main'
include { SPADES                                       } from '../modules/nf-core/spades/main'
include { GAWK                                         } from '../modules/nf-core/gawk/main'
include { VSEARCH_CLUSTER                              } from '../modules/nf-core/vsearch/cluster/main'
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB2      } from '../modules/nf-core/blast/makeblastdb/main'
include { BLAST_TBLASTX                                } from '../modules/local/blast/tblastx/main'
include { BITSCOREFILTER                               } from '../modules/local/bitscorefilter/main'
include { BEDTOOLS_GETFASTA                            } from '../modules/nf-core/bedtools/getfasta/main'
include { GUNZIP                                       } from '../modules/nf-core/gunzip/main'
include { BLAST_TBLASTX as BLAST_TBLASTX2              } from '../modules/local/blast/tblastx/main'
include { ORTHOLOGFILTER                               } from '../modules/local/orthologfilter/main'
include { BEDTOOLS_GETFASTA as GETORTHOLOGS_PROBE      } from '../modules/nf-core/bedtools/getfasta/main'
include { BEDTOOLS_GETFASTA as GETORTHOLOGS_FULL       } from '../modules/nf-core/bedtools/getfasta/main'
include { CLEANHEADERS_FULL                            } from '../modules/local/cleanheaders_full/main'
include { CLEANHEADERS_PROBE                           } from '../modules/local/cleanheaders_probe/main'
include { MINIMAP2_ALIGN                               } from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_COVERAGE                            } from '../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_DEPTH                               } from '../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_STATS                               } from '../modules/nf-core/samtools/stats/main'
include { BCFTOOLS_MPILEUP                             } from '../modules/nf-core/bcftools/mpileup/main'
include { GATHERSTATS                                  } from '../modules/local/gatherstats/main'

workflow TARGETASSEMBLY {

    take:
    ch_samplesheet
    
    main:

    ch_reference = file(params.reference)
    ch_probes = file(params.probes)
    ch_versions = Channel.empty()

    /* -----------------------------------
    Prepare reference, probes, and samples
    ----------------------------------- */

    // make reference genome blast db
    BLAST_MAKEBLASTDB ([[id: ch_reference.baseName], ch_reference])
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    // blastn probes to reference genome
    if (params.probe2ref_blast_method == 'blastn') {
        BLAST_BLASTNPROBE2REF ([[id: ch_probes.baseName], ch_probes], BLAST_MAKEBLASTDB.out.db)
        ch_versions = ch_versions.mix(BLAST_BLASTNPROBE2REF.out.versions)
        probe2ref = BLAST_BLASTNPROBE2REF
    } else if (params.probe2ref_blast_method == 'tblastx') {
        BLAST_TBLASTXPROBE2REF ([[id: ch_probes.baseName], ch_probes], BLAST_MAKEBLASTDB.out.db)
        ch_versions = ch_versions.mix(BLAST_TBLASTXPROBE2REF.out.versions)
        probe2ref = BLAST_TBLASTXPROBE2REF
    }

    // sort probe/reference hits by bitscore and keep only best probe/reference hit by bitscore
    // FEATURE: add filter for paralogy, low quality hits
    PARALOGFILTER (probe2ref.out.txt, params.paralog_ratio_threshold)
    ch_versions = ch_versions.mix(PARALOGFILTER.out.versions)

    // filter reads, gather sequencing qc with fastp
    FASTP (ch_samplesheet, [], false, false, false)
    ch_versions = ch_versions.mix(FASTP.out.versions)

    /* -----
    Assembly
    ----- */

    // assemble scaffolds with SPAdes
    SPADES (FASTP.out.reads.map{[it[0], it[1], [], []]}, [], [])
    ch_versions = ch_versions.mix(SPADES.out.versions)

    // filter scaffolds by minimum kmer coverage and length
    GAWK (SPADES.out.scaffolds, [], false)
    // ch_versions = ch_versions.mix(GAWK.out.versions_gawk)

    // collapse similar scaffolds and unzip output
    // VSEARCH_CLUSTER (GAWK.out.output)
    // ch_versions = ch_versions.mix(VSEARCH_CLUSTER.out.versions)

    // GUNZIP (VSEARCH_CLUSTER.out.centroids)
    // ch_versions = ch_versions.mix(GUNZIP.out.versions)

    /* -------------
    Orthology filter
    ------------- */

    // make blast db from scaffolds
    BLAST_MAKEBLASTDB2 (GAWK.out.output)
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB2.out.versions)

    // tblastx probes (query) to scaffolds (db)
    BLAST_BLASTNPROBE2SCAFFOLD (BLAST_MAKEBLASTDB2.out.db.map{[it[0], ch_probes]}, BLAST_MAKEBLASTDB2.out.db)
    ch_versions = ch_versions.mix(BLAST_BLASTNPROBE2SCAFFOLD.out.versions)

    // keep probe/scaffold hits with bit scores at least 80% of top bit score per probe
    // and transform probe/scaffold blast output to bed for bedtools
    BITSCOREFILTER (BLAST_BLASTNPROBE2SCAFFOLD.out.txt)
    ch_versions = ch_versions.mix(BITSCOREFILTER.out.versions)

    // pull regions of scaffold sequences (putative orthologs) that had tblastx probe hits
    together = BITSCOREFILTER.out.bed.join(GUNZIP.out.gunzip)
    BEDTOOLS_GETFASTA ( together.map{it[0..1]}, together.map{it[2]} )
    ch_versions = ch_versions.mix(BEDTOOLS_GETFASTA.out.versions)

    // tblastx putative orthologs (query) to reference genome (db)
    BLAST_BLASTNSCAFFOLD2REF ( BEDTOOLS_GETFASTA.out.fasta, BLAST_MAKEBLASTDB.out.db )
    ch_versions = ch_versions.mix(BLAST_BLASTNSCAFFOLD2REF.out.versions)

    // keep only top ortholog/reference hit by bitscore for each ortholog
    GNU_SORT2 ( BLAST_BLASTNSCAFFOLD2REF.out.txt )
    ch_versions = ch_versions.mix(GNU_SORT2.out.versions)

    // make sure putative orthologs intersect the same coordinates as the probe/reference blast
    ORTHOLOGFILTER ( GNU_SORT2.out.sorted, PARALOGFILTER.out.sorted )
    ch_versions = ch_versions.mix(ORTHOLOGFILTER.out.versions)

    // pull full and probe orthologs
    full_input = ORTHOLOGFILTER.out.full.join(GUNZIP.out.gunzip)
    probe_input = ORTHOLOGFILTER.out.probe.join(GUNZIP.out.gunzip)
    GETORTHOLOGS_FULL ( full_input.map{it[0..1]}, full_input.map{it[2]} )
    GETORTHOLOGS_PROBE ( probe_input.map{it[0..1]}, probe_input.map{it[2]} )
    ch_versions = ch_versions.mix(GETORTHOLOGS_FULL.out.versions)
    ch_versions = ch_versions.mix(GETORTHOLOGS_PROBE.out.versions)

    /* -----------
    QC and cleanup
    ----------- */

    // remove spades information and keep only locus name in fasta headers
    // calculate some stats as well
    CLEANHEADERS_FULL ( GETORTHOLOGS_FULL.out.fasta )
    ch_versions = ch_versions.mix(CLEANHEADERS_FULL.out.versions)

    CLEANHEADERS_PROBE ( GETORTHOLOGS_PROBE.out.fasta )
    ch_versions = ch_versions.mix(CLEANHEADERS_PROBE.out.versions)

    // align fastqs with orthologs to get coverage stats
    minimap2_input = FASTP.out.reads.join(CLEANHEADERS_FULL.out.fasta)
    MINIMAP2_ALIGN ( minimap2_input.map{[it[0], it[1]]}, minimap2_input.map{[it[0], it[2]]}, true, 'bai', true, false )

    // get alignment coverage and depth stats
    samtools_input = MINIMAP2_ALIGN.out.bam.join(CLEANHEADERS_FULL.out.fasta)
    SAMTOOLS_COVERAGE ( samtools_input.map{[it[0], it[1], []]}, samtools_input.map{[it[0], it[2]]}, [[], []])
    SAMTOOLS_DEPTH ( MINIMAP2_ALIGN.out.bam, [[], []] )
    SAMTOOLS_STATS ( MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index), [[id: ch_reference.baseName], ch_reference] )

    // collect stats for all samples into summary.csv
    GATHERSTATS ( CLEANHEADERS_PROBE.out.general.collect().map{[[id: 'all_samples'], it]} )

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

/* archived processes
include { VSEARCH_SORT                                 } from '../modules/nf-core/vsearch/sort/main'
include { FASTQTOOLS_SORT                              } from '../modules/local/fastqtools_sort/main'

include { SAMTOOLS_FAIDX                               } from '../modules/nf-core/samtools/faidx/main'
include { GETLOCUSLENGTH as GETLOCUSLENGTH_PROBE       } from '../modules/local/getlocuslength/main'
include { GETLOCUSLENGTH as GETLOCUSLENGTH_FULL        } from '../modules/local/getlocuslength/main'
include { COLLECTSTATS                                 } from '../modules/local/collectstats/main'

include { GETSPADESCOVERAGE                            } from '../modules/local/getspadescoverage/main'

// sort scaffolds
VSEARCH_SORT (VSEARCH_CLUSTER.out.centroids, '--sortbylength')
ch_versions = ch_versions.mix(VSEARCH_SORT.out.versions)

// sort fastqs to avoid errors with SPAdes
FASTQTOOLS_SORT ( FASTP.out.reads )

// align fastqs with orthologs to get coverage stats
minimap2_input = FASTQTOOLS_SORT.out.fastq.join(CLEANHEADERS_PROBE.out.fasta)
MINIMAP2_ALIGN ( minimap2_input.map{[it[0], it[1]]}, minimap2_input.map{[it[0], it[2]]}, true, 'bai', true, false )

// get alignment coverage and depth stats
samtools_input = MINIMAP2_ALIGN.out.bam.join(CLEANHEADERS_PROBE.out.fasta)
SAMTOOLS_COVERAGE ( samtools_input.map{[it[0], it[1], []]}, samtools_input.map{[it[0], it[2]]}, [[], []])
SAMTOOLS_DEPTH ( MINIMAP2_ALIGN.out.bam, [[], []] )

// collect statistics
// locus lengths
// mapping coverage
// spades kmer coverage
COLLECTSTATS (
    CLEANHEADERS_PROBE.out.fasta.map{it[1]}.collect().map{[[id: 'all_samples'], it]},
    CLEANHEADERS_FULL.out.fasta.map{it[1]}.collect().map{[[id: 'all_samples'], it]},
    SAMTOOLS_COVERAGE.out.coverage.map{it[1]}.collect().map{[[id: 'all_samples'], it]},
    GETSPADESCOVERAGE.out.txt.map{it[1]}.collect().map{[[id: 'all_samples'], it]},
    SAMTOOLS_DEPTH.out.tsv.map{it[1]}.collect().map{[[id: 'all_samples'], it]},
    FASTP.out.json.map{it[1]}.collect().map{[[id: 'all_samples'], it]},
    ch_probes
)
// pull spades coverage information from the orthologs
GETSPADESCOVERAGE ( ORTHOLOGS_PROBEGETFASTA.out.fasta )
*/
