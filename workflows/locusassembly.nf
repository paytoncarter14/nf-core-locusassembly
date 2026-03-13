include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_locusassembly_pipeline'

include { BLAST_MAKEBLASTDB                            } from '../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN as BLAST_BLASTNPROBE2REF        } from '../modules/nf-core/blast/blastn/main'
include { BLAST_BLASTN as BLAST_BLASTNPROBE2SCAFFOLD   } from '../modules/nf-core/blast/blastn/main'
include { BLAST_BLASTN as BLAST_BLASTNSCAFFOLD2REF     } from '../modules/nf-core/blast/blastn/main'
include { PARALOGFILTER                                } from '../modules/local/paralogfilter/main'
include { GNU_SORT                                     } from '../modules/local/gnu_sort/main'
include { FASTP                                        } from '../modules/nf-core/fastp/main'
include { SPADES                                       } from '../modules/nf-core/spades/main'
include { SPADESFILTER                                 } from '../modules/local/spadesfilter/main'
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB2      } from '../modules/nf-core/blast/makeblastdb/main'
include { BITSCOREFILTER                               } from '../modules/local/bitscorefilter/main'
include { BEDTOOLS_GETFASTA                            } from '../modules/nf-core/bedtools/getfasta/main'
include { ORTHOLOGFILTER                               } from '../modules/local/orthologfilter/main'
include { BEDTOOLS_GETFASTA as GETORTHOLOGS_PROBE      } from '../modules/nf-core/bedtools/getfasta/main'
include { BEDTOOLS_GETFASTA as GETORTHOLOGS_FULL       } from '../modules/nf-core/bedtools/getfasta/main'
include { CLEANHEADERS_FULL                            } from '../modules/local/cleanheaders_full/main'
include { CLEANHEADERS_PROBE                           } from '../modules/local/cleanheaders_probe/main'
include { MINIMAP2_ALIGN                               } from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_STATS                               } from '../modules/nf-core/samtools/stats/main'
include { GATHERSTATS                                  } from '../modules/local/gatherstats/main'
include { CONTAMINATIONCHECK                           } from '../modules/local/contaminationcheck/main'
include { MULTIQC                                      } from '../modules/nf-core/multiqc/main'

workflow LOCUSASSEMBLY {

    take:
    ch_samplesheet
    
    main:

    ch_reference = file(params.reference)
    ch_probes = file(params.probes)
    ch_versions = channel.empty()

    /* -----------------------------------
    Prepare reference, probes, and samples
    ----------------------------------- */

    // make reference genome blast db
    BLAST_MAKEBLASTDB ([[id: ch_reference.baseName], ch_reference])
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    // blastn probes to reference genome
    BLAST_BLASTNPROBE2REF ([[id: ch_probes.baseName], ch_probes], BLAST_MAKEBLASTDB.out.db)
    ch_versions = ch_versions.mix(BLAST_BLASTNPROBE2REF.out.versions)

    // sort probe/reference hits by bitscore, keep only best probe/reference hit by bitscore,
    // and remove loci with ambiguous hits
    PARALOGFILTER (BLAST_BLASTNPROBE2REF.out.txt)
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
    SPADESFILTER ( SPADES.out.scaffolds )

    /* -------------
    Orthology filter
    ------------- */

    // make blast db from scaffolds
    BLAST_MAKEBLASTDB2 (SPADESFILTER.out.scaffolds)
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB2.out.versions)

    // tblastx probes (query) to scaffolds (db)
    BLAST_BLASTNPROBE2SCAFFOLD (BLAST_MAKEBLASTDB2.out.db.map{[it[0], ch_probes]}, BLAST_MAKEBLASTDB2.out.db)
    ch_versions = ch_versions.mix(BLAST_BLASTNPROBE2SCAFFOLD.out.versions)

    // keep probe/scaffold hits with bit scores at least 80% of top bit score per probe
    // and transform probe/scaffold blast output to bed for bedtools
    BITSCOREFILTER (BLAST_BLASTNPROBE2SCAFFOLD.out.txt)
    ch_versions = ch_versions.mix(BITSCOREFILTER.out.versions)

    // pull regions of scaffold sequences (putative orthologs) that had tblastx probe hits
    together = BITSCOREFILTER.out.bed.join(SPADESFILTER.out.scaffolds)
    BEDTOOLS_GETFASTA ( together.map{it[0..1]}, together.map{it[2]} )
    ch_versions = ch_versions.mix(BEDTOOLS_GETFASTA.out.versions)

    // tblastx putative orthologs (query) to reference genome (db)
    BLAST_BLASTNSCAFFOLD2REF ( BEDTOOLS_GETFASTA.out.fasta, BLAST_MAKEBLASTDB.out.db )
    ch_versions = ch_versions.mix(BLAST_BLASTNSCAFFOLD2REF.out.versions)

    // keep only top ortholog/reference hit by bitscore for each ortholog
    GNU_SORT ( BLAST_BLASTNSCAFFOLD2REF.out.txt )
    ch_versions = ch_versions.mix(GNU_SORT.out.versions)

    // make sure putative orthologs intersect the same coordinates as the probe/reference blast
    ORTHOLOGFILTER ( GNU_SORT.out.sorted, PARALOGFILTER.out.sorted )
    ch_versions = ch_versions.mix(ORTHOLOGFILTER.out.versions)

    // pull full and probe orthologs
    full_input = ORTHOLOGFILTER.out.full.join(SPADESFILTER.out.scaffolds)
    probe_input = ORTHOLOGFILTER.out.probe.join(SPADESFILTER.out.scaffolds)
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
    SAMTOOLS_STATS ( MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index), [[id: ch_reference.baseName], ch_reference] )

    // check for contamination
    CONTAMINATIONCHECK ( MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index).join(CLEANHEADERS_FULL.out.fasta) )

    // collect stats for all samples into summary.csv
    GATHERSTATS ( CLEANHEADERS_PROBE.out.general.collect().map{[[id: 'all_samples'], it]}, CONTAMINATIONCHECK.out.contam_index.map{it[1]}.collect().map{[[id: 'all_samples'], it]} )

    // send all stats to MultiQC
    multiqc_input = FASTP.out.json.map{it[1]}.mix(SAMTOOLS_STATS.out.stats.map{it[1]}, GATHERSTATS.out.summary).collect().map{[[id: 'all_samples'], it, [], [], [], []]}
    MULTIQC ( multiqc_input )

    // Collate and save software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'locusassembly_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}
