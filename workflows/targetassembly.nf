/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_targetassembly_pipeline'

include { BLAST_MAKEBLASTDB     } from '../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN          } from '../modules/nf-core/blast/blastn/main'
include { GNU_SORT              } from '../modules/nf-core/gnu/sort/main'
include { GNU_SORT as GNU_SORT2 } from '../modules/nf-core/gnu/sort/main'

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

    //
    // Collate and save software versions
    //
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
