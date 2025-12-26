#!/usr/bin/env python3

import platform

# input fasta:
# >L1062__tanypteryx_R:105603_P008_WB01_25925147_Calopterygidae_Mnais_andersoni::NODE_36_length_541_cov_29.899177:181-411
# CTGTCGGATACCGCTGTTATGGATATGATGATTTCCAACTTGCAACAACAAAGGCAAGTGACAGAGCAGCTCAGAAGGGAGGCTGCTATAAAACGTATAAGTGTGTCCCAAGCAGTGCA

# output:
# fasta without SPAdes information
# csvs of per locus: kmer coverage, full ortholog length, and probe ortholog length
# single line of averages of above stats in that order

def average(list):
    return sum(list) / len(list)

with open('$fasta') as fasta_input, \
open('$prefix' + '.orthologs.probe.fasta', 'w') as fasta_output, \
open('$prefix' + '.stats.kmer_coverage.csv', 'w') as kmer_coverage_output, \
open('$prefix' + '.stats.length_full.csv', 'w') as full_length_output, \
open('$prefix' + '.stats.length_probe.csv', 'w') as probe_length_output, \
open('$prefix' + '.stats.summary.csv', 'w') as general_output:
    full_lengths, probe_lengths, kmer_coverages = ([] for _ in range(3))
    for line in fasta_input.readlines():
        line = line.strip()
        if line[0] == '>':
            line = line.split(':')
            locus = line[0][1:]
            spades_info = line[3].split('_')
            full_length = int(spades_info[3])
            kmer_coverage = float(spades_info[5])
        else:
            sequence = line
            fasta_output.write('>' + locus + '\\n')
            fasta_output.write(sequence + '\\n')

            probe_length = len(sequence)
            probe_length_output.write(locus + ',' + str(probe_length) + '\\n')
            probe_lengths.append(probe_length)

            full_length_output.write(locus + ',' + str(full_length) + '\\n')
            full_lengths.append(full_length)

            kmer_coverage_output.write(locus + ',' + str(round(kmer_coverage, 2)) + '\\n')
            kmer_coverages.append(kmer_coverage)

    general_output.write('$meta.id,' + ','.join([str(round(average(x), 2)) for x in [kmer_coverages, full_lengths, probe_lengths]]) + '\\n')

with open('versions.yml', 'w') as f:
    f.write("\"$task.process\":\\n")
    f.write(f"    orthologfilter: {platform.python_version()}\\n")