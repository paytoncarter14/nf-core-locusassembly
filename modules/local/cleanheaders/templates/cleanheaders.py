#!/usr/bin/env python3

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
open('$prefix.fasta', 'w') as fasta_output, \
open('$prefix.kmer_coverage.csv', 'w') as kmer_coverage_output, \
open('$prefix.full_length.csv', 'w') as full_length_output, \
open('$prefix.probe_length.csv', 'w') as probe_length_output, \
open('$prefix.general.csv', 'w') as general_output:
    full_lengths, probe_lengths, kmer_coverages = ([] for _ in range(3))
    for line in fasta_input.readlines():
        line = line.strip()
        locus, full_length, probe_length, kmer_coverage = ('',) * 4
        if line[0] == '>':
            line = line.split(':')
            locus = line[0][1:]
            spades_info = line[3].split('_')
            full_length = spades_info[4]
            kmer_coverage = float(spades_info[6])
        else:
            sequence = line
            fasta_output.write('>' + locus + '\\n')
            fasta_output.write(sequence)

            probe_length = len(sequence)
            probe_length_output.write(locus + ',' + str(probe_length) + '\\n')
            probe_lengths.append(probe_length)

            full_length_output.write(locus + ',' + str(full_length) + '\\n')
            full_lengths.append(full_length)

            kmer_coverage_output.write(locus + ',' + kmer_coverage + '\\n')
            kmer_coverages.append(kmer_coverage)

    general_output.write(','.join([str(x) for x in ['$meta', average(kmer_coverages), average(full_lengths), average(probe_lengths)]]))
