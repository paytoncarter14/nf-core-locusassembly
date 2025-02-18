#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--probe_blast', required=True, help='Results of probe sequence blast to reference genome, specified in blast 6 format, one hit per locus')
parser.add_argument('--assembly_blast', required=True, help='Results of assembly scaffold blast to reference genome, specified in blast 6 format')
parser.add_argument('--output_prefix', required=False, help='Prefix that output files will be named with')

args = parser.parse_args()

# qseqid                sseqid          pident  length  mismatch    gapopen qstart  qend    sstart      send        evalue  bitscore
# L100__lampyridae_R    NW_022170509.1  100.000 738     0           0       1       738     21766757    21767494    0.0     1363
with open(args.probe_blast, 'r') as f:
    probes = {}
    for line in f.readlines():
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.strip().split('\t')
        if sstart > send:
            direction = 0
        else:
            direction = 1
        probes[qseqid] = {'start': sstart, 'end': send, 'direction': direction}

with open(args.assembly_blast, 'r') as f:
    scaffolds = {}
    for line in f.readlines():
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.strip().split('\t')
        if sstart > send:
            direction = 0
        else:
            direction = 1 
        scaffolds[qseqid] = {'locus': qseqid.split(':')[0], 'start': sstart, 'end': send, 'direction': direction}

keep = set()
for probe, probe_values in probes.items():
    print()
    print(f'Evaluating locus {probe}')
    scaffolds_eval = {key: value for key, value in scaffolds.items() if value['locus'] == probe}
    print(f'{len(scaffolds_eval)} scaffolds to evaluate')
    for scaffold, scaffold_values in scaffolds_eval.items():
        scaffold_min = min(scaffold_values['start'], scaffold_values['end'])
        scaffold_max = max(scaffold_values['start'], scaffold_values['end'])
        probe_min = min(probe_values['start'], probe_values['end'])
        probe_max = max(probe_values['start'], probe_values['end'])
        print(f'Scaffold: {scaffold_min}-{scaffold_max}')
        print(f'Probe: {probe_min}-{probe_max}')
        # if probe_values['direction'] == scaffold_values['direction']:
        if scaffold_min <= probe_max or scaffold_max >= probe_min:
            keep.add(scaffold)
            print('Same location, will keep')
        else:
            print('Different location')
        # else:
        #     print('Different direction')

if args.output_prefix:
    with open(args.output_prefix + '.ortho_probe.bed', 'w') as ortho_probe:
        with open(args.output_prefix + '.ortho_full.bed', 'w') as ortho_full:
            for x in keep:
                locus, _, scaffold, coords = x.split(':')
                length = scaffold.split('_')[3]
                probe_start, probe_end = coords.split('-')
                ortho_probe.write('\t'.join([scaffold, probe_start, probe_end, f'{locus}:{args.output_prefix}']) + '\n')
                ortho_full.write('\t'.join([scaffold, '0', f'{int(length)-1}', f'{locus}:{args.output_prefix}']) + '\n')