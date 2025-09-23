#!/usr/bin/env python3
import glob

locus_lengths = {}

with open('${probe_reference}') as f:
	loci_reference = []
	locus_lengths['reference'] = {}
	for line in f.readlines():
		line = line.strip()
		if line[0] == '>':
			locus = line
			loci_reference.append(locus)
		else:
			locus_lengths['reference'][locus] = len(line)

#for file in glob.glob('probe_fasta/*'):

with open('locus_lengths.txt', 'w') as f:
	counts = []
	f.write(f'sample,{",".join(loci_reference)}\n')
	for sample, loci in locus_lengths.items():
		for locus in loci_reference:
			counts.append(loci.get[locus])
		f.write(f'sample,{",".join(counts)}\n')
