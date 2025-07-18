#!/usr/bin/env python3
import glob, os

# ------------- #
# Locus lengths #
# ------------- #

locus_lengths = {}

with open('${probe_reference}') as f:
  loci_reference = []
  locus_lengths['reference'] = {}
  for line in f.readlines():
    line = line.strip()
    if line[0] == '>':
      locus = line[1:]
      loci_reference.append(locus)
    else:
      locus_lengths['reference'][locus] = {'probe': len(line), 'full': 'NA'}
  loci_reference = sorted(loci_reference)

for probe_file in glob.glob('probe_fasta/*'):
  sample_name = os.path.basename(probe_file).split('.')[0]
  locus_lengths[sample_name] = {}
  with open(probe_file) as f:
    for line in f.readlines():
      line = line.strip()
      if line[0] == '>':
        locus = line[1:]
      else:
        if not locus_lengths[sample_name].get(locus):
          locus_lengths[sample_name][locus] = {'probe': None, 'full': None}
        locus_lengths[sample_name][locus]['probe'] = len(line)

for full_file in glob.glob('full_fasta/*'):
  sample_name = os.path.basename(full_file).split('.')[0]
  with open(full_file) as f:
    for line in f.readlines():
      line = line.strip()
      if line[0] == '>':
        locus = line[1:]
      else:
        if not locus_lengths[sample_name].get(locus):
          locus_lengths[sample_name][locus] = {'probe': None, 'full': None}
        locus_lengths[sample_name][locus]['full'] = len(line)

print(locus_lengths)

with open('locus_lengths.csv', 'w') as f:
  f.write(f'sample,{",".join(loci_reference)}\\n')
  for sample, loci in locus_lengths.items():
    counts = []
    for locus in loci_reference:
      c = loci.get(locus)
      if c:
        counts.append(f'{c['probe']}/{c['full']}')
    f.write(f'{sample},{",".join([str(x) for x in counts])}\\n')

# --------------------- #
# Mean mapping coverage #
# --------------------- #

mapping_coverages = {}
for file in glob.glob('coverage_txt/*.txt'):
  with open(file) as f:
    next(f)
    sample_name = os.path.basename(file).split('.')[0]
    mapping_coverages[sample_name] = {}
    for line in f.readlines():
      line = line.strip().split('\\t')
      locus = line[0]
      count = line[6]
      mapping_coverages[sample_name][locus] = count


with open('mean_mapping_coverage.csv', 'w') as f:
  f.write('sample,' + ','.join(loci_reference) + '\\n')
  for sample, loci in mapping_coverages.items():
    counts = []
    for locus in loci_reference:
      counts.append(loci[locus])
    f.write(sample + ',' + ','.join(counts) + '\\n')

# -------------------- #
# SPAdes kmer coverage #
# -------------------- #

kmer_coverages = {}
for file in glob.glob('kmer_txt/*.txt'):
  with open(file) as f:
    next(f)
    sample_name = os.path.basename(file).split('.')[0]
    kmer_coverages[sample_name] = {}
    for line in f.readlines():
      line = line.strip().split(',')
      locus = line[0]
      count = line[2]
      kmer_coverages[sample_name][locus] = count


with open('spades_kmer_coverage.csv', 'w') as f:
  f.write('sample,' + ','.join(loci_reference) + '\\n')
  for sample, loci in kmer_coverages.items():
    counts = []
    for locus in loci_reference:
      counts.append(loci[locus])
    f.write(sample + ',' + ','.join(counts) + '\\n')