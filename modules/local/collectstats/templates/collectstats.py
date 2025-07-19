#!/usr/bin/env python3
import glob, os, json

# ------------- #
# Locus lengths #
# ------------- #

locus_lengths = {}
locus = ''

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
mapped_reads = {}
for file in glob.glob('coverage_txt/*.txt'):
  with open(file) as f:
    next(f)
    sample_name = os.path.basename(file).split('.')[0]
    mapping_coverages[sample_name] = {}
    mapped_reads[sample_name] = 0
    for line in f.readlines():
      line = line.strip().split('\\t')
      locus = line[0]
      num_reads = line[3]
      count = line[6]
      mapping_coverages[sample_name][locus] = count
      mapped_reads[sample_name] += int(num_reads)


with open('mean_mapping_coverage.csv', 'w') as f:
  f.write('sample,' + ','.join(loci_reference) + '\\n')
  for sample, loci in mapping_coverages.items():
    counts = []
    for locus in loci_reference:
      counts.append(loci[locus])
    f.write(sample + ',' + ','.join(counts) + '\\n')

# ---------------------------- #
# Percent at 80% mean coverage #
# ---------------------------- #

depths = {}
for file in glob.glob('depth_txt/*.tsv'):
    with open(file) as f:
        sample_name = os.path.basename(file).split('.')[0]
        depths[sample_name] = {}
        mean_covs = mapping_coverages[sample_name]
        for line in f.readlines():
            line = line.strip().split('\\t')
            locus = line[0]
            if not depths[sample_name].get(locus):
                depths[sample_name][locus] = {'length': 0, 'above_80': 0}
            cov = int(line[2])
            mean_cov = float(mean_covs[locus])
            depths[sample_name][locus]['length'] += 1
            if cov > mean_cov * 0.8:
                depths[sample_name][locus]['above_80'] += 1

with open('pct_at_80_pct_mean_cov.csv', 'w') as f:
    f.write('sample,' + ','.join(loci_reference) + '\\n')
    for sample, loci in depths.items():
      counts = []
      for locus in loci_reference:
        counts.append(loci[locus]['above_80'] / loci[locus]['length'])
      f.write(sample + ',' + ','.join([str(x) for x in counts]) + '\\n')

# ----------------- #
# Percent on target #
# ----------------- #

total_reads = {}
for file in glob.glob('fastp_json/*.json'):
    sample_name = os.path.basename(file).split('.')[0]
    with open(file) as f:
        j = json.loads(f.read())
    total_reads[sample_name] = j['summary']['before_filtering']['total_reads']

with open('pct_on_target.csv', 'w') as f:
    f.write('sample,pct_on_target\\n')
    for sample, total_read in total_reads.items():
      pct_on_target = mapped_reads[sample] / total_read
      f.write(sample + ',' + str(pct_on_target) + '\\n')

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
