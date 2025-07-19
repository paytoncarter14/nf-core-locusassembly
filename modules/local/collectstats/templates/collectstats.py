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

with open('locus_lengths_probe_mqc.csv', 'w') as f:
    f.write('# plot_type: "table"\\n')
    f.write(f'sample,{",".join(loci_reference)}\\n')
    for sample, loci in locus_lengths.items():
        counts = []
        for locus in loci_reference:
            to_append = loci.get(locus)
            if to_append:
                counts.append(to_append['probe'])
            else:
                counts.append('NA')
        f.write(f'{sample},{",".join([str(x) for x in counts])}\\n')

with open('locus_lengths_full_mqc.csv', 'w') as f:
    f.write('# plot_type: "table"\\n')
    f.write(f'sample,{",".join(loci_reference)}\\n')
    for sample, loci in locus_lengths.items():
        counts = []
        for locus in loci_reference:
            to_append = loci.get(locus)
            if to_append:
                counts.append(to_append['full'])
            else:
                counts.append('NA')
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


with open('mean_mapping_coverage_mqc.csv', 'w') as f:
    f.write('# plot_type: "table"\\n')
    f.write('sample,' + ','.join(loci_reference) + '\\n')
    for sample, loci in mapping_coverages.items():
        counts = []
        for locus in loci_reference:
            counts.append(loci.get(locus, 'NA'))
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

with open('pct_at_80_pct_mean_cov_mqc.csv', 'w') as f:
        f.write('# plot_type: "table"\\n')
        f.write('sample,' + ','.join(loci_reference) + '\\n')
        for sample, loci in depths.items():
            counts = []
            for locus in loci_reference:
                if loci.get(locus):
                    counts.append(loci[locus]['above_80'] / loci[locus]['length'])
                else:
                    counts.append('NA')
            f.write(sample + ',' + ','.join([str(x) for x in counts]) + '\\n')

# ----------------- #
# Percent on target #
# ----------------- #

pct_on_targets = {}
for file in glob.glob('fastp_json/*.json'):
        sample_name = os.path.basename(file).split('.')[0]
        with open(file) as f:
                j = json.loads(f.read())
        pct_on_targets[sample_name] = mapped_reads[sample_name] / j['summary']['before_filtering']['total_reads']

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


with open('spades_kmer_coverage_mqc.csv', 'w') as f:
    f.write('# plot_type: "table"\\n')
    f.write('sample,' + ','.join(loci_reference) + '\\n')
    for sample, loci in kmer_coverages.items():
        counts = []
        for locus in loci_reference:
            counts.append(loci.get(locus, 'NA'))
        f.write(sample + ',' + ','.join(counts) + '\\n')

# ------------- #
# General stats #
# ------------- #

with open('summary_mqc.csv', 'w') as out:
    out.write('''# plot_type: "generalstats"
# headers:
#   - Num Loci:
#       min: 0
#   - Avg SPAdes kmer Coverage:
#       min: 0
''')
    out.write('Sample,Num Loci,Avg Probe Length,Avg Full Length,Avg SPAdes kmer Coverage,Avg Mapping Coverage,Avg Pct Bases at 80 Pct Mean Coverage,Pct Reads Mapped\\n')
    for probe_file in glob.glob('probe_fasta/*'):

        sample_name = os.path.basename(probe_file).split('.')[0]

        locus_count = 0
        with open(probe_file) as f:
            for line in f.readlines():
                if line[0] == '>': locus_count += 1

        def mean(inlist):
            return sum([float(x) for x in inlist]) / len(inlist)

        mean_probe_locus_length = round(mean([x['probe'] for x in locus_lengths[sample_name].values()]))
        mean_full_locus_length = round(mean([x['full'] for x in locus_lengths[sample_name].values()]))
        mean_kmer_coverage = round(mean([x for x in kmer_coverages[sample_name].values()]), 2)
        mean_mapping_coverage = round(mean([x for x in mapping_coverages[sample_name].values()]), 2)
        mean_pct_80_pct_mean_coverage = round(mean([x['above_80'] / x['length'] for x in depths[sample_name].values()]), 3)
        pct_on_target = round(pct_on_targets[sample_name], 3)

        output = [str(x) for x in [sample_name, locus_count, mean_probe_locus_length, mean_full_locus_length, mean_kmer_coverage, mean_mapping_coverage, mean_pct_80_pct_mean_coverage, pct_on_target]]
        out.write(','.join(output) + '\\n')
