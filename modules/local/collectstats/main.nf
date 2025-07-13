process COLLECTSTATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'biocontainers/python:3.12'

    input:
    tuple val(meta) , path(probe_fasta, stageAs: 'probe_fasta/*')
    tuple val(meta2) , path(full_fasta, stageAs: 'full_fasta/*')
    tuple val(meta3), path(coverage_txt, stageAs: 'coverage_txt/*')
    path(probe_reference)

    output:
    tuple val(meta), path('locus_lengths.csv'), emit: locus_lengths

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
#!/usr/bin/env python3
import glob, os

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
      print(c)
      if c:
        counts.append(f'{c['probe']}/{c['full']}')
    f.write(f'{sample},{",".join([str(x) for x in counts])}\\n')
    print(sample)
    print(loci)
    print(len(loci))
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version | head -1 | sed 's|sed (GNU sed) ||g')
    END_VERSIONS
    """
}
