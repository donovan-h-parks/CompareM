#!/usr/bin/env python

# This is a VERY ROUGH scripts to calculate the core gene set
# between 3 genomes.

# This script should be improved to handle any set of 3 genomes
# and eventually generalized to handle cliques of an arbitrary size.

import os
from collections import defaultdict

genes = defaultdict(list)
primary_genes = set()
for f in os.listdir('.'):
  if f.endswith('.tsv'):
    with open(f) as fin:
      fin.readline()

      for line in fin:
        line_split = line.split('\t')

        gene_a = line_split[0]
        gene_b = line_split[1]

        genes[gene_a].append(gene_b)
        genes[gene_b].append(gene_a)

        genome_id_a = line_split[0].split('~')[1].strip()
        genome_id_b = line_split[1].split('~')[1].strip()
        if genome_id_a == 'bathy.b10.gc_47.genes':
          primary_genes.add(gene_a)
        elif genome_id_b == 'bathy.b10.gc_47.genes':
          primary_genes.add(gene_b)

print len(primary_genes)
core = 0
for gene_id in primary_genes:
  gene_list = genes[gene_id]
  if len(gene_list) == 2:
    if gene_list[1] in genes[gene_list[0]]:
      print '%s\t%s\t%s' % (gene_id, gene_list[0], gene_list[1])
      core += 1

  if len(gene_list) > 2:
    print 'what?'

print core
