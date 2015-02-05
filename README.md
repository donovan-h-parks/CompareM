# CompareM

<b>[This project is in active development and not currently recommended for public use.]</b>

CompareM is a comparative genomics toolbox. It provides statistics across groups of genomes (e.g., amino acid identity; core, dispensible, and unique gene sets) and for individual genomes (e.g., GC content, coding density). Emphasis has been placed on providing parallelized implementations for calculating statistics in order to allow scalability to tens of thousands of genomes. The functionality currently planned is:

<i>Comparative genomic statistics:</i>
* calculation of the amino acid identity between genomes
* identification of core, dispensible, and unique gene sets

<i>Single genome statistics:</i>
* GC content
* coding density
* codon and amino acid usage
* automatic identification of translation table
* N50; maximum and mean scaffold/contig size; no. of scaffolds/contigs

A number of auxillary tools are also provides which are often helpful within comparative genomic studies:
* back-translation of amino acid alignments to nucleotides
* identification of homologous genes followed by alignment and tree inference 

## Install

The simplest way to install this package is through pip:
> sudo pip install comparem

This package requires numpy to be installed and makes use of the follow bioinformatic packages:

CheckM relies on several other software packages:

* [prodigal](http://prodigal.ornl.gov/): Hyatt D, Locascio PF, Hauser LJ, Uberbacher EC. 2012. Gene and translation initiation site prediction in metagenomic sequences. <i>Bioinformatics</i> 28: 2223-2230.
* [blast+](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download): Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. 2009. BLAST+: architecture and applications. <i>BMC Bioinformatics</i> 10:421: doi: 10.1186/1471-2105-10-421.
* [pplacer](http://matsen.fhcrc.org/pplacer/): Matsen FA, Kodner RB, Armbrust EV. 2010. pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. BMC Bioinformatics 11: doi:10.1186/1471-2105-11-538.
* [muscle](http://www.drive5.com/muscle/): Edgar RC. 2004. MUSCLE: multiple sequence alignment with high accuracy and throughput. <i>Nucleic Acids Research</i> 32: 1792-1797. doi: 10.1093/nar/gkh340

## Cite

If you find this package useful, please cite this git repository (https://github.com/dparks1134/CompareM)

## Copyright

Copyright © 2014 Donovan Parks. See LICENSE for further details.
