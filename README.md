# CompareM

<b>[This project is in active development, you are welcomed to use this software though it may be unstable.]</b>

[![version status](https://img.shields.io/pypi/v/comparem.svg)](https://pypi.python.org/pypi/comparem)
[![downloads](https://img.shields.io/pypi/dm/comparem.svg)](https://pypi.python.org/pypi/comparem)

CompareM is a software toolkit which supports performing large-scale comparative genomic analyses. It provides statistics across sets of genomes (e.g., amino acid identity) and for individual genomes (e.g., codon usage). Parallelized implementations are provided for computationally intensive tasks in order to allow scalability to thousands of genomes. Common workflows are provided as single methods to support easy adoption by users, and a more granular interface provided to allow experienced users to exploit specific functionality. CompareM is open source and released under the GNU General Public License (Version 3). 

<i>Comparative genomic statistics:</i>
* average amino acid identity (AAI) between genomes
* taxonomic classification by calculating AAI between query genomes and a reference database

<i>Genomic usage patterns:</i>
* codon usage
* amino acid usage
* kmer usage for k <= 8 (e.g., tetranucleotide)
* stop codon usage

<i>Other:</i>
* di-nucleotide and codon usage patterns for identifying LGT
* visualization of usage patterns with PCA plots and hierarchical clustering dendrograms


## Install

The simplest way to install this package is through pip:
```
> sudo pip install comparem
```

This package uses the numpy and biolib python packages, and requires the follow bioinformatic programs to be on your system path:

* [prodigal](http://prodigal.ornl.gov/) >= 2.6.2: Hyatt D, Locascio PF, Hauser LJ, Uberbacher EC. 2012. Gene and translation initiation site prediction in metagenomic sequences. <i>Bioinformatics</i> 28: 2223-2230.
* [diamond](http://ab.inf.uni-tuebingen.de/software/diamond/) >= 0.7.12: Buchfink B, Xie C, Huson DH. 2015. Fast and sensitive protein alignment using DIAMOND. <i>Nature Methods</i> 12: 59–60 doi:10.1038/nmeth.3176.


## Quick Start

The functionality provided by CompareM can be accessed through the help menu:
```
> comparem -h
```

Usage information about specific functions can also be accessed through the help menu, e.g.:
```
> comparem aa_usage –h
```

## Common Workflow

The most common task performed with CompareM is the calculation of amino acid identities (AAI) between a set of genomes. This can be performed using the <i>aai_wf</i> command:
```
> comparem aai_wf <input_files> <output_dir>
```

The <input_file> argument indicates the set of genomes to compare and can either i) a text file where each line indicating the location of a genome, or ii) a directory containing all genomes to be compared.  The genomic nucleotide sequences of genomes must be in FASTA format. The <output_dir> indicates where the desired directory for all output files. A typical use of this command would be:
```
> comparem --cpus 32 aai_wf my_genomes aai_output
```

where the directory <i>my_genomes</i> contains a set of genomes in FASTA format, the results are to be written to a directory called <i>aai_output</i>, and 32 processors should be used to calculate the results.

A number of optional arguments can also be specified. This includes the sequence similarity parameters used to define reciprocal best hits between genomes(i.e., homologs). By default the e-value (<i>--evalue</i>), percent sequence identity (<i>--per_identity</i>), and percent alignment length (<i>--per_aln_len</i>) parameters are set to 1e-5, 30%, and 70%. When specifying a directory of genomes to process, CompareM only processes files with a <i>fna</i> extension. This can be changes with the <i>--file_ext</i> argument. In addition, if genomes are already represented by amino acid protein sequences (as opposed to genomic nucleotide sequences), this must be specified with the <i>--proteins</i> flag. Otherwise, genes will be identified <i>de novo</i> using the Prodigal gene caller. The time to compute all pairwise AAI values can be substantially reduced by using multiple processors as specified with the <i>--cpus</i> argument. Other arguments are for specialized uses and are discussed in the User's Guide.


## Program Usage

Detailed information regarding the use of CompareM can be found in the User's Guide (user_guide.pdf).


## Cite

If you find this package useful, please cite this git repository (https://github.com/dparks1134/CompareM)


## Copyright

Copyright © 2014 Donovan Parks. See LICENSE for further details.
