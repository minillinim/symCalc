# symCalc

## Overview

Toolbox for calculating number of reads to extract from files when making simulated datasets

## Dependencies

	biom >= 1.3.1
	ParseM >= 0.0.1
	json >= 2.0

## Installation

Download and run:

    python setup.py install [--prefix=<path>]

## Example usage

symCalc does not create reads. It works out how many reads should be generated from each reference sequence for each sample.

symCalc only calculates the number of reads needed to be generated from all sequences in each source file. So if a file contains multiple sequences then it is up to the caller to work out the number of reads to be taken from each sequence based on length ratios.

That being said. symCalc needs to know a few things:

    1. The location of the OTU table used to derive relative abundances.
    2. The names of the files containing reference sequences.
    3. The total number of bases in each reference sequence file.
    4. The total number of reads to be generated for each sample.

Use symCalc like this:

	$ symCalc otu_table file_info sample1_volume [sample2_volume, ...]

otu_table can be in classic or biom format.

file_info is a tab separated file with the following format:

    URL1	TOTAL_BASES1
    URL2	TOTAL_BASES2
    URL3	TOTAL_BASES3
    URL4	TOTAL_BASES4

If you supply a 1 column file then symCalc will work out the total number of bases for you (this may take some time). symCalc will mimic the abundances of the first n OTUs in the otu_table where n is the number of rows in file_info. symCalc will throw an error if there are fewer OTUs in the OTU_FILE than URLs supplied in file_info.

"sampleX_volume" is the total number of reads to be generated from the Xth sample. symCalc will use the first n sites in the otu_table where n is the number of sample_volume entries. symCalc will throw an error if there are fewer sites in the OTU_FILE than sample_volumes supplied on the command line.

symCalc prints all output to stdout.

## Help

If you experience any problems using symCalc, open an [issue](https://github.com/minillinim/symCalc/issues) on GitHub and tell us about it.

## Licence and referencing

Project home page, info on the source tree, documentation, issues and how to contribute, see http://github.com/minillinim/symCalc

This software is currently unpublished

## Copyright

Copyright (c) 2013 Michael Imelfort. See LICENSE.txt for further details.
