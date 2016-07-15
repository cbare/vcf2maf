vcf2maf: a VCF to MAF conversion utility
========================================

**This code is no longer maintained**. Please have a look at: [mskcc's vcf2maf](https://github.com/mskcc/vcf2maf).

Files in the MAF (Mutation Annotation Format) are used in TCGA to track DNA variants/mutations.

Vcf2maf is a tool for converting files in [Variant Call Format (VCF)](https://wiki.nci.nih.gov/display/TCGA/TCGA+Variant+Call+Format+%28VCF%29+1.1+Specification) to [MAF](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format) format.

Use the [DCC Archive Validator](https://wiki.nci.nih.gov/display/TCGA/TCGA+Archive+Validator) to check the integrity of a MAF file.

Usage
-----

Basic usage of the vcf2maf converter takes a vcf file as input and produces a corresponding maf file:

    python vcf2maf.py -v --maf example1.maf examples/example1.vcf

List the samples in a VCF file, like this:

    python vcf2maf.py -v --list-samples examples/example1.vcf

Vcf2maf can also take input from STDIN and output to STDOUT, enabling command line piping and redirection. One way this can be useful is to take input from compressed vcf files:

    gunzip --stdout TCGA-AA-A024_W_IlluminaGA-DNASeq_exome.vcf.gz | python vcf2maf.py -v > TCGA-AA-A024_W_IlluminaGA-DNASeq_exome.maf

Further help can be displayed using the -h or --help option:

    python vcf2maf.py -h

Dependencies
------------
  * [PyVCF](https://github.com/jamescasbon/PyVCF): A Variant Call Format reader for Python.

License and Copyright
---------------------

&copy; Copyright 2012 Sage Bionetworks

This software is licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).
