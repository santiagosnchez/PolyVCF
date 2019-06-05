# PolyVCF

This program estimates several metrics of nucleotide polymorphism from a VCF file on a sliding window basis.

## Requirements

It strictly depends on the Python `pysam` [library](https://pysam.readthedocs.io/en/latest/api.html), which can be easily installed with pip

    pip install --upgrade pip
    pip install pysam
    
or conda.

    conda update conda
    conda update conda-build
    conda install pysam

## VCF preprocessing

Because `pysam` uses tabix indexed files to efficiently access intervals within the VCF file it is important to first use `bgzip` and `tabix`.

In a unix shell, simply do:

   bgzip myVariants.vcf
   tabix myVariants.vcf.gz

This will create an index (.tbi) file of your VCF.

## Getting the code

    wget https://raw.githubusercontent.com/santiagosnchez/PolyVCF/master/PolyVCF.py

or

    git clone https://github.com/santiagosnchez/PolyVCF.git

## Running the code

Once requirements are installed and preprocessing is ready, run it first with the `-h` help argument so it will show some of the options.

    python PolyVCF.py -h
    usage: PolyVCF.py [-h] [--vcf VCF] [--group [pop1,pop2,pop3]]
                      [--window WINDOW] [--ploidy PLOIDY]
    
        Fast estimator of nucleotide diversity (pi), theta, and Tajimas's D 
        based on the site frequency spectrum.
    
    optional arguments:
      -h, --help            show this help message and exit
      --vcf VCF, -v VCF     tabix indexed and bgzip compressed VCF file.
      --group [pop1,pop2,pop3], -g [pop1,pop2,pop3]
                            split alignment by populations. A comma-separated list of strings that are found in the sequence headers.
      --window WINDOW, -w WINDOW
                            window size [integer]. (default: 1000 bp)
      --ploidy PLOIDY, -p PLOIDY
                            ploidy [integer]. This is used as a factor to get the total number of chromosomes (default: 2)
    
        IMPORTANT:
    
        It is necesary to first generate a TABIX index using using a UNIX shell:
    
        $ bgzip myFile.vcf
        $ tabix myFile.vcf.gz
    
        Examples:
        python PolyVCF.py -v myFile.vcf.gz -g pop1,pop2 -w 1000 > diversity_w1000.txt
        "-p pop1,pop2" assumes that the VCF has samples that are labeled in VCF header:
    
        #CHROM	(...)	pop1_ind1_XXX	pop1_ind2_XXX	pop2_ind1_XXX	pop2_ind2_XXX
        
        The format is not strict but the identifier (e.g. pop1) needs to be somewhere in the header name.
    
        The output should be tab-delimited.
