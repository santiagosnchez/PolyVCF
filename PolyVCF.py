# PolyVCF
# Santiago Sanchez-Ramirez
# email santiago.snchez@gmail.con

import sys
import math
import copy
import argparse
try:
    import pysam
except ImportError:
    print "Try installing pysam first: pip install pysam"
    sys.exit()

def main():
    # parse arguments
    parser = argparse.ArgumentParser(prog="PolyVCF.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="""
    Fast estimator of nucleotide diversity (pi), theta, and Tajimas\'s D 
    based on the site frequency spectrum.\n""",
    epilog="""

    IMPORTANT:

    It is necesary to first generate a TABIX index using using a UNIX shell:

    $ bgzip myFile.vcf
    $ tabix myFile.vcf.gz

    Examples:
    python PolyVCF.py -v myFile.vcf.gz -g pop1,pop2 -w 1000 > diversity_w1000.txt
    \"-p pop1,pop2\" assumes that the VCF has samples that are labeled in VCF header:

    #CHROM\t(...)\tpop1_ind1_XXX\tpop1_ind2_XXX\tpop2_ind1_XXX\tpop2_ind2_XXX
    
    The format is not strict but the identifier (e.g. pop1) needs to be somewhere in the header name.

    The output should be tab-delimited.\n""")
    parser.add_argument(
    '--vcf', '-v', type=str, default="",
    help='tabix indexed and bgzip compressed VCF file.')
    parser.add_argument(
    '--group', '-g', nargs="?", default=False, metavar='pop1,pop2,pop3', type=str,
    help='split alignment by populations. A comma-separated list of strings that are found in the sequence headers.')
    parser.add_argument(
    '--window', '-w', type=int, default=1000,
    help='window size [integer]. (default: 1000 bp)')
    parser.add_argument(
    '--ploidy', '-p', type=int, default=2,
    help='ploidy [integer]. This is used as a factor to get the total number of chromosomes (default: 2)')

    args = parser.parse_args()

    # read vcf and extract info
    vcf = pysam.VariantFile(args.vcf)
    contigs = list((vcf.header.contigs))
    vfirst = vcf.copy().next()
    samples = vfirst.samples.keys()
    pops = args.group.split(",")
    if args.group:
        group = [ filter(lambda x: k in x, samples) for k in pops ]
    else:
        group = [ samples ]
        pops = [ "all" ]
    if any([ len(x) == 0 for x in group ]):
        print "samples: ",samples
        print "pops: ",pops
        parser.error("There is one or more groups that were not found in the samples.")
    w = args.window
    ploidy = args.ploidy

    # generate an empty SFS vector for each group with the number of chromosomes 

    sfs_base = [ ( len(k)*ploidy, [ 0 for i in range((len(k)*ploidy)/2) ]) for k in group ]
    
    # loop through all the contigs

    print "chr\tstart\tend\tpop\tsegsites\tpi\ttheta\ttajimasD" # header
    for c in contigs:
        stop = 0 # initialize stopping rule
        try:
            vfirst = vcf.fetch(c).next() # get first element of VCF for each contig
        except StopIteration:
            pass
        else:
            fpos = rounddown(vfirst.pos, w) + 1  # first position
            lpos = fpos + (w-1) # last position
            vlast = getLast(vcf.fetch(c)) # get last element of VCF for each contig
            while stop != 1:
                # first check if there is data in the interval
                v_check = vcf.fetch(c, fpos, lpos)
                try:
                    v_check.next()
                except StopIteration:
                    # if not slide the window
                    fpos = lpos + 1
                    lpos = fpos + w
                else:
                    # loop through intervals that have data
                    v_good = vcf.fetch(c, fpos, lpos)
                    sfs = copy.deepcopy(sfs_base) # deep copy of the sfs list
                    sfs = getCounts(sfs, v_good, group) # count frequency of variants
                    printVar(sfs, group, c, fpos, lpos, pops, w) # print result
                    # check if there is data in the next interval
                    v_check = vcf.fetch(c, lpos+1, lpos+w)
                    try:
                        v_check.next()
                    except StopIteration:
                        if lpos+1 >= vlast.pos: # stop if the next interval is greater than the last position
                            stop = 1
                        else:
                            # slide window
                            fpos = lpos + 1
                            lpos = lpos + w
                    else:
                        # slide window
                        fpos = lpos + 1
                        lpos = lpos + w

# functions

# prints to screen for each group
def printVar(sfs, group, c, fpos, lpos, pops, w):
    for grp in range(len(group)):
        poly = polymorphism(sfs[grp],w) 
        print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(c, fpos, lpos, pops[grp], poly[0], poly[1], poly[2], poly[3])

def getCounts(sfs, v, group):
    for line in v: # iterate over subset of variants
        for grp in range(len(group)): # for each group
            line_var = []
            for samp in group[grp]:
                # put genotype into a list
                line_var += line.samples.get(samp).get('GT')
            line_var = filter(lambda x: x != '.', line_var) # remove missing data
            if len(set(line_var)) == 2: # only take biallelic sites
                count = sum([ x != 0 for x in line_var ])
                # construct the folded SFS
                if count > len(sfs[grp][1]):
                    sfs[grp][1][ (count - len(sfs[grp][1])) - 1 ] += 1
                else:
                    sfs[grp][1][ count - 1 ] += 1
    return sfs

def getLast(x): # get last item in iterator
    i = None
    for i in x:
        pass
    return i

def rounddown(x, w): # get the closest integer to the first variant                               
    return int(math.floor(x / float(w))) * w

def polymorphism(sfs,seqlen): # estimate polymorphism metrics
    N,sfs = sfs
    if sum(sfs) == 0:
        return 0,0,0,"NA"
    else:
        ss = sum(sfs)
        a1,dv = Dvar(N,ss)
        pi = (2.0/(N*(N-1.0)))*sum( map(lambda i: sfs[i]*(i+1)*(N-(i+1)), range(len(sfs))) )
        th = float(ss)/a1
        try:
            D = (pi-th)/dv
        except ZeroDivisionError:
            D = "NA"
        return ss,pi/seqlen,th/seqlen,D

def Dvar(N,ss): # Needed for Tajima's D and Theta
    a1 = sum(map(lambda x: 1.0/x, range(1,N)))
    a2 = sum(map(lambda x: 1.0/(x**2), range(1,N)))
    b1 = (N+1.0)/(3.0*(N-1.0))
    b2 = (2.0*((N**2.0)+N+3.0))/(9.0*N*(N-1.0))
    c1 = b1 - (1/a1)
    c2 = b2 - ((N+2)/(a1*N)) + (a2/(a1**2))
    e1 = c1/a1
    e2 = c2/((a1**2)+a2)
    Dv = math.sqrt((e1*ss)+(e2*ss*(ss-1)))
    return a1,Dv

if __name__ == '__main__':
    main()
