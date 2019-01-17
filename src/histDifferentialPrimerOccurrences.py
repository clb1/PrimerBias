#!/usr/bin/env python

import matplotlib.pyplot as plt
from collections import defaultdict
from itertools import product
import gzip
import numpy as np
import sys


def initializeCounter(sequence_length):
    L = [['A','C','G','T']] * sequence_length
    D = dict(map(lambda y: (''.join(y), set()), [x for x in product(*L)]))
    return D


def accumulateCounts(degenerate_region_sequences, primer_counts, primer_len):
    li = primer_len + 1
    ip = gzip.open(degenerate_region_sequences, 'rb')
    for line in ip:
        line = line.strip()
        if (line != ""):
            readID, UMI, primer_region = line.strip().split("\t")
            primer = primer_region[-li:-1]
            signature = (UMI, primer_region[0:-li], primer_region[-1])
            try:
                primer_counts[primer].add(signature)
            except KeyError, ke:
                if (not 'N' in primer):
                    raise ke
    ip.close()


def computeFoldChangeOverBackground(background_primer_counts, bias_primer_counts):
    sum_background_counts = float(sum(map(len, background_primer_counts.values())))
    sum_bias_counts = float(sum(map(len, bias_primer_counts.values())))

    fold_change_over_background = {}
    min_background_frac = 1.0
    leftover = []
    num_no_counts_background, num_no_counts_bias = 0, 0
    for primer_seq in bias_primer_counts.keys():
        if (len(background_primer_counts[primer_seq]) == 0):
            num_no_counts_background += 1
        if (len(bias_primer_counts[primer_seq]) == 0):
            num_no_counts_bias += 1

        background_frac = float(len(background_primer_counts[primer_seq]))/sum_background_counts
        if (background_frac != 0):
            min_background_frac = min(background_frac, min_background_frac)
            bias_frac = float(len(bias_primer_counts[primer_seq]))/sum_bias_counts
            if (bias_frac >= background_frac):
                FC = bias_frac/background_frac
                fold_change_over_background[primer_seq] = FC
            elif (bias_frac != 0):
                FC = -background_frac/bias_frac
                fold_change_over_background[primer_seq] = FC
            else:
                print >> sys.stderr, "WARNING: no bias counts for primer sequence %s" % primer_seq

        elif (bias_frac != 0):
            leftover.append( (primer_seq, bias_frac) )

    for primer_seq, bias_frac in leftover:
        FC = bias_frac/background_frac if (bias_frac >= background_frac) else -background_frac/bias_frac
        fold_change_over_background[primer_seq] = FC

    print >> sys.stderr, "INFO: Number of primers with no background counts = %d" % num_no_counts_background
    print >> sys.stderr, "INFO: Number of primers with no bias counts = %d" % num_no_counts_bias
    
    return fold_change_over_background


def makeHist(fold_change_over_background, primer_len, output_filename):
    L1 = list(.05 * np.array(range(20,301),'f'))
    L2 = L1[:]
    L2.reverse()
    L = [-100] + list(-1.0 * np.array(L2[0:-1])) + list(L1) + [100]
    bins = np.array(L, 'f')
    
    bin_values, bin_edges = np.histogram(a=fold_change_over_background.values(), bins=bins, density=False)

    num_primers_accounted = sum(bin_values)
    if (num_primers_accounted != 4**primer_len):
        print >> sys.stderr, "WARNING: histogram only accounts for %d of of the %d possible primers." % (num_primers_accounted, 4**primer_len)
    
    op = open(output_filename, 'w')
    op.write("Bias(Fold change over background)\tNumber of primers seen that many times\n")
    for fold_change, number_of_primers in zip(bin_edges, bin_values[0:-1]):
        if (fold_change < 1.0 and fold_change > 0):
            fold_change = -1.0/fold_change
        op.write("%4.2f :  %d\n" % (round(fold_change,2), number_of_primers))
    op.close()
    
    # Make figure file
    figure_filename = output_filename.rsplit('.',1)[0] + ".png"
    fig = plt.figure()
    ax = fig.add_subplot(111)
    L = list(-1.0 * np.array(L2[0:-1])) + list(L1)
    bins = np.array(L, 'f')
    n, bins, patches = ax.hist(x=fold_change_over_background.values(), bins=bins, normed=0, range=(-15,15), facecolor='white', edgecolor='black', log=True)
    plt.xlabel("Primer sequence abundance relative to background occurrence (Fold Change)")
    plt.ylabel("Number of primers")
    plt.savefig(figure_filename)


if (__name__ == "__main__"):
    background_degenerate_region_sequences, bias_degenerate_region_sequences, primer_len, output_filename = sys.argv[1:]

    primer_len = int(primer_len)
    background_primer_counts = initializeCounter(primer_len)
    bias_primer_counts = initializeCounter(primer_len)

    accumulateCounts(background_degenerate_region_sequences, background_primer_counts, primer_len)
    accumulateCounts(bias_degenerate_region_sequences, bias_primer_counts, primer_len)

    fold_change_over_background = computeFoldChangeOverBackground(background_primer_counts, bias_primer_counts)
    makeHist(fold_change_over_background, primer_len, output_filename)
    
    sys.exit(0)
