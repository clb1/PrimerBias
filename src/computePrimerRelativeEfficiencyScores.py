#!/usr/bin/env python

from collections import defaultdict, Counter
from itertools import product
from operator import itemgetter
import gzip
import numpy as np
import sys

from histDifferentialPrimerOccurrences import *

import pdb


def computeRelativeEfficiencies(fold_change_over_background):
    relative_efficiencies = {}

    max_FC = max(fold_change_over_background.values())
    for primer_seq, FC in fold_change_over_background.items():
        if (FC < 0):
            score = (-1.0/FC)/max_FC
        else:
            score = FC/max_FC
        primer_seq_5p_3p = primer_seq[::-1]
        relative_efficiencies[primer_seq_5p_3p] = (FC, max_FC, score)

    return relative_efficiencies


def writeOutput(relative_efficiencies, output_filename):
    op = open( output_filename, 'w')
    op.write("# Sequences are from the priming strand, written in the sense (ie 5'->3') direction\n")
    op.write("# Priming_Sequence\tFold_Change\tMax_Fold_Change\tRelative_Efficiency\n")
    primers_and_scores = relative_efficiencies.items()
    primers_and_scores = sorted(primers_and_scores, key=lambda x:x[1][2], reverse=True)
    for primer_seq, (FC, max_FC, score) in primers_and_scores:
        op.write("%s\t%f\t%f\t%f\n" % (primer_seq, FC, max_FC, score))
    op.close()
    

if (__name__ == "__main__"):
    background_degenerate_region_sequences, bias_degenerate_region_sequences, primer_len, output_filename = sys.argv[1:]

    primer_len = int(primer_len)
    background_primer_counts = initializeCounter(primer_len)
    bias_primer_counts = initializeCounter(primer_len)

    accumulateCounts(background_degenerate_region_sequences, background_primer_counts, primer_len)
    accumulateCounts(bias_degenerate_region_sequences, bias_primer_counts, primer_len)

    fold_change_over_background = computeFoldChangeOverBackground(background_primer_counts, bias_primer_counts)

    relative_efficiencies = computeRelativeEfficiencies(fold_change_over_background)
    
    writeOutput(relative_efficiencies, output_filename)
    
    sys.exit(0)
