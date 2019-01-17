#!/usr/bin/env python

from collections import defaultdict, Counter
from itertools import product
import gzip
import numpy as np
import sys

import pdb

def initializeCounter(sequence_length):
    L = [['A','C','G','T']] * sequence_length
    D = dict(map(lambda y: (''.join(y), []), [x for x in product(*L)]))
    return D


def accumulateCounts(degenerate_region_sequences, primer_counts, primer_len):
    ip = gzip.open(degenerate_region_sequences, 'rb')
    for line in ip:
        line = line.strip()
        if (line != ""):
            readID, UMI, primer_region = line.strip().split("\t")
            primer = primer_region[-primer_len:]
            signature = (UMI, primer_region[0:-primer_len], primer_region[-1])
            primer_counts[primer].append(signature)
    ip.close()


def makeHistograms(primer_counts, unique_occurrences_hist, duplication_rate_hist):

    unique_primer_occurrence_counts = []
    duplicates_occurrence_counts = Counter()
    for primer_seq, signatures in primer_counts.items():
        signature_counts = Counter(signatures)
        unique_primer_occurrence_counts.append(len(signature_counts.keys()))
        duplicates_occurrence_counts.update( signature_counts.values() )

    bin_values, bin_edges = np.histogram(a=unique_primer_occurrence_counts, bins=100, density=False)
    op = open(unique_occurrences_hist, 'w')
    op.write("Number of times a primer seen\tNumber of primers seen that many times\n")
    for num_times_seen, number_of_primers in zip(bin_edges, bin_values[0:-1]):
        op.write("%d :  %d\n" % (num_times_seen, number_of_primers))
    op.close()

    op = open(duplication_rate_hist, 'w')
    op.write("Number of PCR duplicates\tNumber of times seen\n")
    for number_of_duplicates, num_times_seen in duplicates_occurrence_counts.items():
        op.write("%d :  %d\n" % (number_of_duplicates, num_times_seen))
    op.close()
    

if (__name__ == "__main__"):
    degenerate_region_sequences, primer_len, unique_occurrences_hist, duplication_rate_hist = sys.argv[1:]

    primer_len = int(primer_len)
    primer_counts = initializeCounter(primer_len)
    
    accumulateCounts(degenerate_region_sequences, primer_counts, primer_len)

    makeHistograms(primer_counts, unique_occurrences_hist, duplication_rate_hist)
    
    sys.exit(0)
