#!/usr/bin/env python

import sys


def readPrimerScores(primer_scores_file):
    primer_scores = {}
    ip = open(primer_scores_file, 'r')
    for line in ip:
        if (line[0] != '#'):
            primer_seq, primer_score = line.strip().split("\t")
            primer_scores[primer_seq] = primer_score
    ip.close()

    return primer_scores


def readPrimersToScore(primers_file_tsv):
    header = None
    primers_to_score = []

    ip = open(primers_file_tsv, 'r')
    for i, line in enumerate(ip):
        fields = line.strip().split("\t")
        if (i == 0):
            header = fields
            assert (header[-3] == "FP")
            assert (header[-2] == "RP")
        else:
            primers_to_score.append(fields)
    ip.close()
    
    return header, primers_to_score


def scorePrimers(primers_to_score, header, DNAPol_primer_scores, RT_primer_scores):
    primers_w_scores = []
    header.extend( ["FP_Efficiency_Score", "RP_Efficiency_Score"] )
    for tup in primers_to_score:
        DNAPol_primer_seq = tup[-3][-6:]
        DNAPol_primer_score = DNAPol_primer_scores[ DNAPol_primer_seq ]

        RT_primer_seq = tup[-2][-6:]
        RT_primer_score = RT_primer_scores[ RT_primer_seq ]

        new_tup = tup + [DNAPol_primer_score, RT_primer_score]
        primers_w_scores.append( new_tup )

    return header, primers_w_scores


def writeScoredPrimers(primers_w_scores, header, output_tsv):
    op = open(output_tsv, 'w')
    op.write("%s\n" % "\t".join(header))
    for tup in primers_w_scores:
        op.write("%s\n" % "\t".join(tup))
    op.close()
    

if (__name__ == "__main__"):
    RT_primer_scores_file, DNAPol_primer_scores_file, primers_file_tsv, output_tsv = sys.argv[1:]

    RT_primer_scores = readPrimerScores(RT_primer_scores_file)
    DNAPol_primer_scores = readPrimerScores(DNAPol_primer_scores_file)

    header, primers_to_score = readPrimersToScore(primers_file_tsv)
    header, primers_w_scores = scorePrimers(primers_to_score, header, DNAPol_primer_scores, RT_primer_scores)

    writeScoredPrimers(primers_w_scores, header, output_tsv)

    sys.exit(0)
