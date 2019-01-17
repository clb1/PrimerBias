#!/usr/bin/env python

from operator import itemgetter
import sys


def readScores(scores_file):
    scores = {}
    header_lines = []

    with open(scores_file, 'r') as ip:
        for line in ip:
            if (line[0] == '#'):
                header_lines.append(line)
            else:
                primer_seq, score = line.strip().split("\t")
                if ('N' not in primer_seq):
                    scores[primer_seq] = float(score)

    return header_lines, scores


def combineScores(scores1, scores2):
    combined_scores = {}

    assert (len(scores1.keys()) == len(scores2.keys()))

    for primer_seq in scores1.keys():
        score1 = scores1[primer_seq]
        score2 = scores2[primer_seq]
        combined_scores[primer_seq] = (score1 + score2)/2.0

    return combined_scores

    

def writeCombinedScores(header_lines1, header_lines2, combined_scores, output_file):
    assert (len(header_lines1) == len(header_lines2))

    with open(output_file, 'w') as op:
        for i in xrange(len(header_lines1)):
            assert (header_lines1[i] == header_lines2[i])
            op.write(header_lines1[i])

        primers_and_scores = sorted(combined_scores.items(), key=itemgetter(1), reverse=True)

        for primer_seq, score in primers_and_scores:
            op.write("%s\t%f\n" % (primer_seq, score))
            
            
if (__name__ == "__main__"):
    scores1_file, scores2_file, output_file = sys.argv[1:]

    header_lines1, scores1 = readScores(scores1_file)
    header_lines2, scores2 = readScores(scores2_file)
    
    combined_scores = combineScores(scores1, scores2)
    writeCombinedScores(header_lines1, header_lines2, combined_scores, output_file)
    
    sys.exit(0)
    
