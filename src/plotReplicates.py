#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import sys


def readScores(scores_file):
    scores = {}

    with open(scores_file, 'r') as ip:
        for line in ip:
            if (line[0] != '#'):
                primer_seq, FC, max_FC, score = line.strip().split("\t")
                if ('N' not in primer_seq):
                    scores[primer_seq] = float(score)
    return scores


def readFCs(scores_file):
    scores = {}

    with open(scores_file, 'r') as ip:
        for line in ip:
            if (line[0] != '#'):
                primer_seq, FC, max_FC, score = line.strip().split("\t")
                if ('N' not in primer_seq):
                    scores[primer_seq] = float(FC)
    return scores


def matchScores(scores1, scores2):
    scores1_list, scores2_list = [], []

    assert (len(scores1.keys()) == len(scores2.keys()))

    for primer_seq in scores1.keys():
        score1 = scores1[primer_seq]
        score2 = scores2[primer_seq]
        scores1_list.append(score1)
        scores2_list.append(score2)

    return (scores1_list, scores2_list)


# http://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
def plotScores(mode, scores1_list, scores2_list, output_png):
    x, y = np.array(scores1_list,'f'), np.array(scores2_list,'f')
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)

    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    fig, ax = plt.subplots()
    ax.scatter(x, y, c=z, s=50, edgecolor='')
    if (mode == "fold_change"):
        ax.set_title("Primer sequence bias relative to background, Fold Change")
        ax.set_xlabel("Fold Change, Replicate 1")
        ax.set_ylabel("Fold Change, Replicate 2")
    plt.savefig(output_png)
    

if (__name__ == "__main__"):
    mode, scores1_file, scores2_file, output_png = sys.argv[1:]
    assert (mode in ["scores", "fold_change"])

    if (mode == "score"):
        scores1 = readScores(scores1_file)
        scores2 = readScores(scores2_file)
    else:
        scores1 = readFCs(scores1_file)
        scores2 = readFCs(scores2_file)

    scores1_list, scores2_list = matchScores(scores1, scores2)
    plotScores(mode, scores1_list, scores2_list, output_png)
    
    sys.exit(0)
    
