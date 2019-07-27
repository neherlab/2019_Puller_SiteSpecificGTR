import os, gzip, glob, pickle
from Bio import AlignIO
import numpy as np
from scipy.stats import pearsonr, spearmanr
from collections import defaultdict
import matplotlib.pyplot as plt

from treetime.treeanc import TreeAnc
from treetime.gtr_site_specific import GTR_site_specific
from treetime.gtr import GTR
from treetime.seq_utils import seq2prof, profile_maps, alphabets

from estimation import *
from generate_toy_data import save_model, load_model
from analyze_HIV_tree import *



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    hiv_out = 'HIV_inferences/'
    parser.add_argument("--prefix", type=str, help="data set")
    parser.add_argument("--gene", type=str, default='pol')
    parser.add_argument("--aa", action='store_true', help="use amino acid alphabet")
    parser.add_argument("--pc", type=float, default=1.0, help="pseudocount for reconstruction")
    args=parser.parse_args()

    alphabet='aa' if args.aa else 'nuc_nogap'
    pc = args.pc

    aln = args.prefix +  '_aligned.fasta'
    tree = args.prefix + '_tree.nwk'
    model_name = hiv_out + os.path.basename(aln)[:-5]+'%1.2f_inferred_model.npz'%args.pc
    aln_freq_name = hiv_out + os.path.basename(aln[:-5])+'alignment_frequencies'

    aln_freq = {}
    gtr = {}
    fitness = {}
    for st in ['B', 'C']:
        aln_freq[st] = np.loadtxt(aln_freq_name.replace('#st#', st))
        gtr[st] = load_model(model_name.replace('#st#', st))
        fitness[st] = calculate_in_out_ratio(gtr[st], fabio_fitness.loc[:,"consensus"])

    overlap_gtr = np.sum(gtr['B'].Pi*gtr['C'].Pi, axis=0).mean()
    overlap_afreq = np.sum(aln_freq['B']*aln_freq['C'], axis=0).mean()

