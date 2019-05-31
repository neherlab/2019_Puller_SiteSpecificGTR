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

def convert_vals(x):
    try:
        return float(x)
    except:
        if '<' in x:
            return 0.0005
        else:
            return 0.2

def load_fitness_landscape(gene, aa=False, subtype='B'):
    import pandas as pd
    d = pd.read_csv('HIV_fitness_landscape/data/fitness_pooled/nuc_%s_selection_coeffcients_%s.tsv'%(gene, subtype),
                    sep='\t', skiprows=1)
    for col in d.columns[3:6]:
        d.loc[:,col] = d.loc[:,col].apply(convert_vals)

    return d


def calculate_in_out_ratio(gtr, major_allele):
    '''
    e^f_a = sum_j pi_ia W_ij/sum W_ij pi_ja
    where i is the major allele at position i
    '''
    major_index = np.zeros_like(major_allele, dtype=int)
    for ni, n in enumerate(gtr.alphabet):
        major_index[major_allele==n]=ni

    r_ind = np.arange(gtr.Pi.shape[-1])
    return np.einsum('a,aj->a', gtr.Pi[major_index, r_ind], gtr.W[major_index])/\
           np.einsum('aj,ja->a', gtr.W[major_index], gtr.Pi)


def diversity(p):
    return 1.0-np.sum(p**2, axis=0)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    hiv_out = 'HIV_inferences/'
    parser.add_argument("--prefix", type=str, help="data set")
    parser.add_argument("--gene", type=str, default='pol')
    parser.add_argument("--subtype", type=str, default='B')
    parser.add_argument("--aa", action='store_true', help="use amino acid alphabet")
    parser.add_argument("--pc", type=float, default=1.0, help="pseudocount for reconstruction")
    parser.add_argument("--redo", action='store_true', default=False, help="re-estimate")
    args=parser.parse_args()

    alphabet='aa' if args.aa else 'nuc'
    pc = args.pc

    aln = args.prefix +  '_aligned.fasta'
    tree = args.prefix + '_tree.nwk'
    model_name = hiv_out + os.path.basename(aln)[:-5]+'%1.2f_inferred_model.npz'%args.pc
    aln_freq_name = hiv_out + os.path.basename(aln[:-5])+'alignment_frequencies'

    gtr=None
    aln_freq=None
    if not args.redo:
        try:
            aln_freq = np.loadtxt(aln_freq_name)
            gtr = load_model(model_name)
        except:
            pass
    if gtr is None:
        tt = TreeAnc(tree=tree, aln = aln, compress=False, alphabet=alphabet, verbose=3)
        tt.optimize_tree(branch_length_mode='marginal', max_iter=10,
                         infer_gtr=True, site_specific_gtr=True, pc=pc)

        gtr = tt.gtr
        save_model(gtr, model_name)

    if aln_freq is None:
        aln_freq = p_from_aln(AlignIO.read(aln, 'fasta'), alphabet=alphabet)
        np.savetxt(aln_freq_name, aln_freq)

    fabio_fitness = load_fitness_landscape(args.gene, subtype=args.subtype, aa=args.aa)

    fitness_estimates_gtr = calculate_in_out_ratio(gtr, fabio_fitness.loc[:,"consensus"])

    ccoef=spearmanr
    div_rec = diversity(gtr.Pi)
    cc_rec = ccoef(fabio_fitness['median'], div_rec)
    plt.figure()
    plt.scatter(fabio_fitness['median'], div_rec, label='reconstructed, r2=%1.2f'%cc_rec[0]**2)
    plt.legend()
    plt.xlim([3e-4, 3e-1])
    plt.xscale('log')

    div_aln = diversity(aln_freq)
    cc_aln = ccoef(fabio_fitness['median'], div_aln)
    plt.figure()
    plt.scatter(fabio_fitness['median'], div_aln, label='alignment, r2=%1.2f'%cc_aln[0]**2)
    plt.legend()
    plt.xlim([3e-4, 3e-1])
    plt.xscale('log')


    plt.figure()
    plt.scatter(fabio_fitness['median'], fitness_estimates_gtr, label='fitness vs fitness')
    plt.legend()
    plt.xlim([3e-4, 3e-1])
    plt.xscale('log')
    plt.yscale('log')