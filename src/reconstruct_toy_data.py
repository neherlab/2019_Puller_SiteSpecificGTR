import os, gzip, glob, pickle
import numpy as np
from Bio import Phylo, AlignIO
from collections import defaultdict

from treetime.treeanc import TreeAnc
from treetime.gtr_site_specific import GTR_site_specific
from treetime.gtr import GTR
from treetime.seq_utils import seq2prof, profile_maps, alphabets

from generate_toy_data import *
from filenames import *

def p_from_aln(in_prefix, params, alphabet='nuc_nogap'):
    with gzip.open(alignment_name(in_prefix, params), 'rt') as fh:
        aln = AlignIO.read(fh, 'fasta')

    alpha = alphabets[alphabet]
    aln_array = np.array(aln)
    af = []
    for a in alpha:
        af.append(np.mean(aln_array==a, axis=0))
    return np.array(af)


def reconstruct_counts(in_prefix, params, gtr='JC69', alphabet='nuc_nogap', marginal=False, reconstructed_tree=False):
    with gzip.open(alignment_name(in_prefix, params), 'rt') as fh:
        aln = AlignIO.read(fh, 'fasta')
    tree_fname =  reconstructed_tree_name(in_prefix, params) if reconstructed_tree else tree_name(in_prefix, params)
    myTree = TreeAnc(gtr=gtr, alphabet=alphabet,
                     tree=tree_fname, aln=aln,
                     reduce_alignment=False, verbose = 0)

    if type(gtr)==str:
        myTree.infer_ancestral_sequences(marginal=True, infer_gtr=True, normalized_rate=False)
    else:
        myTree.infer_ancestral_sequences(marginal=True)

    mutation_counts = get_mutation_count(myTree, alphabet=myTree.gtr.alphabet, marginal=marginal)
    return mutation_counts, myTree.sequence_LH(), myTree


def estimate_GTR(mutation_counts, pc=0.1, single_site=False):
    n_ija, T_ia, root_sequence = mutation_counts
    root_prof = seq2prof(root_sequence, profile_maps['nuc_nogap']).T

    if single_site:
        inferred_model = GTR.infer(n_ija.sum(axis=-1), T_ia.sum(axis=-1),
                          root_state=root_prof.sum(axis=-1), alphabet='nuc_nogap')
    else:
        inferred_model = GTR_site_specific.infer(n_ija, T_ia, pc=pc,
                          root_state=root_prof, alphabet='nuc_nogap')
    return inferred_model



def KL(p,q):
    return np.sum(p*log(p/q))

def chisq(p,q):
    return np.sum((p-q)**2)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    parser.add_argument("-m", type=float, help="simulated mutation rate")
    parser.add_argument("-n", type=int, help="number of taxa")
    parser.add_argument("-L", type=int, help="length of sequence")
    args=parser.parse_args()

    prefix = '2018-12-17_simulated_data'
    mask = "/L{L}_n{n}_m{mu}_*fasta.gz".format(L=args.L or '*', n=args.n or '*', mu=args.m or "*")
    files = glob.glob(prefix+mask)

    mu_dist = defaultdict(lambda: defaultdict(list))
    p_dist = defaultdict(lambda: defaultdict(list))
    W_dist = defaultdict(lambda: defaultdict(list))

    mu_vals = set()
    n_vals = set()

    pc=0.1

    analysis_types = ['naive', 'single', 'dressed', 'regular', 'phylo', 'marginal', 'iterative']

    for fname in files:
        print(fname)

        params = parse_alignment_name(fname)
        mu_vals.add(params['m'])
        n_vals.add(params['n'])

        true_model = load_model(model_name(prefix, params))
        true_mut_counts = load_mutation_count(mutation_count_name(prefix, params))

        dset = (params['L'], params['n'], params['m'])
        for ana in analysis_types:
            if ana=='naive':
                naive = p_from_aln(prefix, params)
                p_dist[ana][dset].append(np.mean([chisq(naive[:,i],true_model.Pi[:,i]) for i in range(params['L'])]))
            else:
                if ana=='dressed':
                    mc = (true_mut_counts, None, None)
                else:
                    mc = reconstruct_counts(prefix, params, gtr=model if ana=='iterative' else 'JC69',
                                        alphabet='nuc_nogap', marginal=ana=='marginal')
                model = estimate_GTR(mc[0], pc=pc, single_site=ana=='single')

                if ana!='single':
                    p_dist[ana][dset].append(np.mean([chisq(model.Pi[:,i],true_model.Pi[:,i]) for i in range(params['L'])]))
                    mu_dist[ana][dset].append(chisq(model.mu/np.mean(model.mu), true_model.mu/np.mean(true_model.mu)))

                np.fill_diagonal(model.W,0)
                W_dist[ana][dset].append(chisq(model.W.flatten(), true_model.W.flatten()))

    out_prefix = prefix + '_results/'
    if not os.path.isdir(out_prefix):
        os.mkdir(out_prefix)

    out_fname = out_prefix + "_".join(["{name}{val}".format(name=n, val=args.__getattribute__(n)) for n in ['L', 'n', 'm'] if args.__getattribute__(n)]) + '.pkl'
    with open(out_fname, 'wb') as fh:
        pickle.dump((sorted(mu_vals), sorted(n_vals), dict(p_dist), dict(mu_dist), dict(W_dist)), fh)
