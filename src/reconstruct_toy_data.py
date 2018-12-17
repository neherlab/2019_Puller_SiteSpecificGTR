import os, gzip, glob
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


def reconstruct_counts(in_prefix, params, gtr='JC69', alphabet='nuc_nogap', marginal=False):
    with gzip.open(alignment_name(in_prefix, params), 'rt') as fh:
        aln = AlignIO.read(fh, 'fasta')
    myTree = TreeAnc(gtr=gtr, alphabet=alphabet,
                     tree=tree_name(in_prefix, params), aln=aln,
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
    prefix = 'simulated_data/'
    files = glob.glob(prefix+'*fasta.gz')

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


    from matplotlib import pyplot as plt
    L=100
    n_vals = [100, 300, 1000, 3000]
    for n in n_vals:
        plt.figure()
        mu_vals = sorted(mu_vals)
        for label, data in p_dist.items():
            plt.errorbar(mu_vals, [np.mean(data[(L,n,mu)]) for mu in mu_vals],
                        [np.std(data[(L,n,mu)]) for mu in mu_vals], label=label)

        plt.yscale('log')
        plt.legend()
        plt.xlabel('average rate')


        plt.figure()
        mu_vals = sorted(mu_vals)
        for label, data in mu_dist.items():
            plt.errorbar(mu_vals, [np.mean(data[(L,n,mu)]) for mu in mu_vals],
                        [np.std(data[(L,n,mu)]) for mu in mu_vals], label=label)
        for label, data in W_dist.items():
            plt.errorbar(mu_vals, [np.mean(data[(L,n,mu)]) for mu in mu_vals],
                        [np.std(data[(L,n,mu)]) for mu in mu_vals], label=label)

        plt.yscale('log')
        plt.legend()
        plt.xlabel('average rate')


    plt.figure()
    for li, (label, data) in enumerate(p_dist.items()):
        for n in n_vals:
            d = []
            for mu in mu_vals:
                d.append((n*mu, np.mean(data[(L,n,mu)]), np.std(data[(L,n,mu)])))

            d = np.array(sorted(d, key=lambda x:x[0]))
            plt.errorbar(d[:,0], d[:,1], d[:,2], label=label if n==n_vals[0] else '', c='C%d'%(li+1))

    plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.xlabel('tree length')
