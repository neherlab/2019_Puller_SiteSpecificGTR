import os, gzip, glob
import numpy as np
from Bio import Phylo, AlignIO
from collections import defaultdict

from treetime.treeanc import TreeAnc
from treetime.gtr_site_specific import GTR_site_specific
from treetime.gtr import GTR
from treetime.seq_utils import seq2prof, profile_maps, alphabets

from generate_toy_data import *

def p_from_aln(in_prefix, params, alphabet='nuc_nogap'):
    with gzip.open(alignment_name(in_prefix, params), 'rt') as fh:
        aln = AlignIO.read(fh, 'fasta')

    alpha = alphabets[alphabet]
    aln_array = np.array(aln)
    af = []
    for a in alpha:
        af.append(np.mean(aln_array==a, axis=0))
    return np.array(af)


def reconstruct_counts(in_prefix, params, gtr='JC69', alphabet='nuc_nogap'):
    with gzip.open(alignment_name(in_prefix, params), 'rt') as fh:
        aln = AlignIO.read(fh, 'fasta')
    myTree = TreeAnc(gtr=gtr, alphabet=alphabet,
                     tree=tree_name(in_prefix, params), aln=aln,
                     reduce_alignment=False, verbose = 2)

    myTree.infer_ancestral_sequences(marginal=True)

    return get_mutation_count(myTree.tree, alphabet=myTree.gtr.alphabet), myTree.sequence_LH()


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

    p_dressed_dist = defaultdict(list)
    p_dist = defaultdict(list)
    p_iter_dist = defaultdict(list)
    p_aln_dist = defaultdict(list)

    mu_corr = defaultdict(list)
    mu_dist = defaultdict(list)
    mu_dressed_dist = defaultdict(list)
    LH = defaultdict(list)

    W_corr = defaultdict(list)
    W_dist = defaultdict(list)
    W_ss_dist = defaultdict(list)
    mu_vals = set()

    pc=0.1

    for fname in files:
        params = parse_alignment_name(fname)
        mu_vals.add(params['m'])

        true_model = load_model(model_name(prefix, params))
        true_mut_counts = load_mutation_count(mutation_count_name(prefix, params))
        dressed = estimate_GTR(true_mut_counts, pc=pc)

        inferred_mut_counts, _ = reconstruct_counts(prefix, params, gtr='JC69', alphabet='nuc_nogap')
        inferred = estimate_GTR(inferred_mut_counts, pc=pc)
        inferred_ss = estimate_GTR(inferred_mut_counts, pc=pc, single_site=True)

        naive = p_from_aln(prefix, params)

        dset = (params['L'], params['n'], params['m'])

        LH[dset].append((reconstruct_counts(prefix, params, gtr=inferred)[1], reconstruct_counts(prefix, params, gtr=inferred_ss)[1]))
        p_dist[dset].append(np.mean([chisq(inferred.Pi[:,i],true_model.Pi[:,i]) for i in range(params['L'])]))
        p_dressed_dist[dset].append(np.mean([chisq(dressed.Pi[:,i],true_model.Pi[:,i]) for i in range(params['L'])]))
        p_aln_dist[dset].append(np.mean([chisq(naive[:,i],true_model.Pi[:,i]) for i in range(params['L'])]))

        mu_corr[dset].append(np.corrcoef(inferred.mu, true_model.mu)[0,1])
        mu_dist[dset].append(chisq(inferred.mu/np.mean(inferred.mu), true_model.mu/np.mean(true_model.mu)))
        mu_dressed_dist[dset].append(chisq(dressed.mu/np.mean(dressed.mu), true_model.mu/np.mean(true_model.mu)))

        np.fill_diagonal(inferred_ss.W,0)
        W_corr[dset].append(np.corrcoef(inferred.W[inferred.W>0], true_model.W[inferred.W>0])[0,1])
        W_dist[dset].append(chisq(inferred.W.flatten(), true_model.W.flatten()))
        W_ss_dist[dset].append(chisq(inferred_ss.W.flatten(), true_model.W.flatten()))


    from matplotlib import pyplot as plt
    L=100
    n_vals = [100, 300, 1000, 3000]
    for n in n_vals:
        plt.figure()
        mu_vals = sorted(mu_vals)
        for label, data in [('reconstructed', p_dist), ('dressed', p_dressed_dist),
                            ('alignment', p_aln_dist)]:
            plt.errorbar(mu_vals, [np.mean(data[(L,n,mu)]) for mu in mu_vals],
                        [np.std(data[(L,n,mu)]) for mu in mu_vals], label=label)

        plt.yscale('log')
        plt.legend()
        plt.xlabel('average rate')


        plt.figure()
        mu_vals = sorted(mu_vals)
        for label, data in [('mu reconstructed', mu_dist), ('mu dressed', mu_dressed_dist)]:
            plt.errorbar(mu_vals, [np.mean(data[(L,n,mu)]) for mu in mu_vals],
                        [np.std(data[(L,n,mu)]) for mu in mu_vals], label=label)
        for label, data in [('W site specific', W_dist),('W single site specific', W_ss_dist)]:
            plt.errorbar(mu_vals, [np.mean(data[(L,n,mu)]) for mu in mu_vals],
                        [np.std(data[(L,n,mu)]) for mu in mu_vals], label=label)

        plt.yscale('log')
        plt.legend()
        plt.xlabel('average rate')


    plt.figure()
    for li, (label, data) in enumerate([('reconstructed', p_dist), ('dressed', p_dressed_dist),
                        ('alignment', p_aln_dist)]):
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
