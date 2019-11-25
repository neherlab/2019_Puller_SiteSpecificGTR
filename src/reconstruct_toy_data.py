import os, gzip, glob, pickle
import numpy as np
from collections import defaultdict
import pandas as pd

from treetime.treeanc import TreeAnc
from treetime.gtr_site_specific import GTR_site_specific
from treetime.gtr import GTR
from treetime.seq_utils import seq2prof, profile_maps, alphabets

from generate_toy_data import *
from filenames import *
from estimation import *

def KL(p,q):
    """Kullback-Leibler divergence between two discrete distributions p and q

    Parameters
    ----------
    p : array
        discrete probability distribution sum(p)=1
    q : array
        discrete probability distribution sum(q)=1

    Returns
    -------
    float
        KL divergenge
    """
    if np.isclose(p.sum(),1) and np.isclose(q.sum(), 1):
        return np.sum(p*log(p/q))
    else:
        raise ValueError("discrete probability distributions have to sum to 1")


def chisq(p,q):
    """Euclidian distance between two probability distributions

    Parameters
    ----------
    p : array
        discrete distribution
    q : array
        discrete distribution

    Returns
    -------
    float
        squared distance
    """
    return np.sum((p-q)**2)


def assess_reconstruction(model, true_model):
    """Calculate a number of distance metrics between two models
    """
    true_model_average_rate = true_model.average_rate().mean()
    res = { "average_mu": model.average_rate().mean(),
            "true_average_mu": true_model_average_rate,
            "chisq_W":chisq(model.W.flatten()/model.W.sum(), true_model.W.flatten()/true_model.W.sum())}

    if len(model.Pi.shape)==2: # site specific model
        res.update({ "chisq_mu": chisq(model.mu/model.average_rate().mean(), true_model.mu/true_model_average_rate),
            "r_mu": np.corrcoef(model.mu, true_model.mu)[0,1],
            "chisq_p": np.mean([chisq(model.Pi[:,i],true_model.Pi[:,i]) for i in range(params['L'])]),
            "model_entropy": -np.mean(np.sum(model.Pi*np.log(model.Pi), axis=0)),
            "true_entropy": -np.mean(np.sum(true_model.Pi*np.log(true_model.Pi), axis=0))})
    return res


def get_LH(tname,aln, gtr):
    """reconstruct ancestral sequences given a GTR model and report the likelihood
    """
    tt = TreeAnc(tree=tname, aln=aln, gtr=gtr, compress=False)
    tt.infer_ancestral_sequences(marginal=True)
    return tt.sequence_LH()


def analyze(ana, tree, aln, alphabet, prefix, params, true_model):
    """
    infer a model as specified by analysis type "ana" and compare the true model
    """
    T = Phylo.read(tree, format="newick")

    # if the true tree is used, rescale its branches with the mutation rate
    # to convert them to evolutionary distance
    if tree == tree_name(prefix, params):
        for n in T.find_clades():
            n.branch_length *= params['m']

    tt = TreeAnc(tree=T, aln = aln, compress=False, alphabet=alphabet)

    if ana=="iterative":
        # this performs iterative estimation of the model and ancestral sequences
        # using marginal ancestral reconstruction
        tt.infer_gtr_iterative(normalized_rate=False, site_specific=True, pc=pc, max_iter=10)
    elif ana=="unspecific":
        tt.infer_gtr(marginal=False, normalized_rate=False, site_specific=False, pc=pc)
    elif ana=="ml_reconstruction":
        tt.infer_gtr(marginal=False, normalized_rate=False, site_specific=True, pc=pc)
    elif ana=="marginalize":
        tt.infer_gtr(marginal=True, normalized_rate=False, site_specific=True, pc=pc)
    elif ana=='optimize_tree':
        tt.optimize_tree(branch_length_mode='marginal', max_iter=10,
                         infer_gtr=True, site_specific_gtr=True, pc=pc)

    # calculate likelihood
    tt.infer_ancestral_sequences(marginal=True)
    model = tt.gtr
    np.fill_diagonal(model.W,0)
    acc = {'LH':tt.sequence_LH(), 'tree_length': tt.tree.total_branch_length()}
    acc.update(assess_reconstruction(model, true_model))

    return acc


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    parser.add_argument("--output", type=str, help="folder to save data")
    parser.add_argument("--aa", action='store_true', help="use amino acid alphabet")
    parser.add_argument("--pc", type=float, default=1.0, help="pseudocount for reconstruction")
    parser.add_argument("--files", nargs="+", type=str, help="alignments to analyze")
    args=parser.parse_args()

    output = args.output
    files = sorted(args.files)
    results = []

    pc = args.pc
    analysis_types_tree = [
        'unspecific', # estimate a single model for all sites
        'ml_reconstruction', # reconstruct counts using JC model
        'marginalize', # use expected number of counts
        'iterative', # iterate model estimation and ancestral reconstruction
        'optimize_tree', #estimate model and branch length interatively
    ]

    alphabet='aa_nogap' if args.aa else 'nuc_nogap'

    for fname in files:
        print(fname)
        data_dir = os.path.dirname(fname)
        params = parse_alignment_name(fname)
        params['pc'] = pc

        # load true model, and aln
        true_model = load_model(model_name(data_dir, params))
        true_model_average_rate = true_model.average_rate().mean()
        true_mut_counts = load_mutation_count(mutation_count_name(data_dir, params))

        with gzip.open(alignment_name(data_dir, params), 'rt') as fh:
            aln = AlignIO.read(fh, 'fasta')

        true_LH = get_LH(tree_name(data_dir, params), aln=aln, gtr=true_model)
        true_tree_length = Phylo.read(tree_name(data_dir, params), 'newick').total_branch_length()*params['m']

        ## tree less methods
        # alignment frequencies
        naive = p_from_aln(aln, alphabet=alphabet)
        tmp = {'method':'naive',
               'chisq_p':np.mean([chisq(naive[:,i],true_model.Pi[:,i]) for i in range(params['L'])]),
               'model_entropy':-np.mean(np.sum(naive*np.log(naive+1e-10), axis=0)),
               'true_entropy' :-np.mean(np.sum(true_model.Pi*np.log(true_model.Pi), axis=0))}
        tmp.update(params)
        results.append(tmp)

        # use true counts, no tree
        model = estimate_GTR(true_mut_counts, pc=pc, single_site=False, alphabet=alphabet)
        accuracy = assess_reconstruction(model, true_model)
        accuracy['method'] = 'dressed'
        accuracy.update(params)
        results.append(accuracy)

        for tree in [tree_name(data_dir, params), reconstructed_tree_name(data_dir, params)]:
            for ana in analysis_types_tree:
                accuracy = analyze(ana, tree, aln, alphabet, data_dir, params, true_model)
                s = '_true' if tree == tree_name(data_dir, params) else ''
                ana_t = ana+s
                accuracy['method'] = ana_t
                accuracy['true_LH'] = true_LH
                accuracy['true_tree_length'] = true_tree_length
                accuracy.update(params)
                results.append(accuracy)


    out_dir = os.path.dirname(output) or '.'
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    df = pd.DataFrame(results)
    df.to_csv(output, sep='\t')
