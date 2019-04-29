import os, gzip, glob, pickle
import numpy as np
from collections import defaultdict

from treetime.treeanc import TreeAnc
from treetime.gtr_site_specific import GTR_site_specific
from treetime.gtr import GTR
from treetime.seq_utils import seq2prof, profile_maps, alphabets

from generate_toy_data import *
from filenames import *
from estimation import *

def KL(p,q):
    return np.sum(p*log(p/q))

def chisq(p,q):
    return np.sum((p-q)**2)

def assess_reconstruction(model, true_model):
    return {"chisq_p": np.mean([chisq(model.Pi[:,i],true_model.Pi[:,i]) for i in range(params['L'])]),
            "model_entropy": -np.mean(np.sum(model.Pi*np.log(model.Pi), axis=0)),
            "true_entropy": -np.mean(np.sum(true_model.Pi*np.log(true_model.Pi), axis=0)),
            "chisq_mu": chisq(model.mu/model.average_rate().mean(), true_model.mu/true_model_average_rate),
            "chisq_W":chisq(model.W.flatten(), true_model.W.flatten())}


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    parser.add_argument("-m", type=float, help="simulated mutation rate")
    parser.add_argument("-n", type=int, help="number of taxa")
    parser.add_argument("-L", type=int, help="length of sequence")
    parser.add_argument("--prefix", type=str, help="folder to save data")
    parser.add_argument("--aa", action='store_true', help="use amino acid alphabet")
    parser.add_argument("--pc", type=float, default=1.0, help="pseudocount for reconstruction")
    args=parser.parse_args()

    prefix = args.prefix
    mask = "/L{L}_n{n}_m{mu}_*fasta.gz".format(L=args.L or '*', n=args.n or '*', mu=args.m or "*")
    files = glob.glob(prefix+mask)

    mu_dist = defaultdict(lambda: defaultdict(list))
    p_dist = defaultdict(lambda: defaultdict(list))
    p_entropy = defaultdict(lambda: defaultdict(list))
    W_dist = defaultdict(lambda: defaultdict(list))
    delta_LH = defaultdict(lambda: defaultdict(list))
    avg_rate = defaultdict(lambda: defaultdict(list))

    mu_vals = set()
    n_vals = set()

    pc = args.pc
    analysis_types_tree = [
        'single', # estimate a single model for all sites
        'regular', # reconstruct counts using JC model
        'marginal', # use expected number of counts
        'iterative', # iterate model estimation and reconstruction
        'optimize_tree', #estimate model and branch length interatively
    ]

    alphabet='aa_nogap' if args.aa else 'nuc_nogap'

    for fname in files:
        print(fname)

        params = parse_alignment_name(fname)
        dset = (params['L'], params['n'], params['m'])
        mu_vals.add(params['m'])
        n_vals.add(params['n'])

        # load true model, and aln
        true_model = load_model(model_name(prefix, params))
        true_model_average_rate = true_model.average_rate().mean()
        true_mut_counts = load_mutation_count(mutation_count_name(prefix, params))

        with gzip.open(alignment_name(prefix, params), 'rt') as fh:
            aln = AlignIO.read(fh, 'fasta')

        tt_true = TreeAnc(tree=tree_name(prefix, params), aln=aln, gtr=true_model, reduce_alignment=False)
        tt_true.infer_ancestral_sequences(marginal=True)
        true_tree_length = tt_true.tree.total_branch_length()

        ## tree less methods
        # alignment frequencies
        naive = p_from_aln(aln, alphabet=alphabet)
        p_dist["naive"][dset].append(np.mean([chisq(naive[:,i],true_model.Pi[:,i]) for i in range(params['L'])]))
        p_entropy["naive"][dset].append([-np.mean(np.sum(naive*np.log(naive+1e-10), axis=0)),
                                     -np.mean(np.sum(true_model.Pi*np.log(true_model.Pi), axis=0))])

        # use true counts, no tree
        model = estimate_GTR(true_mut_counts, pc=pc, single_site=False, alphabet=alphabet)
        accuracy = assess_reconstruction(true_model, model)
        p_dist["dressed"][dset].append(accuracy["chisq_p"])
        p_entropy["dressed"][dset].append([accuracy["model_entropy"], accuracy["true_entropy"]])
        mu_dist["dressed"][dset].append(accuracy["chisq_mu"])
        W_dist["dressed"][dset].append(accuracy["chisq_W"])
        avg_rate["dressed"][dset].append((true_model_average_rate, model.average_rate().mean()))

        for tree in [tree_name(prefix, params), reconstructed_tree_name(prefix, params)]:
            for ana in analysis_types_tree:
                T = Phylo.read(tree, format="newick")
                s = ''
                if tree == tree_name(prefix, params):
                    s = '_true'
                    for n in T.find_clades():
                        n.branch_length *= params['m']

                tt = TreeAnc(tree=T, aln = aln, reduce_alignment=False, alphabet=alphabet)

                if ana=="iterative":
                    # this performs iterative estimation of the model and ancestral sequences
                    # using marginal ancestral reconstruction
                    tt.infer_gtr_iterative(normalized_rate=False, site_specific=True, pc=pc)
                elif ana=="single":
                    tt.infer_gtr(marginal=False, normalized_rate=False, site_specific=False, pc=pc)
                elif ana=="regular":
                    tt.infer_gtr(marginal=False, normalized_rate=False, site_specific=True, pc=pc)
                elif ana=="marginal":
                    tt.infer_gtr(marginal=True, normalized_rate=False, site_specific=True, pc=pc)
                elif ana=='optimize_tree':
                    tt.optimize_tree(branch_length_mode='marginal', max_iter=10,
                                     infer_gtr=True, site_specific_gtr=True, pc=pc)

                # calculate likelihood
                tt.infer_ancestral_sequences(marginal=True)

                model = tt.gtr
                np.fill_diagonal(model.W,0)
                avg_rate[ana+s][dset].append((true_model_average_rate, model.average_rate().mean()))
                delta_LH[ana+s][dset].append( (tt_true.sequence_LH(), tt.sequence_LH() ))

                if ana!='single':
                    accuracy = assess_reconstruction(true_model, model)
                    p_dist[ana+s][dset].append(accuracy["chisq_p"])
                    p_entropy[ana+s][dset].append([accuracy["model_entropy"], accuracy["true_entropy"]])
                    mu_dist[ana+s][dset].append(accuracy["chisq_mu"])
                    W_dist[ana+s][dset].append(accuracy["chisq_W"])
                else:
                    W_dist[ana][dset].append(chisq(model.W.flatten(), true_model.W.flatten()))


    out_prefix = prefix + '_results_pc_%1.2f/'%pc
    if not os.path.isdir(out_prefix):
        os.mkdir(out_prefix)

    out_fname = out_prefix + "_".join(["{name}{val}".format(name=n, val=args.__getattribute__(n)) for n in ['L', 'n', 'm'] if args.__getattribute__(n)]) + '.pkl'
    with open(out_fname, 'wb') as fh:
        pickle.dump((sorted(mu_vals), sorted(n_vals), dict(p_dist), dict(p_entropy), dict(mu_dist),
                     dict(W_dist), dict(delta_LH), dict(avg_rate)), fh)
