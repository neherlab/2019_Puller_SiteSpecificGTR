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
import cProfile, pstats, io

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
    parser.add_argument("--prefix", type=str, help="folder to save data")
    parser.add_argument("--aa", action='store_true', help="use amino acid alphabet")
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

    pc=0.1
    niter=5

    analysis_types = ['naive', 'single', 'dressed', 'regular', 'regular_true',
                      'marginal','marginal_true', 'iterative', 'iterative_true']

    analysis_types = ['iterative']

    alphabet='aa_nogap' if args.aa else 'nuc_nogap'
    for fname in files:
        print(fname)

        params = parse_alignment_name(fname)
        mu_vals.add(params['m'])
        n_vals.add(params['n'])

        true_model = load_model(model_name(prefix, params))
        true_model_average_rate = true_model.average_rate().mean()
        true_mut_counts = load_mutation_count(mutation_count_name(prefix, params))
        dset = (params['L'], params['n'], params['m'])

        _mc, true_LH, true_tree = reconstruct_counts(prefix, params, gtr=true_model, alphabet=alphabet,
                                              marginal=True, reconstructed_tree=False)
        true_tree_length = true_tree.tree.total_branch_length()
        for ana in analysis_types:
            if ana=='naive':
                naive = p_from_aln(prefix, params, alphabet=alphabet)
                p_dist[ana][dset].append(np.mean([chisq(naive[:,i],true_model.Pi[:,i]) for i in range(params['L'])]))
                p_entropy[ana][dset].append([-np.mean(np.sum(naive*np.log(naive+1e-10), axis=0)),
                                            -np.mean(np.sum(true_model.Pi*np.log(true_model.Pi), axis=0))])
            else:
                if ana=='dressed':
                    mc = (true_mut_counts, None, None)
                    bl = [n.branch_length for n in _t.tree.find_clades() if n!=_t]
                elif ana in ['iterative', 'iterative_true']:
                    model = 'JC69'
                    pr = cProfile.Profile()
                    pr.enable()
                    for i in range(niter):
                        mc = reconstruct_counts(prefix, params, gtr=model,
                                                alphabet=alphabet, marginal=True,
                                                reconstructed_tree=ana=='iterative')
                        model = estimate_GTR(mc[0], pc=pc, single_site=False, alphabet=alphabet)
                        print(i, np.mean([chisq(model.Pi[:,i],true_model.Pi[:,i]) for i in range(params['L'])]))
                    pr.disable()

                    s = io.StringIO()

                    ps = pstats.Stats(pr, stream=s).sort_stats('cumtime')
                    ps.print_stats()
                    print(s.getvalue())
                elif ana=='branch_length':
                    model = 'JC69'
                    kappa=1.0
                    print(true_model.mu[:5])
                    for i in range(niter):
                        mc = reconstruct_counts(prefix, params, gtr=model,
                                                alphabet=alphabet, marginal=True,
                                                reconstructed_tree=False)
                        model = estimate_GTR(mc[0], pc=pc, single_site=False, tt=mc[-1] if i else None, alphabet=alphabet)
                        mc[-1].set_gtr(model)
                else:
                    rec_model = {'true_model':true_model}
                    mc = reconstruct_counts(prefix, params, gtr=rec_model.get(ana, 'JC69'),
                                            alphabet=alphabet, marginal=ana in ['marginal', 'true_model', 'iterative'],
                                            reconstructed_tree=('true' not in ana))

                if ana not in ['iterative', 'branch_length']:
                    model = estimate_GTR(mc[0], pc=pc, single_site=ana=='single', alphabet=alphabet)
                avg_rate[ana][dset].append((true_model_average_rate, model.average_rate().mean()))
                _mc, model_LH, _t = reconstruct_counts(prefix, params, gtr=model, alphabet=alphabet,
                                                       marginal=True, reconstructed_tree=True)
                delta_LH[ana][dset].append( (true_LH, model_LH) )

                if ana!='single':
                    p_dist[ana][dset].append(np.mean([chisq(model.Pi[:,i],true_model.Pi[:,i]) for i in range(params['L'])]))
                    p_entropy[ana][dset].append([-np.mean(np.sum(model.Pi*np.log(model.Pi), axis=0)),
                                                -np.mean(np.sum(true_model.Pi*np.log(true_model.Pi), axis=0))])
                    mu_dist[ana][dset].append(chisq(model.mu/model.average_rate().mean(), true_model.mu/true_model_average_rate))

                np.fill_diagonal(model.W,0)
                W_dist[ana][dset].append(chisq(model.W.flatten(), true_model.W.flatten()))

    out_prefix = prefix + '_results_pc_%1.1f/'%pc
    if not os.path.isdir(out_prefix):
        os.mkdir(out_prefix)

    out_fname = out_prefix + "_".join(["{name}{val}".format(name=n, val=args.__getattribute__(n)) for n in ['L', 'n', 'm'] if args.__getattribute__(n)]) + '.pkl'
    with open(out_fname, 'wb') as fh:
        pickle.dump((sorted(mu_vals), sorted(n_vals), dict(p_dist), dict(p_entropy), dict(mu_dist),
                     dict(W_dist), dict(delta_LH), dict(avg_rate)), fh)
