import os, gzip, glob, pickle
import numpy as np
from collections import defaultdict

from scipy.optimize import minimize

from treetime.treeanc import TreeAnc
from treetime.gtr_site_specific import GTR_site_specific
from treetime.gtr import GTR
from treetime.seq_utils import seq2prof, profile_maps, alphabets

from generate_toy_data import *
from filenames import *
from estimation import *
from reconstruct_toy_data import *

def lh(x,tt):
    for ni, n in enumerate(tt.tree.find_clades()):
        n.branch_length = x[ni]**2

    tt.infer_ancestral_sequences(marginal=True)
    print(tt.sequence_LH(), np.sum(x**2))
    return -tt.sequence_LH()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    parser.add_argument("-m", type=float, help="simulated mutation rate")
    parser.add_argument("-n", type=int, help="number of taxa")
    parser.add_argument("-L", type=int, help="length of sequence")
    parser.add_argument("--prefix", type=str, help="folder to save data")
    args=parser.parse_args()

    prefix = args.prefix
    mask = "/L{L}_n{n}_m{mu}_*fasta.gz".format(L=args.L or '*', n=args.n or '*', mu=args.m or "*")
    files = glob.glob(prefix+mask)

    pc=1
    niter=5

    for fname in files[:1]:
        print(fname)

        params = parse_alignment_name(fname)
        true_model = load_model(model_name(prefix, params))
        true_model_average_rate = true_model.average_rate().mean()
        true_mut_counts = load_mutation_count(mutation_count_name(prefix, params))
        dset = (params['L'], params['n'], params['m'])

        _mc, true_LH, true_tree = reconstruct_counts(prefix, params, gtr=true_model, alphabet='nuc_nogap',
                                              marginal=True, reconstructed_tree=False)
        true_tree_length = true_tree.tree.total_branch_length()

        model = 'JC69'
        for i in range(niter):
            mc = reconstruct_counts(prefix, params, gtr=model,
                                    alphabet='nuc_nogap', marginal=True,
                                    reconstructed_tree=False)
            model = estimate_GTR(mc[0], pc=pc, single_site=False)
            print(i, np.mean([chisq(model.Pi[:,i],true_model.Pi[:,i]) for i in range(params['L'])]))

        tt = mc[-1]
        tt.set_gtr(model)
        tt.infer_ancestral_sequences(marginal=True)
        print(tt.sequence_LH())
        y=np.array([n.branch_length for n in tt.tree.find_clades() if n.up])

        from matplotlib import pyplot as plt
        plt.figure()
        current = []
        for ni,n in enumerate(tt.tree.get_terminals()):
            if n.up is None:
                continue
            if ni%10==0:
                tmp = np.array([[t, tt.gtr.prob_t_profiles(tt.marginal_branch_profile(n), tt.multiplicity, t, return_log=True)]
                            for t in np.linspace(0,2*n.branch_length, 21)])
                current.append([n.branch_length, tt.gtr.prob_t_profiles(tt.marginal_branch_profile(n), tt.multiplicity, n.branch_length, return_log=True)])
                plt.plot(tmp[:,0], tmp[:,1]-current[-1][1])
        current = np.array(current)
        # plt.scatter(current[:,0], current[:,1])

        print(tt.sequence_LH())
        for i in range(10):
            x = []
            for n in tt.tree.find_clades():
                if n.up is None:
                    continue
                x.append([tt.optimal_marginal_branch_length(n), n.branch_length])
                n.branch_length = x[-1][0]

            x = np.array(x)
            print(x.sum(axis=0))
            tt.infer_ancestral_sequences(marginal=True)
            print(tt.sequence_LH())

        # x0 = np.array([np.sqrt(n.branch_length) for n in tt.tree.find_clades()])

        # minimize(lh, x0, args=(tt,), method='BFGS', options={'eps':0.001})
