import glob
from matplotlib import pyplot as plt
from Bio import Phylo, AlignIO
from treetime import TreeAnc

from generate_toy_data import *
from filenames import *
from estimation import *
from reconstruct_toy_data import *


def optimize_branch_length(tree, aln, model, n_iter=5, damping=0.5):
    tt = TreeAnc(tree=tree, aln=aln, gtr = model, reduce_alignment=False)
    for i in range(n_iter):
        tt.infer_ancestral_sequences(marginal=True)
        for n in tt.tree.find_clades():
            if n.up is None:
                continue
            new_val = tt.optimal_marginal_branch_length(n, tol=1e-8 + 0.01**(i+1))
            n.branch_length = new_val*(1-damping**(i+1)) + n.branch_length*damping**(i+1)

    return tt

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    parser.add_argument("-m", type=float, help="simulated mutation rate")
    parser.add_argument("-n", type=int, help="number of taxa")
    parser.add_argument("-L", type=int, help="length of sequence")
    parser.add_argument("--pc", type=float, help="pseudo count")
    parser.add_argument("--aa",action='store_true', help="assume amino acid alphabet")
    parser.add_argument("--prefix", type=str, required=True, help="folder to retrieve data from")
    args=parser.parse_args()

    prefix = args.prefix
    mask = "/L{L}_n{n}_m{mu}_*fasta.gz".format(L=args.L or '*', n=args.n or '*', mu=args.m or "*")
    files = glob.glob(prefix+mask)

    n_iter = 5
    alphabet = 'aa_nogap' if args.aa else 'nuc_nogap'
    out_prefix = prefix+'_results_pc_%1.2f'%args.pc
    for fname in files:
        print(fname)
        params = parse_alignment_name(fname)
        true_tree = Phylo.read(tree_name(prefix, params), 'newick')
        true_model = load_model(model_name(prefix, params))
        m = true_model.average_rate().mean()
        true_model.mu/=m
        print(true_tree.total_branch_length(),m)
        for n in true_tree.find_clades():
            n.branch_length *= m
        print(true_tree.total_branch_length(),m)

        iq_tree = Phylo.read(reconstructed_tree_name(prefix, params), 'newick')
        with gzip.open(alignment_name(prefix, params), 'rt') as fh:
            aln = AlignIO.read(fh, 'fasta')

        tt = TreeAnc(tree=true_tree, aln=aln, gtr = 'JC69', reduce_alignment=False, verbose=3, alphabet=alphabet)
        tt.optimize_tree(branch_length_mode='marginal', infer_gtr=True, site_specific_gtr=True, pc=args.pc)
        Phylo.write(tt.tree, reoptimized_tree(out_prefix, params), 'newick')
        print(np.mean([chisq(tt.gtr.Pi[:,i],true_model.Pi[:,i]) for i in range(params['L'])]))
        inferred_model = tt.gtr
        #inferred_model.mu = inferred_model.mu/inferred_model.mu.mean()

        #tt = TreeAnc(tree=reconstructed_tree_name(prefix, params),
        #             aln=aln, gtr = inferred_model, reduce_alignment=False, verbose=3)
        #tt.optimize_tree(branch_length_mode='marginal', infer_gtr=False)

        true_tree = Phylo.read(tree_name(prefix, params), 'newick')
        for n in true_tree.find_clades():
            n.branch_length *= m

        tt = TreeAnc(tree=true_tree,
                     aln=aln, gtr = true_model, reduce_alignment=False, verbose=3)
        tt.optimize_tree(branch_length_mode='marginal', infer_gtr=False)
        Phylo.write(tt.tree, reoptimized_tree_true_model(out_prefix, params), 'newick')
        print(np.mean([chisq(tt.gtr.Pi[:,i],true_model.Pi[:,i]) for i in range(params['L'])]))

