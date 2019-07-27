import glob
from matplotlib import pyplot as plt
from Bio import Phylo, AlignIO
from treetime import TreeAnc

from generate_toy_data import *
from filenames import *
from estimation import *
from reconstruct_toy_data import *



def run_tree(fname, out_prefix, alphabet, true_tree=False, true_model=False, pc=0.1, true_rates=False):
    params = parse_alignment_name(fname)
    m = params['m']
    tree = Phylo.read(tree_name(prefix, params), 'newick') if true_tree else  Phylo.read(reconstructed_tree_name(prefix, params), 'newick')
    tree.root.branch_length = 0.001
    tree.ladderize()
    old_bl = []
    print(np.mean([x for c,x in tree.depths().items() if c.is_terminal()])*(m if true_tree else 1.0))
    print(tree.root.clades[0].branch_length/tree.root.clades[1].branch_length)
    for n in tree.find_clades():
        old_bl.append(n.branch_length)
        if true_tree:
            n.branch_length *= m*(0.6+0.4*np.random.random())
    print(np.sum(old_bl)*(m if true_tree else 1.0),m)
    true_GTR = load_model(model_name(prefix, params))
    if true_model:
        model = true_GTR
        model.mu/=m
    else:
        model='JC69'

    with gzip.open(alignment_name(prefix, params), 'rt') as fh:
        aln = AlignIO.read(fh, 'fasta')

    tt = TreeAnc(tree=tree, aln=aln, gtr = model, compress=False,
                 alphabet=alphabet, verbose=3)

    tt.optimize_tree(branch_length_mode='marginal', max_iter=n_iter,
                     infer_gtr=not true_model, site_specific_gtr=True, pc=pc, damping=0.75)

    if true_rates:
        tt.gtr.mu = true_GTR.mu/m
        tt.optimize_tree(branch_length_mode='marginal', max_iter=n_iter,
                         infer_gtr=False, site_specific_gtr=True, pc=pc, damping=0.75)


    new_bl = []
    for n in tt.tree.find_clades():
        new_bl.append(n.branch_length)


    print(tt.tree.total_branch_length(),tt.gtr.average_rate().mean())
    tt.tree.root_at_midpoint()
    tfname = reoptimized_tree_true_model(out_prefix, params) if args.true_model else reoptimized_tree(out_prefix, params, true_rates=true_rates)
    Phylo.write(tt.tree, tfname, 'newick')

    print(np.mean([x for c,x in tt.tree.depths().items() if c.is_terminal()]), tt.tree.total_branch_length())
    print(tt.tree.root.clades[0].branch_length/tt.tree.root.clades[1].branch_length)
    print(np.corrcoef(old_bl, new_bl)[0,1])


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    parser.add_argument("-m", type=float, help="simulated mutation rate")
    parser.add_argument("-n", type=int, help="number of taxa")
    parser.add_argument("-L", default=1000, type=int, help="length of sequence")
    parser.add_argument("-t", type=int, help="index of tree")
    parser.add_argument("-s", type=int, help="index of sequence")
    parser.add_argument("-g", type=int, help="index of model")
    parser.add_argument("--pc", default=0.1, type=float, help="pseudo count")
    parser.add_argument("--aa", action='store_true', help="assume amino acid alphabet")
    parser.add_argument("--prefix", type=str, required=True, help="folder to retrieve data from")
    parser.add_argument("--true-model", action='store_true', help='assume true model')
    parser.add_argument("--true-rates", action='store_true', help='assume true rates')
    parser.add_argument("--true-tree", action='store_true', help='assume true tree')
    args=parser.parse_args()

    prefix = args.prefix
    mask = "/L{L}_n{n}_m{mu}_*tree{t:02d}_model{g:02d}_seqgen{s:02d}.fasta.gz".format(L=args.L or '*', n=args.n or '*', mu=args.m or "*",
                                                                          t='*' if args.t is None else args.t,
                                                                          g='*' if args.g is None else args.g,
                                                                          s='*' if args.s is None else args.s)
    files = glob.glob(prefix+mask)

    n_iter = 20
    alphabet = 'aa_nogap' if args.aa else 'nuc_nogap'
    if args.true_model:
        out_prefix = prefix+'_results'
    else:
        out_prefix = prefix+'_results_pc_%1.2f'%args.pc

    try:
        if not os.path.isdir(out_prefix):
            os.mkdir(out_prefix)
    except:
        pass

    for fname in files:
        print(fname)
        run_tree(fname, out_prefix, alphabet, args.true_tree, args.true_model, pc=args.pc, true_rates = args.true_rates)
