import glob
from matplotlib import pyplot as plt
from Bio import Phylo, AlignIO
from treetime import TreeAnc

from generate_toy_data import *
from filenames import *
from estimation import *
from reconstruct_toy_data import *



def run_tree(fname, out_prefix, alphabet, true_tree=False,
			 true_model=False, pc=0.1, true_rates=False):
    """
    read a tree and an alignment and optimize its branch length using different types of models.
    the use can specify to either use the true model for optimization, just the true rates, or
    infer the entire model from the data. The either the true or an inferred tree-topology can
    be used.
    """
    params = parse_alignment_name(fname)
    params['pc'] = pc
    prefix = os.path.dirname(fname)
    m = params['m']
    tree = Phylo.read(tree_name(prefix, params), 'newick') if true_tree else  Phylo.read(reconstructed_tree_name(prefix, params), 'newick')
    tree.root.branch_length = 0.001


    old_bl = []
    print(np.mean([x for c,x in tree.depths().items() if c.is_terminal()])*(m if true_tree else 1.0))
    print(tree.root.clades[0].branch_length/tree.root.clades[1].branch_length)

    # randomize branch length of true tree to allow fair comparison
    for n in tree.find_clades():
        old_bl.append(n.branch_length)
        if true_tree:
            # rescale with mutation rate and multiply by a random number between 0.6 and 1.0
            n.branch_length *= m*(0.6+0.4*np.random.random())

    print(np.sum(old_bl)*(m if true_tree else 1.0),m)
    # load true GTR model. Use this for inference if true_tree=True, else start with Jukes Cantor
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

    # run the tree optimization of treetime. the damping parameter slows down the iterative
    # branch length optimization to avoid oscillations and run-away solutions
    # a site-specific GTR model is inferred if true_model is False
    tt.optimize_tree(branch_length_mode='marginal', max_iter=n_iter,
                     infer_gtr=not true_model, site_specific_gtr=True, pc=pc, damping=0.75)

    # if the true raes are to be used, replace those in the model and re-optimize
    if true_rates:
        tt.gtr.mu = true_GTR.mu/m
        tt.optimize_tree(branch_length_mode='marginal', max_iter=n_iter,
                         infer_gtr=False, site_specific_gtr=True, pc=pc, damping=0.75)


    new_bl = []
    for n in tt.tree.find_clades():
        new_bl.append(n.branch_length)


    # save new tree to file
    tt.tree.root_at_midpoint()
    tfname = reoptimized_tree_true_model(out_prefix, params) if args.true_model else reoptimized_tree(out_prefix, params, true_rates=true_rates)
    Phylo.write(tt.tree, tfname, 'newick')

    print(tt.tree.total_branch_length(),tt.gtr.average_rate().mean())
    print(np.mean([x for c,x in tt.tree.depths().items() if c.is_terminal()]), tt.tree.total_branch_length())
    print(tt.tree.root.clades[0].branch_length/tt.tree.root.clades[1].branch_length)
    print(np.corrcoef(old_bl, new_bl)[0,1])


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    parser.add_argument("--files", nargs="+", type=str, help="alignments to analyze")
    parser.add_argument("--output", type=str, help="folder to save data")
    parser.add_argument("--pc", default=0.1, type=float, help="pseudo count")
    parser.add_argument("--aa", action='store_true', help="assume amino acid alphabet")
    parser.add_argument("--true-model", action='store_true', help='assume true model')
    parser.add_argument("--true-rates", action='store_true', help='assume true rates')
    parser.add_argument("--true-tree", action='store_true', help='assume true tree')
    args=parser.parse_args()

    n_iter = 20
    alphabet = 'aa_nogap' if args.aa else 'nuc_nogap'

    try:
        if not os.path.isdir(args.output):
            os.mkdir(args.output)
    except:
        pass

    for fname in args.files:
        print(fname)
        run_tree(fname, args.output, alphabet, args.true_tree, args.true_model,
        		 pc=args.pc, true_rates = args.true_rates)
