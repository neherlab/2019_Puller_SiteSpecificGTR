import glob
from matplotlib import pyplot as plt
from Bio import Phylo, AlignIO
from treetime import TreeAnc
from treetime.gtr_site_specific import GTR_site_specific


from generate_toy_data import *
from filenames import *
from estimation import *
from reconstruct_toy_data import *


def run_tree(fname, alphabet, model, n_iter=20):
    params = parse_alignment_name(fname)
    m = params['m']
    tree = Phylo.read(tree_name(prefix, params), 'newick')
    tree.root.branch_length = 0.001
    tree.ladderize()
    for n in tree.find_clades():
        n.branch_length *= m*(0.6+0.4*np.random.random())

    with gzip.open(alignment_name(prefix, params), 'rt') as fh:
        aln = AlignIO.read(fh, 'fasta')

    tt = TreeAnc(tree=tree, aln=aln, gtr = model, compress=False,
                 alphabet=alphabet, verbose=3)

    tt.optimize_tree(branch_length_mode='marginal', max_iter=n_iter,
                     infer_gtr=False)

    return tt

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    parser.add_argument("--files", nargs="+", type=str, help="alignments to analyze")
    parser.add_argument("--output", type=str, help="folder to save data")
    parser.add_argument("--aa", action='store_true', help="assume amino acid alphabet")
    args=parser.parse_args()

    try:
        dirname = os.path.dirname(args.output)
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
    except:
        pass


    alphabet = 'aa_nogap' if args.aa else 'nuc_nogap'
    res = {}
    eps_vals = list(np.linspace(0.0,1,11))

    for fname in args.files:
        print(fname)
        prefix = os.path.dirname(fname)
        params = parse_alignment_name(fname)
        res[fname] = params

        true_model = load_model(model_name(prefix, params))
        n = true_model.n_states

        tree = Phylo.read(tree_name(prefix, params), 'newick')
        true_ttl = tree.total_branch_length()*params['m']
        res[fname]['true_value'] = true_ttl

        for eps in eps_vals:
            model = load_model(model_name(prefix, params), flatten_p=eps)
            model.mu /= model.average_rate().mean()
            tt = run_tree(fname, alphabet, model)
            res[fname][f'eps_pi_{eps:1.2f}'] = tt.tree.total_branch_length()

        for eps in eps_vals:
            model = load_model(model_name(prefix, params), flatten_mu=eps)
            model.mu /= model.average_rate().mean()
            tt = run_tree(fname, alphabet, model)
            res[fname][f'eps_mu_{eps:1.2f}'] = tt.tree.total_branch_length()

        for eps in eps_vals:
            model = load_model(model_name(prefix, params), flatten_mu=eps, flatten_p=eps, flatten_W=eps)
            model.mu /= model.average_rate().mean()
            tt = run_tree(fname, alphabet, model)
            res[fname][f'eps_all_{eps:1.2f}'] = tt.tree.total_branch_length()

    df = pd.DataFrame(res).T
    df.to_csv(args.output, sep='\t')

