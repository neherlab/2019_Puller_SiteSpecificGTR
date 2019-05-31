import glob
from matplotlib import pyplot as plt
from Bio import Phylo, AlignIO
from treetime import TreeAnc
from treetime.gtr_site_specific import GTR_site_specific


from generate_toy_data import *
from filenames import *
from estimation import *
from reconstruct_toy_data import *


def run_tree(fname, out_prefix, alphabet, model):
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
                 alphabet=alphabet, verbose=1)

    tt.optimize_tree(branch_length_mode='marginal', max_iter=n_iter,
                     infer_gtr=False)

    return tt

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    parser.add_argument("-m", type=float, help="simulated mutation rate")
    parser.add_argument("-n", type=int, help="number of taxa")
    parser.add_argument("-L", type=int, help="length of sequence")
    parser.add_argument("--aa",action='store_true', help="assume amino acid alphabet")
    parser.add_argument("--prefix", type=str, required=True, help="folder to retrieve data from")
    args=parser.parse_args()

    prefix = args.prefix
    mask = "/L{L}_n{n}_m{mu}_*fasta.gz".format(L=args.L or '*', n=args.n or '*', mu=args.m or "*")
    files = glob.glob(prefix+mask)

    n_iter = 20
    alphabet = 'aa_nogap' if args.aa else 'nuc_nogap'
    out_prefix = prefix+'_results_model_dev'


    ttl = {'mu':[], 'pi':[]}
    eps_vals = list(np.linspace(0,.09,10)) + list(np.linspace(0.1,1,10))

    for fname in files:
        params = parse_alignment_name(fname)

        true_model = load_model(model_name(prefix, params))
        n = true_model.n_states

        tree = Phylo.read(tree_name(prefix, params), 'newick')
        true_ttl = tree.total_branch_length()*params['m']

        total_branch_length_pi = []
        total_branch_length_mu = []

        for eps in eps_vals:
            model = GTR_site_specific.custom(mu = true_model.mu/params['m'], W=true_model.W, pi=(true_model.Pi+eps)/1+n*eps, alphabet=alphabet)
            tt = run_tree(fname, out_prefix, alphabet, model)
            print(model.average_rate().mean())
            total_branch_length_pi.append((eps, tt.tree.total_branch_length(), true_ttl, model.average_rate().mean()))

        model = GTR_site_specific.custom(mu = true_model.mu/params['m'], W=true_model.W,
                                         pi=np.ones_like(true_model.Pi, dtype=float)/n, alphabet=alphabet)
        tt = run_tree(fname, out_prefix, alphabet, model)
        total_branch_length_pi.append((np.inf, tt.tree.total_branch_length(), true_ttl, model.average_rate()))

        for eps in eps_vals:
            model = GTR_site_specific.custom(mu = (true_model.mu/params['m'] + eps)/(1+eps), W=true_model.W, pi=true_model.Pi, alphabet=alphabet)
            tt = run_tree(fname, out_prefix, alphabet, model)
            total_branch_length_mu.append((eps, tt.tree.total_branch_length(),true_ttl))

        model = GTR_site_specific.custom(mu = np.ones_like(true_model.mu), W=true_model.W,
                                         pi=true_model.Pi, alphabet=alphabet)
        tt = run_tree(fname, out_prefix, alphabet, model)
        total_branch_length_mu.append((np.inf, tt.tree.total_branch_length(),true_ttl))

        ttl['pi'].append(total_branch_length_pi)
        ttl['mu'].append(total_branch_length_mu)


    for x in ttl['pi']:
        c=np.array(x)
        plt.plot(c[:,0], c[:,1]/c[:,2])

    for x in ttl['mu']:
        c=np.array(x)
        plt.plot(c[:,0], c[:,1]/c[:,2])