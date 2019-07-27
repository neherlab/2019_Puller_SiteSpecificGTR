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
                 alphabet=alphabet, verbose=3)

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
    out_prefix = prefix+'_results'
    muvals =  [0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.35, 0.5]
    out_file = out_prefix+'/model_deviation_n{n}_m{m}.pkl'.format(n=args.n, m=args.m)

    try:
        with open(out_file, 'rb') as fh:
            import pickle
            ttl, ttl_flat = pickle.load(fh)
    except:
        ttl_flat = {'pi':[], 'mu':[], 'all':[]}
        ttl = {'mu':[], 'pi':[], 'all':[]}
        eps_vals = list(np.linspace(0.0,1,11))

        for fname in files:
            print(fname)
            params = parse_alignment_name(fname)

            true_model = load_model(model_name(prefix, params))
            n = true_model.n_states

            tree = Phylo.read(tree_name(prefix, params), 'newick')
            true_ttl = tree.total_branch_length()*params['m']

            total_branch_length_pi = []
            total_branch_length_mu = []
            total_branch_length_all = []

            for eps in eps_vals:
                model = load_model(model_name(prefix, params), flatten_p=eps)
                model.mu /= model.average_rate().mean()
                tt = run_tree(fname, out_prefix, alphabet, model)
                total_branch_length_pi.append((eps, tt.tree.total_branch_length(), true_ttl, model.mu.mean()))
                print(total_branch_length_pi[-1])

            for eps in eps_vals:
                model = load_model(model_name(prefix, params), flatten_mu=eps)
                model.mu /= model.average_rate().mean()
                tt = run_tree(fname, out_prefix, alphabet, model)
                total_branch_length_mu.append((eps, tt.tree.total_branch_length(),true_ttl, model.mu.mean()))

            for eps in eps_vals:
                model = load_model(model_name(prefix, params), flatten_mu=eps, flatten_p=eps, flatten_W=eps)
                model.mu /= model.average_rate().mean()
                tt = run_tree(fname, out_prefix, alphabet, model)
                total_branch_length_all.append((eps, tt.tree.total_branch_length(),true_ttl, model.mu.mean()))

            ttl['pi'].append(total_branch_length_pi)
            ttl['mu'].append(total_branch_length_mu)
            ttl['all'].append(total_branch_length_all)

        for mu in muvals:
            mask = "/L{L}_n{n}_m{mu}_*fasta.gz".format(L=args.L or '*', n=args.n or '*', mu=mu)
            files = glob.glob(prefix+mask)
            total_branch_length_pi = []
            total_branch_length_mu = []
            total_branch_length_all = []

            eps=1.0
            for fname in files:
                print(fname)
                params = parse_alignment_name(fname)

                true_model = load_model(model_name(prefix, params))
                n = true_model.n_states

                tree = Phylo.read(tree_name(prefix, params), 'newick')
                true_ttl = tree.total_branch_length()*params['m']

                model = load_model(model_name(prefix, params), flatten_p=eps)
                model.mu /= model.average_rate().mean()
                tt = run_tree(fname, out_prefix, alphabet, model)
                total_branch_length_pi.append((mu, tt.tree.total_branch_length(), true_ttl, model.mu.mean()))
                print(total_branch_length_pi[-1])

                model = load_model(model_name(prefix, params), flatten_mu=eps)
                model.mu /= model.average_rate().mean()
                tt = run_tree(fname, out_prefix, alphabet, model)
                total_branch_length_mu.append((mu, tt.tree.total_branch_length(),true_ttl, model.mu.mean()))

                model = load_model(model_name(prefix, params), flatten_mu=eps, flatten_p=eps, flatten_W=eps)
                model.mu /= model.average_rate().mean()
                tt = run_tree(fname, out_prefix, alphabet, model)
                total_branch_length_all.append((mu, tt.tree.total_branch_length(),true_ttl, model.mu.mean()))

            ttl_flat['pi'].append(total_branch_length_pi)
            ttl_flat['mu'].append(total_branch_length_mu)
            ttl_flat['all'].append(total_branch_length_all)

        if not os.path.isdir(out_prefix):
            os.mkdir(out_prefix)

        with open(out_file, 'wb') as fh:
            import pickle
            pickle.dump((ttl, ttl_flat), fh)


    fig, axs = plt.subplots(1,2, sharey=True, figsize=(13,7))
    fs=16
    for k, label_str in [('pi', 'frequencies'), ('mu', 'rates'), ('all', 'both')]:
        c=np.array([x for x in ttl[k]])
        x = c[0,:,0]
        val = c[:,:,1]/c[:,:,2]
        m = val.mean(axis=0)
        std = val.std(axis=0)
        axs[0].errorbar(x, m, std, label=label_str)
    axs[0].set_xlabel('fraction true', fontsize=fs)
    axs[0].set_ylabel('rel. total branch length error', fontsize=fs)
    axs[0].tick_params(labelsize=fs*0.8)

    for k, label_str in [('pi', 'frequencies'), ('mu', 'rates'), ('all', 'both')]:
        c=np.array([x for x in ttl_flat[k]])
        x = c[:,0,0]
        val = c[:,:,1]/c[:,:,2]
        m = val.mean(axis=1)
        std = val.std(axis=1)
        axs[1].errorbar(x, m, std, label=label_str)
    axs[1].set_xlabel('evolutionary rate', fontsize=fs)
    axs[1].legend(fontsize=fs)
    axs[1].tick_params(labelsize=fs*0.8)

    plt.savefig('figures/model_deviation_{prefix}_n{n}_m{m}.pdf'.format(prefix=args.prefix, n=args.n, m=args.m))
