import glob
from matplotlib import pyplot as plt
from Bio import Phylo, AlignIO
from generate_toy_data import *
from filenames import *
from estimation import *
from treetime import TreeAnc

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    parser.add_argument("-m", type=float, help="simulated mutation rate")
    parser.add_argument("-n", type=int, help="number of taxa")
    parser.add_argument("--pc", type=float, help="pseudo count")
    parser.add_argument("-L", type=int, help="length of sequence")
    parser.add_argument("--prefix", type=str, required=True, help="folder to retrieve data from")
    args=parser.parse_args()

    prefix = args.prefix
    result_prefix = prefix+'_results_pc_%1.2f'%args.pc
    mask = "/L{L}_n{n}_m{mu}_*reoptimized.nwk".format(L=args.L or '*', n=args.n or '*', mu=args.m or "*")
    files = glob.glob(result_prefix+mask)

    length = []
    depth = []
    for fname in files:
        params = parse_alignment_name(fname)
        true_tree = Phylo.read(tree_name(prefix, params), 'newick')
        true_model = load_model(model_name(prefix, params))
        cf = true_model.average_rate().mean()/params['m']
        for n in true_tree.find_clades():
            n.branch_length *= params['m']

        iq_tree = Phylo.read(reconstructed_tree_name(prefix, params), 'newick')
        iq_tree.root.branch_length=0.001
        for n in iq_tree.find_clades():
            n.branch_length /= cf
        reoptimized_T = Phylo.read(reoptimized_tree(result_prefix, params), 'newick')
        for n in reoptimized_T.find_clades():
            n.branch_length /= cf
        reoptimized_T_true_model = Phylo.read(reoptimized_tree_true_model(result_prefix, params), 'newick')

        tmp_ttl = []
        tmp_depth = []
        for t in [true_tree, iq_tree, reoptimized_T, reoptimized_T_true_model]:
            t.root_at_midpoint()
            tmp_ttl.append(0.5*t.total_branch_length()/params['n'])
            tmp_depth.append(np.mean([x for c,x in t.depths().items() if c.is_terminal()]))

        if tmp_depth[2]>1.2*tmp_depth[-1]:
            print("odd", fname, tmp_depth)
        length.append(tmp_ttl)
        depth.append(tmp_depth)


    length = np.array(length)
    depth = np.array(depth)

    plt.figure(1)
    plt.scatter(length[:,0], length[:,1], label='iqtree')
    plt.scatter(length[:,0], length[:,2], label='inferred model')
    plt.scatter(length[:,0], length[:,3], label='true model')
    plt.plot([0,.25], [0,.25])
    plt.legend()
    plt.savefig("figures/"+result_prefix+"_length.pdf")

    plt.figure(2)
    plt.scatter(depth[:,0], depth[:,1], label='iqtree')
    plt.scatter(depth[:,0], depth[:,2], label='inferred model')
    plt.scatter(depth[:,0], depth[:,3], label='true model')
    plt.plot([0,3.5], [0,3.5])
    plt.legend()
    plt.savefig("figures/"+result_prefix+"_depth.pdf")
