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
    parser.add_argument("-L", type=int, help="length of sequence")
    parser.add_argument("--prefix", type=str, required=True, help="folder to retrieve data from")
    args=parser.parse_args()

    prefix = args.prefix
    result_prefix = prefix+'_results_pc_0.1'
    mask = "/L{L}_n{n}_m{mu}_*reoptimized_true_model.nwk".format(L=args.L or '*', n=args.n or '*', mu=args.m or "*")
    files = glob.glob(result_prefix+mask)

    length = []
    depth = []
    for fname in files:
        print(fname)
        params = parse_alignment_name(fname)
        true_tree = Phylo.read(tree_name(prefix, params), 'newick')
        true_model = load_model(model_name(prefix, params))
        m = true_model.average_rate().mean()
        for n in true_tree.find_clades():
            n.branch_length *= m
        
        iq_tree = Phylo.read(reconstructed_tree_name(prefix, params), 'newick')
        reoptimized_T = Phylo.read(reoptimized_tree(result_prefix, params), 'newick')
        reoptimized_T_true_model = Phylo.read(reoptimized_tree_true_model(result_prefix, params), 'newick')

        tmp_ttl = []
        tmp_depth = []
        for t in [true_tree, iq_tree, reoptimized_T, reoptimized_T_true_model]:
            tmp_ttl.append(0.5*t.total_branch_length()/params['n'])
            tmp_depth.append(np.mean([x for c,x in t.depths().items() if c.is_terminal()]))
            

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

    plt.figure(2)
    plt.scatter(depth[:,0], depth[:,1], label='iqtree')
    plt.scatter(depth[:,0], depth[:,2], label='inferred model') 
    plt.scatter(depth[:,0], depth[:,3], label='true model') 
    plt.plot([0,3.5], [0,3.5])
    plt.legend()
    
