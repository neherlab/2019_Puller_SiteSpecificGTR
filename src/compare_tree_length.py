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
    mask = "/L{L}_n{n}_m{mu}_*fasta.gz".format(L=args.L or '*', n=args.n or '*', mu=args.m or "*")
    files = glob.glob(prefix+mask)

    length = []
    depth = []
    for fname in files:
        print(fname)
        params = parse_alignment_name(fname)
        true_tree = Phylo.read(tree_name(prefix, params), 'newick')
        true_model = load_model(model_name(prefix, params))
        aln = AlignIO.read(alignment_name(prefix, params), 'fasta')
        iq_tree = Phylo.read(reconstructed_tree_name(prefix, params), 'newick')

        iq_tree_ttl = 0.5*iq_tree.total_branch_length()/params['n']
        iq_tree_depth = np.mean([x for c,x in iq_tree.depths().items() if c.is_terminal()])

        m = true_model.average_rate().mean()
        true_model.mu/=m
        tt = TreeAnc(tree=iq_tree, aln=aln, gtr = true_model, reduce_alignment=False)
        tt.infer_ancestral_sequences(marginal=True)
        for i in range(5):
            x = []
            for n in tt.tree.find_clades():
                if n.up is None:
                    continue
                x.append([tt.optimal_marginal_branch_length(n), n.branch_length])
                n.branch_length = 0.5*(x[-1][0] + n.branch_length)

            x = np.array(x)
            print(x.sum(axis=0))
            tt.infer_ancestral_sequences(marginal=True)
            print(tt.sequence_LH())

        length.append([0.5*true_tree.total_branch_length()*m/params['n'],
                       0.5*tt.tree.total_branch_length()/params['n'], iq_tree_ttl])

        depth.append([np.mean([x*m for c,x in true_tree.depths().items() if c.is_terminal()]),
                      np.mean([x for c,x in tt.tree.depths().items() if c.is_terminal()]), iq_tree_depth])


    length = np.array(length)
    depth = np.array(depth)

    plt.figure(1)
    plt.scatter(length[:,0], length[:,1])
    plt.scatter(length[:,0], length[:,2]*0.5)
    plt.plot([0,.25], [0,.25])


    plt.figure(2)
    plt.scatter(depth[:,0], depth[:,1])
    plt.scatter(depth[:,0], depth[:,2])
    plt.plot([0,.5], [0,.5])
