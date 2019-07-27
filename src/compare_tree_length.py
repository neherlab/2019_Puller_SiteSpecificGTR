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
    parser.add_argument("--true-rates", action='store_true', help="include estimates with true rates")
    parser.add_argument("--pc", nargs='+', default=[0.1],type=float, help="pseudo count")
    parser.add_argument("-L", type=int, help="length of sequence")
    parser.add_argument("--prefix", type=str, required=True, help="folder to retrieve data from")
    args=parser.parse_args()

    prefix = args.prefix
    aa ='_aa' in args.prefix

    length = {}
    depth = {}
    terminal_bl = {}

    for dset in ['true', 'ML']+[(pc, 'inferred-rates') for pc in args.pc]:
        length[dset]=[]
        depth[dset]=[]
        terminal_bl[dset]=[]

    if args.true_rates:
        for dset in [(pc, 'true-rates') for pc in args.pc]:
            length[dset]=[]
            depth[dset]=[]
            terminal_bl[dset]=[]


    mask = "/L{L}_n{n}_m{mu}_*.fasta.gz".format(L=args.L or '*', n=args.n or '*', mu=args.m or "*")
    files = glob.glob(prefix+mask)
    for fname in files:
        print(fname)
        params = parse_alignment_name(fname)
        true_tree = Phylo.read(tree_name(prefix, params), 'newick')
        #true_model = load_model(model_name(prefix, params))
        for n in true_tree.find_clades():
            n.branch_length *= params['m']
        true_ttl = 0.5*true_tree.total_branch_length()/params['n']
        true_depth = np.mean([x for c,x in true_tree.depths().items() if c.is_terminal()])
        true_terminal = np.mean([x.branch_length for x in true_tree.get_terminals()])
        true_rates = False
        for dset in depth.keys():
            result_prefix = args.prefix+'_results'
            if type(dset)==tuple:
                result_prefix += '_pc_%1.2f'%dset[0]
                true_rates = dset[1]=='true-rates'

            if dset=='ML':
                tree = Phylo.read(reconstructed_tree_name(prefix, params), 'newick')
            elif dset=='true':
                tree = Phylo.read(reoptimized_tree_true_model(result_prefix, params), 'newick')
            else:
                tree = Phylo.read(reoptimized_tree(result_prefix, params, true_rates=true_rates), 'newick')

            tree.root_at_midpoint()
            tree.root.branch_length=0.001

            terminal_bl[dset].append([true_terminal, np.mean([x.branch_length for x in tree.get_terminals()])])
            length[dset].append([true_ttl, 0.5*tree.total_branch_length()/params['n']])
            depth[dset].append([true_depth, np.mean([x for c,x in tree.depths().items() if c.is_terminal()])])


    for dset in length:
        length_a = np.array(length[dset])
        depth_a = np.array(depth[dset])
        terminal_a = np.array(terminal_bl[dset])
        if dset=='ML':
            label =  'FastTree JTT+CAT20' if aa else 'IQ-tree GTR+R10'
        elif dset=='true':
            label='True model'
        else:
            label=f"Inferred model {'(true rates)' if dset[1]=='true-rates' else ''}, pc={dset[0]:1.1f}"
        plt.figure(1)
        plt.scatter(length_a[:,0], length_a[:,1], label=label)
        plt.figure(2)
        plt.scatter(depth_a[:,0], depth_a[:,1], label=label)
        plt.figure(3)
        plt.scatter(terminal_a[:,0], terminal_a[:,1], label=label)

    plt.figure(1)
    plt.plot([0,0.25], [0,0.25])
    plt.ylabel('inferred average branch length')
    plt.xlabel('true average branch length')
    plt.legend()
    plt.savefig("figures/"+prefix+"_length.pdf")

    plt.figure(2)
    plt.plot([0,3.5], [0,3.5])
    plt.ylabel('inferred average root-to-tip distance')
    plt.xlabel('true average root-to-tip distance')
    plt.legend()
    plt.savefig("figures/"+prefix+"_depth.pdf")

    plt.figure(3)
    plt.plot([0,0.25], [0,0.25])
    plt.ylabel('inferred average terminal branch length')
    plt.xlabel('true average terminal branch length')
    plt.legend()
    plt.savefig("figures/"+prefix+"_terminal.pdf")
    # for dset in all_bl:
    #     plt.figure()
    #     plt.title(str(dset))
    #     for x in all_bl[dset]:
    #         plt.scatter(x[:,0], x[:,1])
    #     plt.plot([0,1],[0,1])
    #     plt.ylim([0,1.5])
    #     plt.xlim([0,1.5])
