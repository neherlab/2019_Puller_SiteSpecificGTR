def rtt(T):
    d = T.depths()
    return np.mean([v for k,v in d.items() if k.is_terminal()])

true_model = load_model(model_name(prefix, params), flatten_mu=0.0)
true_mut_counts = load_mutation_count(mutation_count_name(prefix, params))

for f in [0.0, 0.05, 0.1, 0.2, 0.5, 0.75, 1.0]:
    model = load_model(model_name(prefix, params), flatten_mu=f, flatten_p=f)
    m = model.average_rate().mean()
    model.mu/=m
    print(np.abs(model.Pi - true_model.Pi).sum()/1000, np.sum(model.Pi**2, axis=0).mean())

    true_tree = Phylo.read(tree_name(prefix, params), 'newick')
    for n in true_tree.find_clades():
        n.branch_length *= m

    trtt = rtt(true_tree)
    tt = TreeAnc(tree=true_tree,
                 aln=aln, gtr = model, reduce_alignment=False, verbose=0)
    tt.optimize_tree(branch_length_mode='marginal', infer_gtr=False)
    print(f,tt.tree.total_branch_length(), rtt(tt.tree),trtt)


for f in [0.0, 0.05, 0.1, 0.2, 0.5, 0.75, 1.0]:
    model = load_model(model_name(prefix, params), flatten_p=f, flatten_mu=f, flatten_W=0.0)
    m = model.average_rate().mean()
    model.mu/=m

    iq_tree = Phylo.read(reconstructed_tree_name(prefix, params), 'newick')
    iq_tree.root_at_midpoint()
    bl = iq_tree.total_branch_length()
    iqrtt = rtt(iq_tree)

    tt = TreeAnc(tree=iq_tree,
                 aln=aln, gtr = model, reduce_alignment=False, verbose=3)
    tt.optimize_tree(branch_length_mode='marginal', infer_gtr=False)
    print(f,tt.tree.total_branch_length(), bl, rtt(tt.tree),iqrtt)


tt = TreeAnc(tree=true_tree, aln=aln, gtr = 'JC69', reduce_alignment=False, verbose=3, alphabet=alphabet)
tt.optimize_tree(branch_length_mode='marginal', infer_gtr=True, site_specific_gtr=True, pc=args.pc, max_iter=10)
print(f,tt.tree.total_branch_length())

print(np.abs(true_model.Pi - tt.gtr.Pi).sum()/1000, np.sum(tt.gtr.Pi**2, axis=0).mean())


tt.optimize_tree(branch_length_mode='marginal', infer_gtr=False, max_iter=10)


mu_dev = tt.gtr.mu/tt.gtr.mu.mean() - true_model.mu/true_model.mu.mean()

p_dev = np.sum((tt.gtr.Pi - true_model.Pi)**2, axis=0)

p_dev = 1.0/np.sum(tt.gtr.Pi**2, axis=0) - 1.0/np.sum(true_model.Pi**2, axis=0)
