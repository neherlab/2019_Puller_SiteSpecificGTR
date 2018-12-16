import os

def parse_alignment_name(fname):
    base = os.path.basename(fname)[:-9]
    params = {}
    for x in base.split('_'):
        try:
            params[x[0]]=int(x[1:])
        except:
            try:
                params[x[0]]=float(x[1:])
            except:
                try:
                    params[x[:-2]]=int(x[-2:])
                except:
                    pass

    return params


def model_name(prefix, params):
    return prefix+'/L{L}_n{n}_m{m}_tree{tree:02d}_model{model:02d}.npz'.format(**params)

def mutation_count_name(prefix, params):
    return prefix+'/L{L}_n{n}_m{m}_tree{tree:02d}_model{model:02d}_seqgen{seqgen:02}_mutations.npz'.format(**params)

def alignment_name(prefix, params):
    return prefix+'/L{L}_n{n}_m{m}_tree{tree:02d}_model{model:02d}_seqgen{seqgen:02}.fasta.gz'.format(**params)

def tree_name(prefix, params):
    return prefix+'/L{L}_n{n}_m{m}_tree{tree:02d}.nwk'.format(**params)

def reconstructed_tree_name(prefix, params):
    return prefix+'/L{L}_n{n}_m{m}_tree{tree:02d}_model{model:02d}_seqgen{seqgen:02}_reconstructed.nwk'.format(**params)
