import matplotlib.pyplot as plt
from generate_toy_data import *

def entropy(a):
    eps = 1e-10
    return np.sum(np.log(a+eps)*(a+eps), axis=0)


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    parser.add_argument('--mask', type=str, nargs='+')
    args = parser.parse_args()

    fs = 16
    plt.figure()
    for mask in args.mask
        files = glob.glob(mask)
        allS = []
        for file in files:
            m = load_model(file)
            S = entropy(m.Pi)
            allS.extend(S)

        plt.plot(sorted(allS), np.linspace(0,1,len(allS)))
        print(f"{mask}, average entropy: {np.mean(allS)}")
    plt.ylabel('fraction with less entropy', fontsize=fs)
    plt.xlabel('fraction of sites', fontsize=fs)
    plt.savefig('figures/entropy_stats.pdf')
