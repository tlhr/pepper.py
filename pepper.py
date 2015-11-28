from itertools import product
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt


masses = {'Boc':  101.126,
          'tBu':  73.116,
          'Tfa':  97.017,
          'Pbf':  253.338,
          'Trt':  243.331}


scs = {'Boc': ['K', 'H', 'W'],
       'tBu': ['Y', 'S', 'T', 'D', 'E'],
       'NH2': ['K', 'H', 'W', 'Q', 'N'],
       'Pbf': ['R'],
       'Trt': ['K', 'H', 'W', 'Q', 'N']}


def pepattr(seq):
    permut = {k: 0 for k in scs.keys()}
    # Permute every amino acid
    for aa in seq:
        for sc in scs.keys():
            if aa in scs[sc]:
                permut[sc] += 1
    # Output ranges instead of amounts
    return {k: range(0, 1 + v) for k, v in permut.items()}


def pepcut(seq, N=3):
    seqlist = [seq]
    for n in range(N):
        # Cuts n amino acids from the sequence
        seqlist.extend([seq[:i] + seq[i+n+1:] for i, aa in enumerate(seq)])
    return seqlist


def permute(seq, target_mass, eps=0.9):
    mol = Chem.rdmolfiles.MolFromSequence(seq)
    orig_mass = ExactMolWt(mol)
    for p in pepcut(seq, 3):
        print(p)
        # Create all possible permutations from side chains and sequences
        for pp in product(*pepattr(p).values()):
            test_mass = orig_mass
            for i, ppp in enumerate(pp):
                test_mass += ppp * list(masses.values())[i]
            if abs(test_mass - target_mass) <= eps:
                print("Found possible target: M = {}".format(test_mass))
                print("Permutation: {}".format(dict(zip(masses.keys(), pp))))
                print("Sequence: {}".format(p))


def main():
    seq = 'KLKARIRRLFGGY'
    permute(seq, 4029., 1.0)


if __name__ == '__main__':
    main()
