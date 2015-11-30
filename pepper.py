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
    """
    Create possible side chain mutations from sequence.

    Parameters
    ----------
    seq : string
        One letter code sequence string

    Returns
    -------
    permut : dict
        Dictionary containing possible protecting groups
        and the ranges of possible additions.

    """
    permut = {k: 0 for k in scs.keys()}
    # Permute every amino acid
    for aa in seq:
        for sc in scs.keys():
            if aa in scs[sc]:
                permut[sc] += 1
    # Output ranges instead of amounts
    return {k: range(0, 1 + v) for k, v in permut.items()}


def pepcut(seqlist, n=3):
    """
    Cuts amino acids from the sequence string.

    Parameters
    ----------
    seq : string
        One letter code sequence string
    N : int (default=3)
        Number of amino acids to cut from the sequence

    Returns
    -------
    seqlist : list
        List of all possible permutations

    """
    if n == 0:
        return [seqlist]
    else:
        plist = []
        for seq in seqlist:
            for i, _ in enumerate(seq):
                plist.append(seq[:i] + seq[i+1:])
        return pepcut(plist, n-1)


def permute(seq, target_mass, eps=0.9, N=3):
    """
    Permutes the sequence and compares the permutations to the target mass
    with a given accuracy epsilon.

    Parameters
    ----------
    seq : string
        One letter code sequence string
    target_mass : float
        The target mass to compare to the permutation
    eps : float (default=0.9)
        Maximum deviation of the test mass from the target mass.
        Should be chosen according to the accuracy of the mass
        spectrometer instrument.
    N : int (default=3)
        Number of amino acids to cut from the sequence

    """
    mol = Chem.rdmolfiles.MolFromSequence(seq)
    orig_mass = ExactMolWt(mol)
    for p in pepcut([seq], N):
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
