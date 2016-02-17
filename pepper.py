from itertools import product
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt


orig_masses = {'Boc':  101.126,
               'tBu':  73.116,
               'Tfa':  97.017,
               'Pbf':  253.338,
               'Trt':  243.331,
               'red':  -28.0140 + 2 * 1.008,
               'Fmoc': 223.254}

misc_masses = {'H2O':  18.0153,
               'H':    1.008}

aamw = {'I': 131.1736,
        'L': 131.1736,
        'K': 146.1882,
        'M': 149.2124,
        'F': 165.1900,
        'T': 119.1197,
        'W': 204.2262,
        'V': 117.1469,
        'R': 174.2017,
        'H': 155.1552,
        'A': 89.0935,
        'N': 132.1184,
        'D': 133.1032,
        'C': 121.1590,
        'E': 147.1299,
        'Q': 146.1451,
        'G': 75.0669,
        'P': 115.1310,
        'S': 105.0930,
        'Y': 181.1894,
        'a': 113.1160,
        'b': 206.2050}

masses = {k: v - misc_masses['H'] for k, v in orig_masses.items()}

scs = {'Boc': ['K', 'H', 'W'],
       'tBu': ['Y', 'S', 'T', 'D', 'E'],
       'NH2': ['K', 'H', 'W', 'Q', 'N'],
       'Pbf': ['R'],
       'Trt': ['K', 'H', 'W', 'Q', 'N'],
       'red': ['b']}


def pepweight(seq):
    """
    Calculate peptide mass based on MW table.

    Parameters
    ----------
    seq : string
        One letter code sequence string

    Returns
    -------
    peptide_weight : float
        Exact mol. weight of peptide

    """
    return (sum(aamw[aa] for aa in seq) -
            (len(seq) - 1) * misc_masses['H2O'])


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
    rperm = {k: range(0, 1 + v) for k, v in permut.items()}
    rperm['Fmoc'] = [1]
    return rperm


def pepcut(seqlist, n=3):
    """
    Cuts amino acids from the sequence string.

    Parameters
    ----------
    seqlist : list
        List of one letter code sequence strings
    n : int (default=3)
        Number of amino acids to cut from the sequence

    Returns
    -------
    seqlist : list
        List of all possible permutations

    """
    if n == 0:
        return seqlist
    else:
        plist = []
        for seq in seqlist:
            for i, _ in enumerate(seq):
                plist.append(seq[:i] + seq[i+1:])
        return pepcut(plist, n-1) + seqlist


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
    prev_mass = orig_mass = 0
    for p in pepcut([seq], N):
        prev_mass = orig_mass
        # mol = Chem.rdmolfiles.MolFromSequence(p)
        # orig_mass = ExactMolWt(mol)
        orig_mass = pepweight(seq)
        if prev_mass == orig_mass:
            continue
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
    seq = 'YGGbLRRIRAKaK'
    # mol = Chem.rdmolfiles.MolFromSequence(seq)
    # rdw = ExactMolWt(mol)
    # w = pepweight(seq)
    # print(rdw, w, abs(rdw - w))
    permute(seq, 1762., 9.9, 3)
    # AAKLKDRRRKKD +Boc, +Pbf


if __name__ == '__main__':
    main()
