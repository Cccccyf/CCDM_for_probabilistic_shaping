import numpy as np


def iterDiscreteQuant(pOpt, nSyms):
    """
    Parameters
    ----------
    pOpt:   ndarray
            the probability distribution given
    nSyms:  int
            the symbol num
    Returns:    n_i: ndarray
                    the number of every symbol that conform the optimal distribution
                p_iter_dis_quant: ndarray
                    the discrete probability after discrete quantification
    -------

    """
    m = len(pOpt)  
    ni = np.zeros(m)
    t = np.log(np.reciprocal(pOpt))
    pIterDisQuant = np.copy(t)

    for i in range(nSyms):
        #print(p_iter_dis_quant)
        index = np.argmin(pIterDisQuant)
        cj = ni[index] + 1
        ni[index] = cj
        pIterDisQuant[index] = (cj + 1) * np.log(cj + 1) - cj * np.log(cj) + t[index]
    pIterDisQuant = ni / nSyms
    # print(n_i)
    # print(p_iter_dis_quant)
    return ni, pIterDisQuant


def nChooseKLog2(n, k):
    """

    Parameters
    ----------
    n: n
    k

    Returns
    -------

    """
    nCK = 0
    if k > n - k:
        k = n - k
    for i in range(1, k + 1):
        nCK += np.log((n - (k - i)) / i) / np.log(2)
    return nCK


def nChooseKsLog2(n, ks):
    """
    get the num of combination
    Parameters
    ----------
    n:
    ks: the list of the pick num

    Returns
    -------
    the num of combination
    """
    nn = n
    nCKs = 0
    ks = np.sort(ks)
    for each in ks:
        nCKs += nChooseKLog2(int(nn), int(each))
        nn -= each
    return nCKs


def initialize(pOpt, nSyms):
    ni, pIterDisQuant = iterDiscreteQuant(pOpt, nSyms)
    numInfoBits = np.floor(nChooseKsLog2(nSyms, ni))
    return ni, pIterDisQuant, int(numInfoBits)


if __name__ == "__main__":
    pOpt = [0.1, 0.2, 0.2, 0.5]
    n = 10
    print(initialize(pOpt, n))
