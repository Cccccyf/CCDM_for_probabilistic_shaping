import numpy as np
from utils import grayMapping
from CCDM import encode
from CCDMUtils import initialize



def getMQamConstellations(M):
    """

    Parameters
    ----------
    M

    Returns
    -------

    """
    return grayMapping(M, "qam")


class probabilisticShaping:
    def __init__(self, M, shapingFactor, symbolNum):
        self.M = M
        self.shapingFactor = shapingFactor
        self.symbolNum = symbolNum
        self.const = getMQamConstellations(self.M)
        self.probabilities = self.genProbabilities()
        # enlarge the constellation of probabilistic shaping to make the shaped signals have equal power as signals following uniform distrition 
        self.const = np.sqrt(np.mean(np.abs(self.const) ** 2) / np.sum(np.abs(self.const) ** 2 * self.probabilities)) * self.const
        self.indexConstDict = {}      
        self.constIndexDict = {}      #
        self.initDict()

    def initDict(self):
        for i in range(len(self.probabilities)):
            self.indexConstDict[i + 1] = self.const[i]
            self.constIndexDict[self.const[i]] = i + 1

    def genProbabilities(self):
        """
        the probability of M-QAM probabilistic shaping 's  each constellation points conforming Maxwell-Boltzmann distribution
        Parameters
        ----------

        Returns
        -------

        """
        constAbs = np.abs(self.const)
        probabilities = np.exp(-self.shapingFactor * (constAbs ** 2)) / (
            np.sum(np.exp(-self.shapingFactor * (constAbs ** 2))))
        return probabilities.flatten()

    def genSymbols(self):
        ni, pIterDisQuant, numInfoBits = initialize(self.probabilities, self.symbolNum)
        print(numInfoBits)
        bitsSeq = np.random.randint(0, 2, numInfoBits)
        symbolsSeq = encode(bitsSeq, self.symbolNum, ni, self.indexConstDict)
        return symbolsSeq


if __name__ == "__main__":
    # print(getMQamConstellations(16))
    ps = probabilisticShaping(M=16, shapingFactor=0.1, symbolNum=1000)
    print(np.mean(np.abs(getMQamConstellations(16)) ** 2))
    # print(ps.index_constellation_dict)
    #print(ps.index_constellation_dict.get(1))
    print(np.sum(np.abs(ps.const) ** 2 * ps.probabilities))
    # print(np.abs(get_mqam_constellations(16).flatten()))
