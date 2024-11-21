import numpy as np
from codeCandidate import sourceInterval, codeCandidate, codeCandidateList


def updateSrcInterval(srcInterval, srcProbability, srcSymbols):
    """
    Parameters
    ----------
    srcInterval:   sourceInterval
                    current interval of current bit sequence
    srcProbability:    list or ndarray
                        contains the probability of 0, 1 both equal to 0.5
    srcSymbols:   int
                    the coming bit (0 or 1)

    Returns
    -------
    src_interval_copy:  sourceInterval
                        the updated interval
    """
    srcIntervalCopy = srcInterval
    new_border = srcIntervalCopy.lowerBound + (srcIntervalCopy.upperBound - srcIntervalCopy.lowerBound) * \
                 srcProbability[0]
    if srcSymbols == 0:
        srcIntervalCopy.upperBound = new_border
    else:
        srcIntervalCopy.lowerBound = new_border
    return srcIntervalCopy


def updateNi(Ni, symbols):
    n_i_copy = Ni
    for each in symbols:
        n_i_copy[each - 1] -= 1
    return n_i_copy


def findIdentifiedCodeCandidateIndex(srcInterval, codeCandiList, ni):
    """

    Parameters
    ----------
    srcInterval
    codeCandiList
    ni

    Returns
    -------

    """
    i = 0
    for each in codeCandiList.codeCandidate:
        if (srcInterval.lowerBound <= each.lowerBound) & (srcInterval.upperBound > each.lowerBound) & (ni[i] != 0):
            return i
        i += 1
    return -1


def findCodeIntervalFromCandidate(codeCandiList, searchList):
    """

    Parameters
    ----------
    codeCandiList
    searchList

    Returns
    -------

    """
    code_interval = sourceInterval()
    for each in codeCandiList.codeCandidate:
        cmp_list = each.symbols
        if cmp_list == searchList:
            code_interval.lowerBound = each.lowerBound
            code_interval.upperBound = each.upperBound
            searchList = []
            break
    return code_interval, searchList


def updateCodeCandidates(codeCandiList, ni):
    """

    Parameters
    ----------
    codeCandiList
    ni

    Returns
    -------

    """
    codeCandiList.codeCandidate = []
    n = 0
    k = codeCandiList.k
    sumPi = 0
    pi = np.zeros(k)
    for each in ni:
        n += each
    for i in range(k):
        codeCandi = codeCandidate()
        pi[i] = float(ni[i] / n)
        codeCandi.lowerBound = sumPi
        sumPi += pi[i]
        if i == k - 1:
            codeCandi.upperBound = 1.0
        else:
            codeCandi.upperBound = sumPi
        codeCandi.probability = pi[i]
        codeCandi.symbols.append(i + 1)
        codeCandiList.codeCandidate.append(codeCandi)
    return codeCandiList


def checkForOutputAndRescale(srcInterval, codeCandiList, ni, indexConstDict):
    """
    decide for the output symbol and rescale src_interval
    Parameters
    ----------
    indexConstDict: dict
                    the reflection between indexes to constellation points
    srcInterval:   sourceInterval
                    current interval
    codeCandiList:    codeCandidateList
                            contains all the amplitude and their probability, upper border, lower border
    ni:    int[]
            the remain num of all the amplitude
    Returns
    -------
    code_symbols_new:   list, maybe ndarray later
                        the output sequence when current bit comes, can be empty list
    src_interval:   sourceInterval
                    the rescaled src_interval
    code_candidate_list:    codeCandidateList
                            the updated code_candidate_list

    """
    k = codeCandiList.k
    codeSymbol = None
    codeSymbolsNew = []
    success = 0
    for each in codeCandiList.codeCandidate:
        codeSymbol = each
        if (srcInterval.lowerBound >= each.lowerBound) & (srcInterval.upperBound <= each.upperBound):
            success = 1
            break

    while success:
        srcInterval.lowerBound = (srcInterval.lowerBound - codeSymbol.lowerBound) / (
                codeSymbol.upperBound - codeSymbol.lowerBound)
        srcInterval.upperBound = (srcInterval.upperBound - codeSymbol.lowerBound) / (
                codeSymbol.upperBound - codeSymbol.lowerBound)
        if srcInterval.upperBound > 1.0:
            srcInterval.upperBound = 1.0
        ni = updateNi(ni, codeSymbol.symbols)
        codeSymbolsNew.extend([indexConstDict.get(codeSymbol.symbols[0])])
        codeCandiList = updateCodeCandidates(codeCandiList, ni)

        success = 0
        for each in codeCandiList.codeCandidate:
            codeSymbol = each
            if (srcInterval.lowerBound >= codeSymbol.lowerBound) & (
                    srcInterval.upperBound <= codeSymbol.upperBound):
                success = 1
                break

    return codeSymbolsNew, srcInterval, codeCandiList


##########################################################

def finaliseCodeSymbols(srcInterval, codeCandiList, ni, indexConstDict):
    """

    Parameters
    ----------
    indexConstDict
    codeCandiList
    srcInterval
    ni

    Returns
    -------

    """
    k = codeCandiList.k
    symbolsNew = []
    index = findIdentifiedCodeCandidateIndex(srcInterval, codeCandiList, ni)
    if index == -1:
        return symbolsNew
    symbol = codeCandiList.codeCandidate[index].symbols
    symbolsNew.extend([indexConstDict.get(symbol[0])])
    ni = updateNi(ni, symbol)
    for i in range(k):
        for j in range(int(ni[i])):
            symbolsNew.append(indexConstDict.get(i + 1))
    return symbolsNew


def encode(srcSymbols, n, ni, indexConstDict):
    """
    Parameters
    ----------
    indexConstDict
    srcSymbols:    ndarray or list
                    source_symbols sequence consists of 0, 1
    n:  int
        the length of output sequence
    ni:    int[]
            the times of each symbols' appearance
    Returns
    -------

    """
    k = len(ni)  
    srcProbability = [0.5, 0.5]  # probability of 0, 1
    srcInterval = sourceInterval()  # initialise the interval

    codeCandiList = codeCandidateList()
    codeCandiList.k = k

    niCopy = ni
    sumPi = 0

    codeSymbols = []  # generated sequence
    pi = np.zeros(k)

    # code_candidate_list
    for i in range(k):
        pi[i] = float(niCopy[i] / n)
        codeCandi = codeCandidate()
        codeCandi.lowerBound = sumPi
        sumPi += pi[i]
        if i == k - 1:
            codeCandi.upperBound = 1.0
        else:
            codeCandi.upperBound = sumPi
        codeCandi.symbols.append(i + 1)
        codeCandiList.codeCandidate.append(codeCandi)

    for symbol in srcSymbols:
        srcInterval = updateSrcInterval(srcInterval, srcProbability, symbol)
        newSym, srcInterval, codeCandiList = checkForOutputAndRescale(srcInterval, codeCandiList,
                                                                              niCopy, indexConstDict)
        if len(newSym) != 0:
            codeSymbols.extend(newSym)

    finSymbols = finaliseCodeSymbols(srcInterval, codeCandiList, ni, indexConstDict)  # padding generated sequence
    codeSymbols.extend(finSymbols)
    return codeSymbols


def decode(codeSymbols, ni, m):
    """

    Parameters
    ----------
    codeSymbols
    ni
    m

    Returns
    -------

    """
    n = len(codeSymbols)
    k = len(ni)
    futureNi = np.zeros(k)

    srcInterval = sourceInterval()
    # code_interval = sourceInterval()

    codeCandiList = codeCandidateList()
    codeCandiList.k = k
    codeCandiList = updateCodeCandidates(codeCandiList, ni)

    srcSymbols = np.zeros(m)
    srcSymbolIndex = 0
    codeSymbolsUnprocessed = []

    codeSymbolIndex = 0

    while codeSymbolIndex != len(codeSymbols) - 1:
        codeSymbolsUnprocessed.append(codeSymbols[codeSymbolIndex])
        codeInterval, codeSymbolsUnprocessed = findCodeIntervalFromCandidate(codeCandiList,
                                                                                codeSymbolsUnprocessed)

        if len(codeSymbolsUnprocessed) != 0:
            continue

        codeIntervalIndex = codeSymbolIndex

        futureN = n
        for i in range(k):
            futureNi[i] = ni[i]
        futureNi[codeSymbols[codeSymbolIndex] - 1] -= 1
        futureN -= 1

        codeIntervalIndex += 1

        performedScaling = 0

        while performedScaling == 0:
            newBorder = srcInterval.lowerBound + (srcInterval.upperBound - srcInterval.lowerBound) * 0.5
            while (codeInterval.lowerBound >= newBorder) | (codeInterval.upperBound < newBorder):
                if codeInterval.lowerBound >= newBorder:
                    srcSymbols[srcSymbolIndex] = 1
                    srcSymbolIndex += 1
                    srcInterval.lowerBound = newBorder
                elif codeInterval.upperBound < newBorder:
                    srcSymbols[srcSymbolIndex] = 0
                    srcSymbolIndex += 1
                    srcInterval.upperBound = newBorder

                if srcSymbolIndex >= m:
                    return srcSymbols  

                newBorder = srcInterval.lowerBound + (srcInterval.upperBound - srcInterval.lowerBound) * 0.5

                checkSrcInterval = srcInterval

                checkCode, checkSrcInterval, codeCandiList = checkForOutputAndRescale(checkSrcInterval,
                                                                                               codeCandiList,
                                                                                               ni)

                if (checkSrcInterval.lowerBound != srcInterval.lowerBound) | (
                        checkSrcInterval.upperBound != srcInterval.upperBound):
                    # for i in range(len(check_code)):
                    codeSymbolIndex += 1 * len(checkCode)
                    # code_symbol = code_symbols[code_symbol_index]
                    srcInterval = checkSrcInterval
                    performedScaling = 1
                    break

            if codeIntervalIndex == len(codeSymbols) - 1:
                codeInterval.upperBound = codeInterval.lowerBound + (
                        codeInterval.upperBound - codeInterval.lowerBound) * 0.01
            else:
                codeIntervalSymbol = codeSymbols[codeIntervalIndex]
                buffer = sourceInterval()
                n = 0
                for i in range(k):
                    n += futureNi[i]
                sumNi = 0
                for i in range(codeIntervalSymbol - 1):
                    sumNi += futureNi[i]
                lbound = float(sumNi / n)
                sumNi += futureNi[codeIntervalSymbol - 1]
                ubound = float(sumNi / n)
                # if ubound < lbound:
                #     ubound = 1
                buffer.lowerBound = codeInterval.lowerBound + (
                        codeInterval.upperBound - codeInterval.lowerBound) * lbound
                buffer.upperBound = codeInterval.lowerBound + (
                        codeInterval.upperBound - codeInterval.lowerBound) * ubound

                codeInterval = buffer
                futureNi[codeIntervalSymbol - 1] -= 1
                futureN -= 1
                codeIntervalIndex += 1


if __name__ == "__main__":
    p = [0.1, 0.1, 0.3, 0.5]
    num_info_bits = 12
    n_i = [1, 2, 2, 5]
    p_iter_dis_quant = [0.1, 0.2, 0.2, 0.5]
    src_symbols = np.random.randint(0, 2, 12)
    print("src_symbols: " + str(src_symbols))
    # print(len(src_symbols))
    code_symbols = encode(src_symbols, 10, n_i)
    print("code_symbols: " + str(code_symbols))
    # code_symbols = [4, 2, 4, 4, 3, 3, 4, 2, 4, 1]
    # decoded_src_symbols = decode(code_symbols, n_i, len(src_symbols))
    # print("code_symbols: " + str(decoded_src_symbols))
