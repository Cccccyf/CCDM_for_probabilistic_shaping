import numpy as np
"""
Edson Porto da Silva, Adolfo Fernandes Herbster. "OptiCommPy: Open-source Simulation of Fiber Optic Communications with Python", Journal of Open Source Software, 9(98), 6600, (2024) 
https://doi.org/10.21105/joss.06600
"""

def pamConst(M):
    """
    Generate a Pulse Amplitude Modulation (PAM) constellation.

    Parameters
    ----------
    M : int
        Number of symbols in the constellation. It must be an integer.

    Returns
    -------
    np.array
        1D PAM constellation.
    
    References
    ----------
    [1] Proakis, J. G., & Salehi, M. Digital Communications (5th Edition). McGraw-Hill Education, 2008.
    """
    L = int(M - 1)
    return np.arange(-L, L + 1, 2)

def pskConst(M):
    """
    Generate a Phase Shift Keying (PSK) constellation.

    Parameters
    ----------
    M : int
        Number of symbols in the constellation. It must be a power of 2 positive integer.

    Returns
    -------
    np.array
        Complex M-PSK constellation.

    References
    ----------
    [1] Proakis, J. G., & Salehi, M. Digital Communications (5th Edition). McGraw-Hill Education, 2008.
    """
    # generate complex M-PSK constellation
    pskPhases = np.arange(0, 2 * np.pi, 2 * np.pi / M)
    return np.exp(1j * pskPhases)


def apskConst(M, m1=None, phaseOffset=None):
    """
    Generate an Amplitude-Phase Shift Keying (APSK) constellation.

    Parameters
    ----------
    M : int
        Constellation order.
    m1 : int
        Number of bits used to index the radii of the constellation.

    Returns
    -------
    const : np.array
        APSK constellation

    References
    ----------
    [1] Z. Liu, et al "APSK Constellation with Gray Mapping," IEEE Communications Letters, vol. 15, no. 12, pp. 1271-1273, 2011.

    [2] Proakis, J. G., & Salehi, M. Digital Communications (5th Edition). McGraw-Hill Education, 2008.
    """
    if m1 is None:
        if M == 16:
            m1 = 1
        elif M == 32:
            m1 = 2
        elif M == 64:
            m1 = 2
        elif M == 128:
            m1 = 3
        elif M == 256:
            m1 = 3
        elif M == 512:
            m1 = 4
        elif M == 1024:
            m1 = 4

    nRings = int(2**m1)  # bits that index the rings
    m2 = int(np.log2(M) - m1)  # bits that index the symbols per ring

    symbolsPerRing = int(2**m2)

    const = np.zeros((M,), dtype=np.complex64)

    if phaseOffset is None:
        phaseOffset = np.pi / symbolsPerRing

    for idx in range(nRings):
        radius = np.sqrt(-np.log(1 - ((idx + 1) - 0.5) * symbolsPerRing / M))

        if (idx + 1) % 2 == 1:
            const[idx * symbolsPerRing : (idx + 1) * symbolsPerRing] = radius * np.flip(
                pskConst(symbolsPerRing)
            )
        else:
            const[
                idx * symbolsPerRing : (idx + 1) * symbolsPerRing
            ] = radius * pskConst(symbolsPerRing)

    return const * np.exp(1j * phaseOffset)

def qamConst(M):
    """
    Generate a Quadrature Amplitude Modulation (QAM) constellation.

    Parameters
    ----------
    M : int
        Number of symbols in the constellation. It must be a perfect square.

    Returns
    -------
    const : np.array
        Complex square M-QAM constellation.

    References
    ----------
    [1] Proakis, J. G., & Salehi, M. Digital Communications (5th Edition). McGraw-Hill Education, 2008.
    """
    L = int(np.sqrt(M) - 1)

    # generate 1D PAM constellation
    PAM = np.arange(-L, L + 1, 2)
    PAM = np.array([PAM])

    # generate complex square M-QAM constellation
    const = np.tile(PAM, (L + 1, 1))
    const = const + 1j * np.flipud(const.T)

    for ind in np.arange(1, L + 1, 2):
        const[ind] = np.flip(const[ind], 0)

    return const

def grayCode(n):
    """
    Gray code generator.

    Parameters
    ----------
    n : int
        length of the codeword in bits.

    Returns
    -------
    code : list
           list of binary strings of the gray code.

    """
    code = []

    for i in range(1 << n):
        # Generating the decimal
        # values of gray code then using
        # bitset to convert them to binary form
        val = i ^ (i >> 1)

        # Converting to binary string
        s = bin(val)[2::]
        code.append(s.zfill(n))
    return code

def grayMapping(M, constType):
    """
    Gray Mapping for digital modulations.

    Parameters
    ----------
    M : int
        modulation order
    constType : 'qam', 'psk', 'pam' or 'ook'.
        type of constellation.

    Returns
    -------
    const : np.array
        constellation symbols (sorted according their corresponding
        Gray bit sequence as integer decimal).
    
    References
    ----------
    [1] Proakis, J. G., & Salehi, M. Digital Communications (5th Edition). McGraw-Hill Education, 2008.
    """
    if M != 2 and constType == "ook":
        M = 2

    bitsSymb = int(np.log2(M))

    code = grayCode(bitsSymb)
    if constType == "ook":
        const = np.arange(0, 2)
    elif constType == "pam":
        const = pamConst(M)
    elif constType == "qam":
        const = qamConst(M)
    elif constType == "psk":
        const = pskConst(M)
    elif constType == "apsk":
        const = apskConst(M)

    const = const.reshape(M, 1)
    const_ = np.zeros((M, 2), dtype=complex)           

    for ind in range(M):
        const_[ind, 0] = const[ind, 0]  # complex constellation symbol
        const_[ind, 1] = int(code[ind], 2)  # mapped bit sequence (as integer decimal)

    # sort complex symbols column according to their mapped bit sequence (as integer decimal)
    const = const_[const_[:, 1].real.argsort()]
    const = const[:, 0]

    if constType in ["pam", "ook"]:
        const = const.real
    return const
