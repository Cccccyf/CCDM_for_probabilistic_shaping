class sourceInterval:
    def __init__(self):
        self.upperBound = 1.0
        self.lowerBound = 0.0

    def updateInterval(self, upperBound, lowerBound):
        self.upperBound = upperBound
        self.lowerBound = lowerBound


class codeCandidate:
    def __init__(self):
        self.upperBound = None
        self.lowerBound = None
        self.probability = None  # float
        self.symbols = []


class codeCandidateList:
    def __init__(self):
        self.k = None
        self.codeCandidate = []
