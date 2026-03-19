class RNaseEChecker:
    def __init__(self):
        self.forbidden_sites = []

    def initiate(self):
        # Consensus: RAUUW -> DNA: (A/G)ATT(A/T)
        self.forbidden_sites = ["AATTA", "AATTT", "GATTA", "GATTT"]

    def run(self, dna_seq: str):
        for site in self.forbidden_sites:
            if site in dna_seq:
                return False, site
        return True, None