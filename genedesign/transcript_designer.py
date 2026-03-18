import random
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker

class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and selects an RBS.
    Validated against multiple biological constraints with robust error handling.
    """

    MAX_ATTEMPTS = 10000

    def __init__(self):
        self.codon_table = {}
        self.rbsChooser = RBSChooser()
        self.codonChecker = CodonChecker()
        self.forbiddenChecker = ForbiddenSequenceChecker()
        self.promoterChecker = PromoterChecker()

    def initiate(self) -> None:
        self.rbsChooser.initiate()
        self.codonChecker.initiate()
        self.forbiddenChecker.initiate()
        self.promoterChecker.initiate()
        self.codon_table = self._default_codon_table()

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates peptide to DNA, ensuring it passes all checker constraints.
        """
        for attempt in range(1, self.MAX_ATTEMPTS + 1):
            codons = self._sample_codons(peptide)
            cds = ''.join(codons)

            # FIX: Pass both the string (for structural checks) and list (for codon checks)
            failures = self._check_sequence(cds, codons)
            
            if not failures:
                selected_rbs = self.rbsChooser.run(cds, ignores)
                # Ensure the RBS chooser actually returned a sequence
                if selected_rbs:
                    return Transcript(selected_rbs, peptide, codons)

        raise RuntimeError(f"Failed to find valid sequence after {self.MAX_ATTEMPTS} attempts.")

    def _sample_codons(self, peptide: str) -> list[str]:
        codons = []
        for aa in peptide:
            options = self.codon_table.get(aa)
            if not options:
                raise ValueError(f"Unknown amino acid: '{aa}'")
            
            c_list = [opt[0] for opt in options]
            w_list = [opt[1] for opt in options]
            codons.append(random.choices(c_list, weights=w_list, k=1)[0])

        stop_options = self.codon_table.get('*')
        s_list = [opt[0] for opt in stop_options]
        sw_list = [opt[1] for opt in stop_options]
        codons.append(random.choices(s_list, weights=sw_list, k=1)[0])
        return codons

    def _check_sequence(self, cds: str, codons: list[str]) -> list[str]:
        """
        Runs all checkers, ensuring the correct data format is sent to each.
        """
        failures = []
        
        # 1. Structural Checkers (Expect DNA String)
        forbidden_ok, forbidden_site = self.forbiddenChecker.run(cds)
        if not forbidden_ok:
            failures.append(f"Forbidden: {forbidden_site}")

        promoter_ok, promoter_seq = self.promoterChecker.run(cds)
        if not promoter_ok:
            failures.append(f"Promoter: {promoter_seq}")

        hairpin_ok, hairpin_seq = hairpin_checker(cds)
        if not hairpin_ok:
            failures.append(f"Hairpin: {hairpin_seq}")

        # 2. Codon Usage Checker (Expects List of Codons)
        codon_res = self.codonChecker.run(codons)
        if not codon_res[0]:
            # codon_res[1] is diversity, [2] is rare count, [3] is CAI
            failures.append(f"Codon: Diversity={codon_res[1]:.2f}, RareCount={codon_res[2]}, CAI={codon_res[3]:.2f}")
        
        return failures

    @staticmethod
    def _default_codon_table() -> dict:
        return {
            'A': [('GCT', 0.18), ('GCC', 0.26), ('GCA', 0.23), ('GCG', 0.33)],
            'R': [('CGT', 0.36), ('CGC', 0.40), ('CGA', 0.07), ('CGG', 0.11), ('AGA', 0.04), ('AGG', 0.02)],
            'N': [('AAT', 0.45), ('AAC', 0.55)], 'D': [('GAT', 0.63), ('GAC', 0.37)],
            'C': [('TGT', 0.45), ('TGC', 0.55)], 'Q': [('CAA', 0.34), ('CAG', 0.66)],
            'E': [('GAA', 0.68), ('GAG', 0.32)], 'G': [('GGT', 0.35), ('GGC', 0.40), ('GGA', 0.11), ('GGG', 0.14)],
            'H': [('CAT', 0.57), ('CAC', 0.43)], 'I': [('ATT', 0.49), ('ATC', 0.39), ('ATA', 0.11)],
            'L': [('TTA', 0.14), ('TTG', 0.13), ('CTT', 0.12), ('CTC', 0.10), ('CTA', 0.04), ('CTG', 0.47)],
            'K': [('AAA', 0.74), ('AAG', 0.26)], 'M': [('ATG', 1.00)], 'F': [('TTT', 0.58), ('TTC', 0.42)],
            'P': [('CCT', 0.18), ('CCC', 0.13), ('CCA', 0.20), ('CCG', 0.49)],
            'S': [('TCT', 0.17), ('TCC', 0.15), ('TCA', 0.14), ('TCG', 0.14), ('AGT', 0.16), ('AGC', 0.25)],
            'T': [('ACT', 0.19), ('ACC', 0.40), ('ACA', 0.17), ('ACG', 0.25)],
            'W': [('TGG', 1.00)], 'Y': [('TAT', 0.59), ('TAC', 0.41)],
            'V': [('GTT', 0.28), ('GTC', 0.20), ('GTA', 0.17), ('GTG', 0.35)],
            '*': [('TAA', 0.61), ('TAG', 0.09), ('TGA', 0.30)],
        }