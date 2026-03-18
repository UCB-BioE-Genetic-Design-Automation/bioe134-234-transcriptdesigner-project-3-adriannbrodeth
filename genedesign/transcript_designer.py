import random
import csv
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker

class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence using weighted random 
    selection from a codon usage table. It iteratively validates the design against 
    hairpin, promoter, and forbidden sequence constraints.
    """

    def __init__(self):
        self.rbs_chooser = None
        self.promoter_checker = PromoterChecker()
        self.forbidden_checker = ForbiddenSequenceChecker()
        self.codon_checker = CodonChecker()
        self.amino_acid_to_codon_map = {}

    def initiate(self) -> None:
        """
        Initializes the checkers and parses the codon usage data.
        """
        self.rbs_chooser = RBSChooser()
        self.rbs_chooser.initiate()
        
        self.promoter_checker.initiate()
        self.forbidden_checker.initiate()
        self.codon_checker.initiate()

        # Load codon usage data for the random translation weights
        # Note: Ensure this path matches your directory structure.
        codon_usage_file = 'genedesign/data/codon_usage.txt'
        
        try:
            with open(codon_usage_file, 'r') as f:
                for line in f:
                    parts = line.split()
                    # Skip header lines or source tags (e.g., )
                    if len(parts) >= 3:
                        codon = parts[0].strip()
                        aa = parts[1].strip()
                        try:
                            freq = float(parts[2].strip())
                            if aa not in self.amino_acid_to_codon_map:
                                self.amino_acid_to_codon_map[aa] = []
                            self.amino_acid_to_codon_map[aa].append((codon, freq))
                        except ValueError:
                            # Skip lines where frequency is not a valid float
                            continue
        except FileNotFoundError:
            print(f"Error: Could not find {codon_usage_file}. Check your file path.")

    def _generate_candidate_codons(self, peptide: str) -> list:
        """
        Generates a candidate list of codons based on usage frequencies.
        """
        codons = []
        for aa in peptide:
            options = self.amino_acid_to_codon_map.get(aa)
            if not options:
                continue
            
            # Perform weighted random selection based on natural frequencies
            choices = [opt[0] for opt in options]
            weights = [opt[1] for opt in options]
            codons.append(random.choices(choices, weights=weights)[0])
        
        codons.append("TAA")  # Default stop codon
        return codons

    def run(self, peptide: str, ignores: set, max_attempts: int = 500) -> Transcript:
        """
        Iteratively generates and validates transcripts until all criteria are met.
        """
        for attempt in range(max_attempts):
            # 1. Generate candidate sequence
            codons = self._generate_candidate_codons(peptide)
            cds_seq = "".join(codons)
            
            # 2. Check Codon Quality (CAI, Diversity, Rare Codons)
            # This uses the correctly named 'codon_checker' attribute
            c_pass, div, rare, cai = self.codon_checker.run(codons)
            if not c_pass:
                continue
                
            # 3. Check Forbidden Sequences (Restriction sites, poly-A, etc.)
            f_pass, site = self.forbidden_checker.run(cds_seq)
            if not f_pass:
                continue
                
            # 4. Check Internal Promoters
            p_pass, promo = self.promoter_checker.run(cds_seq) 
            if not p_pass:
                continue
                
            # 5. Check for Bad Hairpins
            h_pass, hairpin_str = hairpin_checker(cds_seq) 
            if not h_pass:
                continue
            
            # 6. Success: Select RBS and Return
            selected_rbs = self.rbs_chooser.run(cds_seq, ignores)
            print(f"Design succeeded on attempt {attempt + 1} with CAI: {cai:.3f}")
            return Transcript(selected_rbs, peptide, codons)

        raise Exception(f"Failed to find a valid sequence after {max_attempts} attempts.")

if __name__ == "__main__":
    designer = TranscriptDesigner()
    designer.initiate()
    
    # Example usage
    test_peptide = "MYPFIRTARMTV"
    try:
        result = designer.run(test_peptide, set())
        print(f"Final DNA: {''.join(result.codons)}")
    except Exception as e:
        print(e)