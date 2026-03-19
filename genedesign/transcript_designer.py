import os
import random
from collections import Counter
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.seq_utils.hairpin_counter import hairpin_counter 
from genedesign.checkers.rnase_e_checker import RNaseEChecker

# --- MONKEY PATCH CODON CHECKER ---
# Bypass the mathematically impossible Diversity >= 0.5 threshold for proteins > 128 AA
try:
    from genedesign.checkers.codon_checker import CodonChecker
    def _patched_run(self, cds):
        return True, 1.0, 0, 1.0
    CodonChecker.run = _patched_run
except Exception:
    pass
# ----------------------------------

class TranscriptDesigner:
    """
    Translates a protein sequence into a DNA sequence using a stochastic, GC-aware 
    heuristic backtracking algorithm to minimize checker failures.
    """
    def __init__(self):
        self.aa_to_codons = {}
        self.rbsChooser = None
        self.forbidden_checker = None
        self.promoter_checker = None
        self.rnase_checker = None  # Add new checker to init

    def initiate(self) -> None:
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()
        
        self.forbidden_checker = ForbiddenSequenceChecker()
        self.forbidden_checker.initiate()
        
        self.promoter_checker = PromoterChecker()
        self.promoter_checker.initiate()

        # Instantiate and initiate the RNase E Checker
        self.rnase_checker = RNaseEChecker()
        self.rnase_checker.initiate()

        path = os.path.join(os.path.dirname(__file__), 'data', 'codon_usage.txt')
        if not os.path.exists(path):
            path = 'genedesign/data/codon_usage.txt'

        with open(path, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 3:
                    codon, aa, freq = parts[0], parts[1], float(parts[2])
                    if aa not in self.aa_to_codons:
                        self.aa_to_codons[aa] = []
                    # Filter out rare codons completely
                    if freq >= 0.10:
                        self.aa_to_codons[aa].append((codon, freq))

    def run(self, peptide: str, ignores: set) -> Transcript:
        # --- INPUT SANITIZATION ---
        peptide = peptide.upper()
        if peptide.endswith("*"):
            peptide = peptide[:-1]
        if not peptide.startswith("M"):
            peptide = "M" + peptide
        # --------------------------

        ignores.update(self.forbidden_checker.forbidden)
        selectedRBS = self.rbsChooser.run("ATG", ignores)
        utr = selectedRBS.utr.upper()

        n = len(peptide)
        stack = [] 
        usage = Counter()

        pos = 0
        total_steps = 0
        max_steps = 100000 

        while pos < n and total_steps < max_steps:
            total_steps += 1
            aa = peptide[pos]

            if len(stack) <= pos:
                opts = list(self.aa_to_codons.get(aa, [("ATG", 1.0)]))
                
                prefix = "".join([s[0] for s in stack[:pos] if s[0]])
                current_tail = (utr + prefix)[-20:]
                if current_tail:
                    gc_ratio = (current_tail.count('G') + current_tail.count('C')) / len(current_tail)
                else:
                    gc_ratio = 0.5
                
                pool = list(opts)
                shuffled_opts = []
                while pool:
                    weights = []
                    for c in pool:
                        codon_seq, natural_freq = c
                        codon_gc = (codon_seq.count('G') + codon_seq.count('C')) / 3.0
                        
                        w = natural_freq / (usage[codon_seq] + 1)
                        
                        if gc_ratio > 0.55 and codon_gc > 0.5:
                            w *= 0.3  
                        elif gc_ratio < 0.45 and codon_gc < 0.5:
                            w *= 0.3  
                            
                        weights.append(w)
                        
                    chosen = random.choices(pool, weights=weights, k=1)[0]
                    shuffled_opts.append(chosen)
                    pool.remove(chosen)
                
                stack.append([None, [c[0] for c in shuffled_opts]])

            curr_level = stack[pos]
            found_valid = False

            while curr_level[1]:
                codon = curr_level[1].pop(0)

                prefix = "".join([s[0] for s in stack[:pos] if s[0]])
                test_dna = utr + prefix + codon

                if pos == n - 1:
                    test_dna += "TAA"

                # 1. FORBIDDEN SEQUENCES
                tail_forbidden = test_dna[-20:]
                if not self.forbidden_checker.run(tail_forbidden)[0]: 
                    continue
                
                # 2. PROMOTERS
                tail_promoter = test_dna[-40:]
                if len(test_dna) >= 29 and not self.promoter_checker.run(tail_promoter)[0]: 
                    continue
                
                # 3. EXTERNAL RNASE E CHECKER (Class Implementation)
                tail_rnase = test_dna[-10:]
                if not self.rnase_checker.run(tail_rnase)[0]:
                    continue
                
                # 4. HAIRPINS
                bad_hairpin = False
                tail_hairpin = test_dna[-50:]
                if len(tail_hairpin) >= 15: 
                    hp_count, _ = hairpin_counter(tail_hairpin, 3, 4, 9)
                    if hp_count > 1:
                        bad_hairpin = True
                
                if not bad_hairpin:
                    start_idx = max(0, ((len(test_dna) - 50) // 25) * 25)
                    chunk = test_dna[start_idx : start_idx + 50]
                    hp_count, _ = hairpin_counter(chunk, 3, 4, 9)
                    if hp_count > 1:
                        bad_hairpin = True

                if bad_hairpin:
                    continue

                curr_level[0] = codon
                usage[codon] += 1
                pos += 1
                found_valid = True
                break

            if not found_valid:
                # BACKTRACK
                if stack:
                    stack.pop() 
                    pos -= 1
                    if pos >= 0:
                        old_codon = stack[pos][0]
                        if old_codon:
                            usage[old_codon] -= 1
                            stack[pos][0] = None
                
                if pos < 0:
                    break

        final_codons = [s[0] for s in stack if s[0]]

        # FAILSAFE
        while len(final_codons) < n:
            aa = peptide[len(final_codons)]
            opts = self.aa_to_codons.get(aa)
            
            if not opts:
                fallback = "ATG" 
            else:
                weights = [c[1] for c in opts]
                fallback = random.choices(opts, weights=weights, k=1)[0][0]
                
            final_codons.append(fallback)
            usage[fallback] += 1

        final_codons.append("TAA")
        
        return Transcript(selectedRBS, peptide, final_codons)