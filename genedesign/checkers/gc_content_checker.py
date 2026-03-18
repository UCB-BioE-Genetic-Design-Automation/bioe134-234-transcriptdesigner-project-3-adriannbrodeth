class GCContentChecker:
    """
    Checks if the GC content of a DNA sequence is within a specified range.
    High GC content can lead to secondary structures, while low GC can decrease 
    mRNA stability.
    """

    def __init__(self):
        self.min_gc = 0.3  # 30%
        self.max_gc = 0.7  # 70%

    def initiate(self, min_gc: float = 0.3, max_gc: float = 0.7) -> None:
        """
        Sets the GC content thresholds.
        """
        self.min_gc = min_gc
        self.max_gc = max_gc

    def run(self, dna: str) -> tuple[bool, str]:
        """
        Calculates GC content and validates against thresholds.
        
        Returns:
            (True, "GC content is [X]%") if within range.
            (False, "GC content [X]% is out of range") if outside range.
        """
        if not dna:
            return False, "Empty sequence"

        # Count G and C occurrences
        g_count = dna.upper().count('G')
        c_count = dna.upper().count('C')
        gc_content = (g_count + c_count) / len(dna)

        formatted_gc = f"{gc_content * 100:.2f}%"

        if self.min_gc <= gc_content <= self.max_gc:
            return True, f"GC content is {formatted_gc}"
        else:
            return False, f"GC content {formatted_gc} is outside of allowed range ({self.min_gc*100}% - {self.max_gc*100}%)"