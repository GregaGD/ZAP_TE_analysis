Motif Enrichment Pipeline – Manual
==================================

Overview
--------
This tool calculates enrichment of nucleotide motifs (e.g. UA, CG, UACG) 
in RNA/DNA sequences.  
It supports two main modes:

1. Internal expectation (within peaks):
   Compares observed motif counts in peaks vs. expected counts based on
   their own nucleotide composition (or shuffled versions).
   → Answers: “Are my peaks enriched for this motif compared to
   randomised versions of the same peaks?”

2. Background comparison (peaks vs mRNAs):
   Compares motif frequencies in peaks vs. frequencies in a background set
   (e.g. random mRNA windows).
   → Answers: “Are my peaks enriched for this motif compared to normal mRNA?”

Input
-----
--peaks       FASTA file with your peak sequences
--background  FASTA file with background sequences (only for background mode)
--motifs      Comma-separated motifs (default: UA,CG,UACG)
--alphabet    RNA (default; converts T→U) or DNA

Output
------
Tab-separated text with:
- Observed counts
- Possible motif positions
- Expected counts (internal mode)
- Frequencies (background mode)
- O/E or enrichment ratios

Internal Expectation Mode
-------------------------

1. Mononucleotide model (--internal-model mono)
   - Observed (O): motif counts in peaks
   - Expected (E):
       f(A) = #A / total bases, same for U,C,G
       P(motif) = product of base frequencies
       E(m) = P(m) × total possible positions
   - O/E = O / E

2. Dinucleotide model for tetranucleotides (--internal-model mono+di4)
   - For 4-mers WXYZ:
       P(WXYZ) ≈ [P(WX) × P(XY) × P(YZ)] / [P(X) × P(Y)]
       E = P(WXYZ) × positions
   - Output shows both mono and di expectations.

3. Empirical shuffle models (--internal-model empirical1 or empirical2)
   - Shuffle each sequence many times
   - empirical1: mono-preserving shuffle
   - empirical2: rough di-preserving shuffle
   - Expected = average motif count across shuffles
   - O/E = observed / expected

Background Comparison Mode (--mode background)
----------------------------------------------
- Requires peaks + background FASTA
- For each dataset:
    O = observed motif count
    Positions = sum(len(seq)-k+1)
    Frequency = O / Positions
- Enrichment ratio = freq(peaks) / freq(background)
- >1 = enriched in peaks, <1 = depleted

Example Commands
----------------
# Internal O/E (mono model)
python3 motif_enrichment.py --peaks peaks.fa --mode internal --internal-model mono

# Internal O/E with mono+di model for tetranucleotides
python3 motif_enrichment.py --peaks peaks.fa --mode internal --internal-model mono+di4

# Internal O/E with empirical shuffle (mono-preserving, 500 shuffles)
python3 motif_enrichment.py --peaks peaks.fa --mode internal --internal-model empirical1 --shuffles 500

# Peaks vs background (RNA alphabet)
python3 motif_enrichment.py --peaks peaks.fa --background background.fa --mode background --alphabet RNA

Practical Notes
---------------
- Match background sequence lengths and context to peaks
- Script outputs counts/frequencies; apply stats tests separately if needed
- Use --alphabet DNA if sequences contain T
