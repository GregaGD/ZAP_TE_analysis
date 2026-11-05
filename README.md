# ZAP_TE_analysis

Analysis scripts for studying the antiviral protein **ZAP** (zinc-finger antiviral protein), which recognises CpG-rich RNA. The main tool, [`motif_enrichment.py`](scripts/motif_enrichment.py), tests whether a set of sequences, such as RNA-binding sites (e.g. iCLIP peaks), is enriched for short motifs like CpG and UpA.

## motif_enrichment.py

Counts overlapping occurrences of nucleotide motifs (default: `UA`, `CG`, `UACG`) and compares them with an expectation:

| Mode | Compares observed counts with… | Answers |
|---|---|---|
| `--mode internal` | expected counts from the sequences' own composition (mononucleotide, dinucleotide or shuffle-based models) | "Is this motif more common than the base composition alone would predict?" |
| `--mode background` | motif frequency in a background set (e.g. random mRNA windows) | "Is this motif more common in my sites than in typical sequences?" |
| `--per-seq` | per-sequence observed/expected values (mononucleotide model) | "Which individual sequences drive the enrichment?" |

Output is a tab-separated table with observed counts, possible positions, expected counts or frequencies, and observed/expected (O/E) or enrichment ratios.

### Requirements

Python 3 only. The script uses the standard library, so there is nothing to install.

### Quick start

```bash
# Observed/expected within your sequences (mononucleotide model)
python3 scripts/motif_enrichment.py --peaks peaks.fa --mode internal --internal-model mono

# Shuffle-based expectation (mononucleotide-preserving, 500 shuffles)
python3 scripts/motif_enrichment.py --peaks peaks.fa --mode internal --internal-model empirical1 --shuffles 500

# Sites vs background
python3 scripts/motif_enrichment.py --peaks peaks.fa --background background.fa --mode background

# Custom motifs, one row per sequence
python3 scripts/motif_enrichment.py --peaks peaks.fa --motifs CG,UACG,CGCG --per-seq > per_sequence.tsv
```

Input sequences can be DNA or RNA: with the default `--alphabet RNA`, T is converted to U.

See the **[full manual](scripts/README.md)** for how each expectation model is calculated and for practical notes on choosing a background set.

## Related publication

Bouton CR†, Gimpelj Domjanič G†, Lista MJ†, et al. **The long isoform of ZAP coordinates multiple enzymes to mediate complete decay of target transcripts.** *Cell Reports* 45(5):117305 (2026). [doi:10.1016/j.celrep.2026.117305](https://doi.org/10.1016/j.celrep.2026.117305)  
† co-first authors
