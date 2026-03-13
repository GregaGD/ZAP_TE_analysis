#!/usr/bin/env python3
import argparse, random, sys
from collections import Counter, defaultdict

# ----------- FASTA utilities ----------- #
def read_fasta(path):
    seqs = []
    cur = []
    with open(path) as fh:
        for line in fh:
            line=line.strip()
            if not line: continue
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur=[]
            else:
                cur.append(line)
        if cur: seqs.append("".join(cur))
    return seqs

def read_fasta_with_ids(path):
    """Yield (id, seq) for per-seq mode."""
    name=None; seq=[]
    with open(path) as fh:
        for line in fh:
            line=line.rstrip("\n")
            if not line: continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].strip()
                seq=[]
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)

def normalize_alphabet(seq, alphabet="RNA"):
    s = seq.upper()
    if alphabet.upper() == "RNA":
        s = s.replace("T","U")
        allowed = set("AUCGN")
    else:
        allowed = set("ATCGN")
    # keep only allowed characters (drop gaps, etc.)
    return "".join([c for c in s if c in allowed])

# ----------- counting utilities ----------- #
def count_overlapping(seq, motif):
    # counts overlapping occurrences
    k = len(motif)
    if k==0 or len(seq)<k: return 0
    c=0; i=0
    while True:
        i = seq.find(motif, i)
        if i==-1: break
        c+=1; i+=1
    return c

def count_motifs_and_bases(seqs, motifs):
    mono = Counter()
    di = Counter()
    obs = Counter()
    total_len = 0
    total_di_positions = 0
    # count observed motifs and mono/di composition
    for s in seqs:
        total_len += len(s)
        mono.update(s)
        # dinucleotide composition
        if len(s)>=2:
            total_di_positions += (len(s)-1)
            for i in range(len(s)-1):
                di[s[i:i+2]] += 1
        # motif observation
        for m in motifs:
            obs[m] += count_overlapping(s, m)
    return mono, di, obs, total_len, total_di_positions

def total_positions(seqs, k):
    return sum(max(0, len(s)-k+1) for s in seqs)

def total_positions_for_seq(seq_len, k):
    return max(0, seq_len - k + 1)

# ----------- expectations ----------- #
def mono_expected(seqs, motifs, mono_counts):
    total_bases = sum(mono_counts[b] for b in mono_counts)
    freqs = defaultdict(float)
    for b in "AUCGT":
        if total_bases>0:
            freqs[b] = mono_counts.get(b,0)/total_bases
    exp = {}
    for m in motifs:
        k=len(m)
        p=1.0
        for ch in m:
            p *= freqs[ch]
        exp[m] = p * total_positions(seqs, k)
    return exp

def dinuc_expected_tetranucleotide(seqs, motif, mono_counts, di_counts, total_di_positions):
    """
    For a tetranucleotide WXYZ, approximate:
      P(WXYZ) ≈ [P(WX) * P(XY) * P(YZ)] / [P(X) * P(Y)]
    Then multiply by number of 4-mer positions.
    Only applied if motif length==4 and we have di + mono frequencies.
    """
    k = len(motif)
    assert k==4
    # mono
    total_bases = sum(mono_counts.values())
    Pmono = {b: (mono_counts.get(b,0)/total_bases if total_bases>0 else 0.0) for b in set("AUCGT")|set(motif)}
    # di
    Pdi = {}
    for b1 in "AUCGT":
        for b2 in "AUCGT":
            d=b1+b2
            Pdi[d] = di_counts.get(d,0)/total_di_positions if total_di_positions>0 else 0.0
    W,X,Y,Z = motif[0], motif[1], motif[2], motif[3]
    numerator = Pdi[W+X]*Pdi[X+Y]*Pdi[Y+Z]
    denom = (Pmono[X]*Pmono[Y]) if (Pmono[X]>0 and Pmono[Y]>0) else 0.0
    p = (numerator/denom) if denom>0 else 0.0
    return p * total_positions(seqs, 4)

def empirical_shuffle_expected(seqs, motifs, alphabet="RNA", k_preserve=1, n_shuffles=100, seed=13):
    """
    Empirical expectation via shuffling sequences.
    k_preserve=1  -> preserves mono composition
    k_preserve=2  -> crude dinuc-preserving shuffle using block permutations (approximate)
    """
    rng = random.Random(seed)
    def shuffle_mono(s):
        arr=list(s); rng.shuffle(arr); return "".join(arr)
    def shuffle_di_block(s):
        # Quick-and-dirty: shuffle pairs then stitch; not a perfect dinuc-preserving shuffle,
        # but better than mono for many cases.
        if len(s)<2: return s
        pairs=[s[i:i+2] for i in range(0, len(s)-1, 2)]
        tail = "" if len(s)%2==0 else s[-1]
        rng.shuffle(pairs)
        t="".join(pairs)+tail
        return t[:len(s)]

    acc = Counter({m:0 for m in motifs})
    for _ in range(n_shuffles):
        shuf_counts = Counter({m:0 for m in motifs})
        for s in seqs:
            if k_preserve==1:
                ss = shuffle_mono(s)
            else:
                ss = shuffle_di_block(s)
            for m in motifs:
                shuf_counts[m] += count_overlapping(ss, m)
        for m in motifs:
            acc[m] += shuf_counts[m]
    return {m: acc[m]/n_shuffles for m in motifs}

# ----------- per-sequence (mono model) ----------- #
def per_seq_mono_table(fasta_path, motifs, alphabet="RNA"):
    """
    Print one row per FASTA record:
      seq_id len f_A f_U f_C f_G  [obs_M pos_M exp_mono_M OE_mono_M for each motif M]
    """
    # header
    dynamic_cols = []
    for m in motifs:
        dynamic_cols += [f"obs_{m}", f"pos_{m}", f"exp_mono_{m}", f"OE_mono_{m}"]
    print("\t".join(["seq_id","len","f_A","f_U","f_C","f_G"] + dynamic_cols))

    for name, raw in read_fasta_with_ids(fasta_path):
        seq = normalize_alphabet(raw, alphabet)
        n = len(seq)
        if n == 0:
            print("\t".join([name, "0", "0","0","0","0"] + ["0"]*len(dynamic_cols)))
            continue

        mono = Counter(seq)
        total_bases = mono.get("A",0)+mono.get("U",0)+mono.get("C",0)+mono.get("G",0)
        if total_bases == 0:
            fA=fU=fC=fG = (0.0,0.0,0.0,0.0)
        else:
            fA = mono.get("A",0)/total_bases
            fU = mono.get("U",0)/total_bases
            fC = mono.get("C",0)/total_bases
            fG = mono.get("G",0)/total_bases

        row_vals = []
        for m in motifs:
            k = len(m)
            obs = count_overlapping(seq, m)
            pos = total_positions_for_seq(n, k)
            # expected (mono model): product of base freqs * positions
            p = 1.0
            for ch in m:
                if ch == "A": p *= fA
                elif ch == "U": p *= fU
                elif ch == "C": p *= fC
                elif ch == "G": p *= fG
                else: p *= 0.0
            exp_mono = p * pos
            if exp_mono > 0:
                oe = obs / exp_mono
            else:
                oe = float('inf') if obs > 0 else 0.0
            row_vals.extend([str(obs), str(pos), f"{exp_mono:.6g}", f"{oe:.6g}"])

        print("\t".join([name, str(n), f"{fA:.6g}", f"{fU:.6g}", f"{fC:.6g}", f"{fG:.6g}"] + row_vals))

# ----------- main ----------- #
def main():
    ap = argparse.ArgumentParser(
        description="Motif enrichment for peaks: internal O/E, peaks vs background, and per-sequence O/E."
    )
    ap.add_argument("--peaks", required=True, help="FASTA with peak sequences")
    ap.add_argument("--background", help="FASTA with background sequences (for peaks-vs-background mode)")
    ap.add_argument("--alphabet", choices=["RNA","DNA"], default="RNA", help="Treat input as RNA or DNA (T->U for RNA)")
    ap.add_argument("--mode", choices=["internal","background"], default="internal",
                    help="internal = O/E within peaks; background = peaks vs background frequency ratio")
    ap.add_argument("--motifs", default="UA,CG,UACG", help="Comma-separated motifs (use RNA letters if --alphabet RNA)")
    ap.add_argument("--internal-model", choices=["mono","mono+di4","empirical1","empirical2"], default="mono",
                    help="internal expectation model: mono (default), mono+di4 (use dinuc model for len-4 motifs), "
                         "empirical1 (mono-preserving shuffle), empirical2 (rough di-preserving shuffle)")
    ap.add_argument("--shuffles", type=int, default=200, help="number of shuffles for empirical models")
    ap.add_argument("--seed", type=int, default=13, help="random seed for empirical models")
    ap.add_argument("--per-seq", action="store_true",
                    help="Output per-sequence table (mono model) instead of global summary (works with --mode internal).")
    args = ap.parse_args()

    motifs = [m.strip().upper() for m in args.motifs.split(",") if m.strip()]
    # Normalize motifs to chosen alphabet
    if args.alphabet.upper()=="RNA":
        motifs = [m.replace("T","U") for m in motifs]
    else:
        motifs = [m.replace("U","T") for m in motifs]

    if getattr(args, "per_seq", False):  # python3.8-safe guard
        if args.mode != "internal":
            sys.exit("Error: --per-seq is only supported with --mode internal (per-sequence O/E uses mono model).")
        per_seq_mono_table(args.peaks, motifs, args.alphabet)
        return

    # ----- GLOBAL MODES BELOW -----
    peaks_raw = [normalize_alphabet(s, args.alphabet) for s in read_fasta(args.peaks)]
    peaks = [s for s in peaks_raw if len(s)>0]

    if args.mode == "background":
        if not args.background:
            sys.exit("Error: --background FASTA is required for mode=background")
        bg_raw = [normalize_alphabet(s, args.alphabet) for s in read_fasta(args.background)]
        back = [s for s in bg_raw if len(s)>0]

        # observed and positions
        _,_,obs_peaks,_,_ = count_motifs_and_bases(peaks, motifs)
        _,_,obs_back,_,_  = count_motifs_and_bases(back,  motifs)
        # per-set total possible positions
        pos_peaks = {m: total_positions(peaks, len(m)) for m in motifs}
        pos_back  = {m: total_positions(back,  len(m)) for m in motifs}

        print("# Mode: peaks vs background (frequency ratio)")
        print("motif\tobs_peaks\tpos_peaks\tfreq_peaks\tobs_bg\tpos_bg\tfreq_bg\tenrichment_peaks_over_bg")
        for m in motifs:
            fp = (obs_peaks[m]/pos_peaks[m]) if pos_peaks[m]>0 else 0.0
            fb = (obs_back[m]/pos_back[m])   if pos_back[m]>0  else 0.0
            enr = (fp/fb) if fb>0 else float('inf') if fp>0 else 0.0
            print(f"{m}\t{obs_peaks[m]}\t{pos_peaks[m]}\t{fp:.6g}\t{obs_back[m]}\t{pos_back[m]}\t{fb:.6g}\t{enr:.6g}")
        return

    # internal mode (O/E within peaks)
    mono_p, di_p, obs_p, total_len, total_di_pos = count_motifs_and_bases(peaks, motifs)

    if args.internal_model == "mono":
        exp = mono_expected(peaks, motifs, mono_p)

    elif args.internal_model == "mono+di4":
        # mononucleotide expectation for all motifs,
        # but for length-4 motifs also compute dinuc expected and report both.
        exp = mono_expected(peaks, motifs, mono_p)
        exp_di4 = {}
        for m in motifs:
            if len(m)==4:
                exp_di4[m] = dinuc_expected_tetranucleotide(peaks, m, mono_p, di_p, total_di_pos)
        print("# Mode: internal O/E (mono model; plus dinuc-based expectation for 4-mers)")
        print("motif\tobs\tpos\tE_mono\tO/E_mono\tE_dinuc4\tO/E_dinuc4")
        for m in motifs:
            pos = total_positions(peaks, len(m))
            Emono = exp[m]
            OEM = (obs_p[m]/Emono) if Emono>0 else float('inf') if obs_p[m]>0 else 0.0
            if len(m)==4:
                Edi = exp_di4[m]
                OED = (obs_p[m]/Edi) if Edi>0 else float('inf') if obs_p[m]>0 else 0.0
                print(f"{m}\t{obs_p[m]}\t{pos}\t{Emono:.6g}\t{OEM:.6g}\t{Edi:.6g}\t{OED:.6g}")
            else:
                print(f"{m}\t{obs_p[m]}\t{pos}\t{Emono:.6g}\t{OEM:.6g}\tNA\tNA")
        return

    elif args.internal_model in ("empirical1","empirical2"):
        k_preserve = 1 if args.internal_model=="empirical1" else 2
        exp = empirical_shuffle_expected(peaks, motifs, alphabet=args.alphabet,
                                         k_preserve=k_preserve, n_shuffles=args.shuffles, seed=args.seed)
        print(f"# Mode: internal O/E (empirical shuffle, k_preserve={k_preserve}, n={args.shuffles})")
        print("motif\tobs\tE_empirical\tO/E_empirical")
        for m in motifs:
            E = exp[m]
            OE = (obs_p[m]/E) if E>0 else float('inf') if obs_p[m]>0 else 0.0
            print(f"{m}\t{obs_p[m]}\t{E:.6g}\t{OE:.6g}")
        return

    else:
        # default mono
        pass

    # If we get here, internal_model == mono
    print("# Mode: internal O/E (mononucleotide model)")
    print("motif\tobs\tpos\tE_mono\tO/E_mono")
    for m in motifs:
        pos = total_positions(peaks, len(m))
        E = exp[m]
        OE = (obs_p[m]/E) if E>0 else float('inf') if obs_p[m]>0 else 0.0
        print(f"{m}\t{obs_p[m]}\t{pos}\t{E:.6g}\t{OE:.6g}")

if __name__ == "__main__":
    main()
