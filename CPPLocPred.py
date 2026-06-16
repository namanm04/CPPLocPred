#!/usr/bin/env python
# CPPLocPred: Prediction of Cell Penetrating Peptides and Subcellular Localization
# Usage: python CPPLocPred.py -i input.fasta -o output.csv -j 1
#        python CPPLocPred.py -i input.fasta -o output.csv -j 2 -l Cytoplasm -c None
# motifs/ directory layout: motifs/<class>/<Location>_motif
#   classes: None, Koolman, Betts-Russell, Rasmol
#   locations: Cytoplasm, Nucleus, Mitochondria, Endo_lysosome, Others

import os, sys, re, pickle, argparse
import pandas as pd

# ── Constants ──────────────────────────────────────────────
STD_AA = "ACDEFGHIKLMNPQRSTVWY"

AAC_COLS = ["AAC_"+a for a in STD_AA]
DDR_COLS = ["DDR_"+a for a in STD_AA]

LOCATIONS = ["Cytoplasm", "Nucleus", "Mitochondria", "Endo_lysosome", "Others"]

CLASSES = ["None", "Koolman", "Betts-Russell", "Rasmol"]



# ── FASTA reader ───────────────────────────────────────────
def read_fasta(path):
    seqs = []
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, "".join(seq)))
                header, seq = line[1:], []
            else:
                seq.append(line.upper())
    if header is not None:
        seqs.append((header, "".join(seq)))
    return seqs

# ── Feature extraction ─────────────────────────────────────
def aac(seq):
    n = len(seq)
    return [round(seq.count(a) / n * 100, 2) for a in STD_AA]

def ddr(seq):
    feats = []
    rev = seq[::-1]
    for a in STD_AA:
        pos     = [i for i, c in enumerate(seq) if c == a]
        rev_pos = [i for i, c in enumerate(rev)  if c == a]
        gaps = []
        for i in range(len(pos) - 1):
            gaps.append(pos[i+1] - pos[i] - 1)
        if pos:
            gaps.insert(0, pos[0])
            gaps.append(rev_pos[0])
        cc1 = sum(gaps) + 1
        cc  = sum(x * x for x in gaps)
        feats.append(round(cc / cc1, 2))
    return feats

def validate(seq):
    bad = set(seq) - set(STD_AA)
    return (False, bad) if bad else (True, None)

# ══════════════════════════════════════════════════════════
# JOB 1 — PREDICTION
# ══════════════════════════════════════════════════════════
def run_prediction(args):
    model_dir = args.model_dir

    # load models
    def load(f): 
        with open(os.path.join(model_dir, f), "rb") as fh:
            return pickle.load(fh)

    cpp_model    = load("CPP_vs_NonCPP_ET_AAC.pkl")
    cyto_model   = load("Cytoplasm_DDR_CatBoost.pkl")
    nuc_model    = load("Nucleus_DDR_CatBoost.pkl")
    mito_model   = load("Mitochondria_DDR_CatBoost.pkl")
    endo_model   = load("Endo_lyso_DDR_CatBoost.pkl")
    other_model  = load("Other_DDR_CatBoost.pkl")
    loc_thresh   = load("Localization_thresholds.pkl")

    cpp_threshold = float(args.threshold)
    seqs = read_fasta(args.input)
    print("Sequences loaded: {}".format(len(seqs)))

    rows = []
    for sid, seq in seqs:
        ok, bad = validate(seq)
        if not ok:
            rows.append([sid, seq, "Invalid AA: {}".format(",".join(sorted(bad))),
                         "Invalid", "", "", "", "", "", "NA"])
            continue

        # Stage 1 — CPP vs Non-CPP
        cpp_prob = cpp_model.predict_proba(
            pd.DataFrame([aac(seq)], columns=AAC_COLS))[0][1]

        if cpp_prob < cpp_threshold:
            rows.append([sid, seq, round(cpp_prob,4),
                         "Non-CPP", "", "", "", "", "", "NA"])
            continue

        # Stage 2 — localization
        ddr_df = pd.DataFrame([ddr(seq)], columns=DDR_COLS)
        probs = {
            "Cytoplasm":    cyto_model.predict_proba(ddr_df)[0][1],
            "Nucleus":      nuc_model.predict_proba(ddr_df)[0][1],
            "Mitochondria": mito_model.predict_proba(ddr_df)[0][1],
            "Endo_lysosome":endo_model.predict_proba(ddr_df)[0][1],
            "Others":       other_model.predict_proba(ddr_df)[0][1],
        }

        predicted = [loc for loc in LOCATIONS
                     if probs[loc] >= loc_thresh.get(loc, 0.5)]

        if not predicted:
            best = max(probs, key=probs.get)
            location = "Uncertain (Best={}, Prob={:.3f})".format(best, probs[best])
        else:
            location = ";".join(predicted)

        rows.append([sid, seq,
                     round(cpp_prob, 4), "CPP",
                     round(probs["Cytoplasm"],    4),
                     round(probs["Nucleus"],       4),
                     round(probs["Mitochondria"],  4),
                     round(probs["Endo_lysosome"], 4),
                     round(probs["Others"],        4),
                     location])

    cols = ["ID","Sequence","CPP_Probability","CPP_Prediction",
            "Cytoplasm_Probability","Nucleus_Probability",
            "Mitochondria_Probability","Endo_lysosome_Probability",
            "Others_Probability","Final_Localization"]

    pd.DataFrame(rows, columns=cols).to_csv(args.output, index=False)
    print("Results saved to: {}".format(args.output))

# ══════════════════════════════════════════════════════════
# JOB 2 — MOTIF SEARCH
# ══════════════════════════════════════════════════════════
def parse_merci_output(merci_file):
    """Parse MERCI motif locator output into a hit list."""
    hits = {}    # seq_id -> [ {motif, start, end, matched}, ... ]
    coverage = {}  # seq_id -> total motif match count

    cur_motif  = ""
    cur_seq_id = ""
    cur_pos    = None
    in_coverage = False

    with open(merci_file) as f:
        for line in f:
            line_s = line.strip()

            if line_s == "COVERAGE":
                in_coverage = True
                continue

            if in_coverage:
                # e.g.  seq_20|label=1 (6 motifs match)
                m = re.match(r'^(\S+)\s+\((\d+)\s+motifs? match\)', line_s)
                if m:
                    coverage[m.group(1)] = int(m.group(2))
                continue

            if re.match(r'^\*+$', line_s):
                continue

            if line_s.startswith("MOTIF:"):
                cur_motif  = line_s[6:].strip()
                cur_seq_id = ""
                cur_pos    = None
                continue

            if line_s.startswith(">"):
                cur_seq_id = line_s[1:].strip()
                cur_pos    = None
                continue

            if re.match(r'^\(sequence \d+\)', line_s):
                continue

            if re.search(r'sequences? contain the motif', line_s):
                continue

            m = re.match(r'^position\(s\):\s*(.+)', line_s)
            if m:
                cur_pos = int(m.group(1).split(",")[0].strip())
                continue

            m = re.match(r'^motif\(s\):\s*(\S+)', line_s)
            if m and cur_seq_id and cur_pos is not None:
                matched = m.group(1)
                start   = cur_pos
                end     = cur_pos + len(matched) - 1
                hits.setdefault(cur_seq_id, []).append({
                    "motif":   cur_motif,
                    "start":   start,
                    "end":     end,
                    "matched": matched,
                })
                cur_pos = None

    return hits, coverage



def run_motif_search(args):
    location   = args.location
    cls        = args.cls
    motif_dir  = args.motif_dir

    # resolve motif file: motif_dir/class/Location_motif
    motif_file = os.path.join(motif_dir, cls, "{}_motif".format(location))
    if not os.path.exists(motif_file):
        sys.exit("Motif file not found: {}".format(motif_file))

    seqs = read_fasta(args.input)
    print("Sequences loaded: {}".format(len(seqs)))

    # write tmp fasta for MERCI
    tmp_fasta = args.output + ".tmp.fasta"
    merci_out = args.output + ".tmp.merci"

    with open(tmp_fasta, "w") as f:
        for sid, seq in seqs:
            f.write(">{}\n{}\n".format(sid, seq))

    # run MERCI
    merci_pl = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "MERCI_motif_locator.pl")
    if not os.path.exists(merci_pl):
        sys.exit("MERCI_motif_locator.pl not found next to this script.")

    cmd = '{} {} -p {} -i {} -o {}'.format(
        args.perl, merci_pl, tmp_fasta, motif_file, merci_out)
    ret = os.system(cmd)
    if ret != 0 or not os.path.exists(merci_out):
        sys.exit("MERCI failed. Return code: {}".format(ret))

    # parse MERCI output
    hit_map, coverage = parse_merci_output(merci_out)

    # ── Summary table: one row per sequence ──
    summary_rows = []
    for sid, seq in seqs:
        hits       = hit_map.get(sid, [])
        total_hits = len(hits)
        motif_pats = ";".join(sorted(set(h["motif"]   for h in hits))) if hits else ""
        positions  = ";".join("{}-{}".format(h["start"], h["end"]) for h in hits)
        matched    = ";".join(h["matched"] for h in hits)
        summary_rows.append([sid, seq, len(seq), total_hits, motif_pats, positions, matched])

    summary_cols = ["ID", "Sequence", "Length",
                    "Total_Hits", "Motif_Patterns",
                    "Hit_Positions", "Matched_Residues"]

    # ── Detailed table: one row per hit ──
    detail_rows = []
    for sid, seq in seqs:
        hits = hit_map.get(sid, [])
        if not hits:
            detail_rows.append([sid, seq, 0, "", "", "", ""])
        else:
            for h in hits:
                detail_rows.append([
                    sid, seq,
                    h["start"], h["end"],
                    h["end"] - h["start"] + 1,
                    h["motif"],
                    h["matched"],
                ])

    detail_cols = ["ID", "Sequence", "Start", "End",
                   "Hit_Length", "Motif_Pattern", "Matched_Residues"]

    # write both sheets to one CSV (summary first, then blank line, then detail)
    summary_df = pd.DataFrame(summary_rows, columns=summary_cols)
    detail_df  = pd.DataFrame(detail_rows,  columns=detail_cols)

    summary_df.to_csv(args.output, index=False)

    detail_file = args.output.replace(".csv", "_hits.csv")
    detail_df.to_csv(detail_file, index=False)

    print("Summary saved to : {}".format(args.output))
    print("Hit detail saved to: {}".format(detail_file))

    # print terminal summary
    print("\n{:<25} {:>6} {:>12}".format("Sequence ID", "Length", "Total Hits"))
    print("-" * 46)
    for row in summary_rows:
        print("{:<25} {:>6} {:>12}".format(row[0], row[2], row[3]))
    print("-" * 46)
    print("Location scanned: {}".format(location))
    print("Motif class     : {}".format(cls))

    # cleanup
    for f in [tmp_fasta, merci_out]:
        if os.path.exists(f):
            os.remove(f)

# ── Argument parser ────────────────────────────────────────
def build_parser():
    p = argparse.ArgumentParser(
        description="CPPLocPred: Cell Penetrating Peptide prediction and motif search",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Examples:
  # Job 1 — prediction
  python CPPLocPred.py -i input.fasta -o results.csv -j 1 -t 0.44

  # Job 2 — motif search (Nucleus)
  python CPPLocPred.py -i input.fasta -o motifs.csv -j 2 -l Nucleus -c Koolman
""")

    p.add_argument("-i", "--input",     required=True,  help="Input FASTA file")
    p.add_argument("-o", "--output",    required=True,  help="Output CSV file")
    p.add_argument("-j", "--job",       type=int, default=1, choices=[1,2],
                   help="Job type: 1=Prediction (default), 2=Motif Search")

    # prediction options
    p.add_argument("-t", "--threshold", default="0.44",
                   help="CPP probability threshold (default: 0.44)")
    p.add_argument("-m", "--model_dir", default="./",
                   help="Directory containing model .pkl files (default: ./)")

    # motif search options
    p.add_argument("-l", "--location",  default="Cytoplasm",
                   choices=LOCATIONS,
                   help="Subcellular location for motif search (default: Cytoplasm)")
    p.add_argument("-c", "--cls",       default="None",
                   choices=CLASSES,
                   help="Motif class subfolder: None, Koolman, Betts-Russell, Rasmol (default: None)")
    p.add_argument("--motif_dir",       default="./motifs",
                   help="Directory containing motif class subfolders (default: ./motifs)")
    p.add_argument("--perl",            default="perl",
                   help="Path to perl executable (default: perl)")

    return p

# ── Main ───────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("  CPPLocPred: CPP Prediction & Subcellular Localization")
    print("=" * 60)

    parser = build_parser()
    args   = parser.parse_args()

    if not os.path.exists(args.input):
        sys.exit("Input file not found: {}".format(args.input))

    if args.job == 1:
        print("Job: Prediction")
        run_prediction(args)
    elif args.job == 2:
        print("Job: Motif Search  |  Location: {}  |  Class: {}".format(
            args.location, args.cls))
        run_motif_search(args)

if __name__ == "__main__":
    main()