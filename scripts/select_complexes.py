#!/usr/bin/env python3
"""Select a deterministic PDBBind PP subset with exact measured Kd entries."""
import argparse
import csv
import math
import re
import json
import urllib.request
from pathlib import Path

from Bio.PDB import PDBParser, is_aa
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT = Path(__file__).resolve().parents[1]
INDEX = ROOT / "pdbbind/index/INDEX_general_PP.2020R1.lst"
STRUCTURES = ROOT / "pdbbind/P-P"
STANDARD = set("ACDEFGHIKLMNPQRSTVWY")
SEED = 20261030
ALIGNER = PairwiseAligner()
ALIGNER.mode = "global"
ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
ALIGNER.open_gap_score = -10
ALIGNER.extend_gap_score = -0.5
KD_RE = re.compile(r"\bKd=(\d+(?:\.\d+)?)(fM|pM|nM|uM|mM|M)\b", re.I)
FACTORS = {"fm": 1e-15, "pm": 1e-12, "nm": 1e-9, "um": 1e-6, "mm": 1e-3, "m": 1.0}


def chain_info(path):
    structure = PDBParser(QUIET=True).get_structure("x", str(path))
    model = next(structure.get_models())
    info = {}
    for chain in model:
        amino_acid_like = [r for r in chain if is_aa(r, standard=False)]
        if any(not is_aa(r, standard=True) for r in amino_acid_like):
            return None
        residues = amino_acid_like
        if not residues:
            continue
        sequence = "".join(seq1(r.resname, custom_map={}) for r in residues)
        if any(x not in STANDARD for x in sequence):
            return None
        info[chain.id] = (len(residues), sequence)
    return info


def parse_index(line):
    m = re.match(r"^(\w{4})\s+(\S+)\s+(\d{4})\s+(\S+)\s+//\s+(.+)$", line)
    if not m:
        return None
    pdb, resolution, year, kd, comment = m.groups()
    km = KD_RE.search(kd)
    if not km or re.search(r"[<>~]|\b(?:Ki|IC50)\s*=", line, re.I):
        return None
    try:
        resolution = float(resolution)
    except ValueError:
        return None
    if resolution > 3.0:
        return None
    if re.search(r"incomplete ligand|different protein in assay|nonstandard assay|\bmutant\b|\bmutation\b|variant|\bnot for\b|\buntreated\b|\btreated\b|\?{1,}|\bKd\s*=\s*\d+(?:\.\d+)?\s*[-–]\s*\d|\bnot\s+\d|\bcode\s+[a-z0-9]{4}\b", comment, re.I):
        return None
    # Exclude comments reporting multiple alternative Kd values.
    if len(re.findall(r"\bKd\s*=", comment, re.I)) > 1:
        return None
    pdf = re.search(r"([^ /]+)\.pdf", comment, re.I)
    parent = pdf.group(1).lower() if pdf else pdb.lower()
    # The index partner notation is the first parenthesized pair: (chain(s)|chain(s)).
    pm = re.search(r"\(([^()|]+)\|([^()|]+)\)", comment)
    if not pm:
        return None
    a, b = pm.groups()
    if len(a) != 1 or len(b) != 1 or a == b:
        return None
    kd_m = float(km.group(1)) * FACTORS[km.group(2).lower()]
    return dict(id=pdb.lower(), parent=parent, resolution=resolution, year=int(year), kd=kd,
                kd_m=kd_m, logkd=math.log10(kd_m), comment=comment,
                chain1=a, chain2=b)


def rcsb_method(pdb_id):
    """Read experimental method from the RCSB core-entry API (standard library only)."""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    request = urllib.request.Request(url, headers={"User-Agent": "bench20-selector/1.0"})
    with urllib.request.urlopen(request, timeout=15) as response:
        payload = json.load(response)
    return payload.get("rcsb_entry_info", {}).get("experimental_method")


def same_family(seq_a, seq_b):
    """Flag close homologs: >=70% identity over >=70% of the shorter sequence."""
    coordinates = ALIGNER.align(seq_a, seq_b)[0].coordinates
    identical = aligned = 0
    for j in range(coordinates.shape[1] - 1):
        a0, a1 = map(int, coordinates[0, j:j+2])
        b0, b1 = map(int, coordinates[1, j:j+2])
        if a1 > a0 and b1 > b0:
            n = min(a1-a0, b1-b0)
            aligned += n
            identical += sum(seq_a[a0+k] == seq_b[b0+k] for k in range(n))
    return aligned / min(len(seq_a), len(seq_b)) >= 0.70 and identical / max(aligned, 1) >= 0.70


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--number", type=int, default=20)
    ap.add_argument("--index", type=Path, default=INDEX)
    ap.add_argument("--structures", type=Path, default=STRUCTURES)
    args = ap.parse_args()
    candidates = []
    for line in args.index.read_text(errors="replace").splitlines():
        row = parse_index(line)
        if not row:
            continue
        path = args.structures / f"{row['id']}_complex.pdb"
        if not path.is_file():
            continue
        chains = chain_info(path)
        if not chains or set(chains) != {row['chain1'], row['chain2']}:
            continue
        l1, s1 = chains[row['chain1']]
        l2, s2 = chains[row['chain2']]
        if min(l1, l2) <= 50 or l1 + l2 > 600:
            continue
        row.update(len1=l1, len2=l2, seq1=s1, seq2=s2)
        candidates.append(row)
    # Deterministic spread across affinity range, followed by fixed-seed shuffle for ties.
    import random
    rng = random.Random(SEED)
    rng.shuffle(candidates)
    candidates.sort(key=lambda x: x["logkd"])
    if len(candidates) < args.number:
        raise SystemExit(f"Only {len(candidates)} structures pass mechanical filters; need {args.number}.")
    # Greedily maximize minimum affinity separation; prefer unique PDBs (already unique).
    chosen = [candidates[0], candidates[-1]]
    used_parents = {x["parent"] for x in chosen}
    rest = [x for x in candidates if x not in chosen and x["parent"] not in used_parents]
    while len(chosen) < args.number:
        options = [x for x in rest if x["parent"] not in used_parents
                   and all(not same_family(sa, sb)
                           for y in chosen
                           for sa in (x["seq1"], x["seq2"])
                           for sb in (y["seq1"], y["seq2"]))]
        if not options:
            raise SystemExit(f"Only {len(chosen)} complexes pass filters with unique parent PDB references and no >=70% identity family match; need {args.number}.")
        nxt = max(options, key=lambda x: min(abs(x['logkd'] - y['logkd']) for y in chosen))
        chosen.append(nxt); used_parents.add(nxt["parent"])
        rest.remove(nxt)
    chosen.sort(key=lambda x: x["id"])
    out = ROOT / "bench20"
    method_cache = out / "experimental_methods.tsv"
    cached = {}
    if method_cache.exists():
        with method_cache.open() as f:
            cached = {r["id"]: r["experimental_method"] for r in csv.DictReader(f, delimiter="\t")}
    methods = {}
    for row in chosen:
        try:
            methods[row["id"]] = cached.get(row["id"]) or rcsb_method(row["id"])
        except Exception as exc:
            raise SystemExit(f"Could not verify X-ray method for {row['id']} via RCSB: {exc}")
        if str(methods[row["id"]]).lower() != "x-ray":
            raise SystemExit(f"{row['id']} has experimental method {methods[row['id']]!r}, not X-ray; revise the selection before continuing.")
    with method_cache.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["id", "experimental_method", "source"])
        w.writerows((row["id"], methods[row["id"]], "RCSB PDB core entry API") for row in chosen)

    out.mkdir(exist_ok=True)
    cols = ["id", "chain1", "chain2", "len1", "len2", "resolution", "year", "Kd", "Kd_M", "log10_Kd", "comment"]
    with (out / "selection.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for x in chosen:
            w.writerow({"id": x["id"], "chain1": x["chain1"], "chain2": x["chain2"],
                        "len1": x["len1"], "len2": x["len2"], "resolution": x["resolution"],
                        "year": x["year"], "Kd": x["kd"], "Kd_M": f"{x['kd_m']:.8g}",
                        "log10_Kd": f"{x['logkd']:.8f}", "comment": x["comment"]})
    for filename in ("comp_list", "clone_list"):
        (ROOT / filename).write_text("".join(x["id"] + "\n" for x in chosen))
    with (ROOT / "expt_kd.csv").open("w", newline="") as f:
        w = csv.writer(f); w.writerow(["complex", "logKd"])
        w.writerows((x["id"], f"{x['logkd']:.8f}") for x in chosen)
    print(f"Selected {len(chosen)}; logKd span {min(x['logkd'] for x in chosen):.2f} to {max(x['logkd'] for x in chosen):.2f}")

if __name__ == "__main__":
    main()
