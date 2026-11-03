#!/usr/bin/env python3
"""Build one-seed AlphaFold Server jobs, FASTA, and native/A-B chain map."""
import csv, json
from pathlib import Path
from Bio.PDB import PDBParser, is_aa
from Bio.SeqUtils import seq1

ROOT = Path(__file__).resolve().parents[1]
SEL = ROOT / "bench20/selection.tsv"
PDBDIR = ROOT / "pdbbind/P-P"
OUT = ROOT / "alphafold3_server"
STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


def atom_sequences(path):
    model = next(PDBParser(QUIET=True).get_structure("s", str(path)).get_models())
    out = {}
    for chain in model:
        residues = [r for r in chain if is_aa(r, standard=True)]
        if residues:
            out[chain.id] = "".join(seq1(r.resname) for r in residues)
    return out


def seqres_sequences(path):
    names = {}
    with path.open() as f:
        for line in f:
            if line.startswith("SEQRES"):
                chain = line[11].strip() or " "
                names.setdefault(chain, []).extend(line[19:70].split())
    return {c: "".join(seq1(x) for x in rs) for c, rs in names.items()}


def main():
    jobs = OUT / "jobs/bench20_ready"
    jobs.mkdir(parents=True, exist_ok=True)
    map_path = OUT / "bench20_chain_mapping.tsv"
    fasta_path = OUT / "bench20_sequences.fasta"
    records, mapping = [], []
    with SEL.open() as f:
        selected = list(csv.DictReader(f, delimiter="\t"))
    for row in selected:
        cid = row["id"]
        pdb = PDBDIR / f"{cid}_complex.pdb"
        atom, seqres = atom_sequences(pdb), seqres_sequences(pdb)
        seqs=[]
        for i, key in enumerate(("chain1", "chain2")):
            native = row[key]
            seq = seqres.get(native) or atom.get(native)
            source = "SEQRES" if native in seqres and seqres[native] else "ATOM"
            if seq and any(a not in STANDARD for a in seq) and atom.get(native):
                seq, source = atom[native], "ATOM (SEQRES nonstandard)"
            if not seq or any(a not in STANDARD for a in seq):
                raise SystemExit(f"{cid} chain {native}: absent or nonstandard {source} sequence")
            if native in atom and seqres.get(native) and seqres[native] != atom[native]:
                discrepancy = f"DIFF seqres={len(seqres[native])} atom={len(atom[native])}"
            else:
                discrepancy = "match" if native in seqres else "SEQRES absent; ATOM used"
            seqs.append(seq)
            af = "A" if i == 0 else "B"
            records.append((f">{cid}_{af} native_chain={native} source={source}", seq))
            mapping.append([cid, i+1, native, af, len(atom.get(native, "")), len(seqres.get(native, "")), source, discrepancy, str(pdb.relative_to(ROOT))])
        payload = [{"name": f"PDBBind_{cid}_bench20", "modelSeeds": [20261030],
                    "sequences": [{"proteinChain": {"sequence": seqs[0], "count": 1}},
                                  {"proteinChain": {"sequence": seqs[1], "count": 1}}],
                    "dialect": "alphafoldserver", "version": 1}]
        dest = jobs / f"{cid}.json"
        if dest.exists():
            raise SystemExit(f"Refusing to overwrite existing job: {dest}")
        dest.write_text(json.dumps(payload, indent=2) + "\n")
    fasta_path.write_text("\n".join(h + "\n" + s for h, s in records) + "\n")
    with map_path.open("w", newline="") as f:
        w=csv.writer(f, delimiter="\t")
        w.writerow(["complex","partner_order","native_chain","af3_chain","atom_length","seqres_length","sequence_source","sequence_check","source_pdb"])
        w.writerows(mapping)
    print(f"Wrote {len(selected)} jobs to {jobs}; FASTA {fasta_path}; mapping {map_path}")

if __name__ == "__main__": main()
