#!/usr/bin/env python3
"""Proton-linkage scoring of pH-dependent binding.

For every ionizable site, PROPKA gives a pKa in the free state and one in the
complex. A site only contributes pH dependence when those differ, and the
contribution follows the Wyman linkage relation

    dG_site(pH) = -RT * ln[ (1 + 10^(pKa_bound - pH)) / (1 + 10^(pKa_free - pH)) ]

The switch strength is the difference between the two assay pH values,

    ddG_switch = dG(pH_on) - dG(pH_off)

which is negative when binding is stronger at pH_on. One fully coupled site can
shift binding by at most 2.303*RT*dpH = 1.23 kcal/mol over 6.5 -> 7.4, i.e. ~8x
in affinity, so a design that must lose detectable binding at pH 7.4 needs three
or more coupled sites.

    python ph_linkage.py --complex design.pdb --binder B --target A
"""
import argparse, math, re, subprocess, sys, tempfile, os, csv

RT = 0.0019872 * 298.15          # kcal/mol
TITRATABLE = {"ASP", "GLU", "HIS", "LYS", "ARG", "CYS", "TYR"}


def write_chains(pdb, chains, out):
    with open(out, "w") as fh:
        for line in open(pdb):
            if line.startswith(("ATOM", "HETATM")) and line[21] in chains:
                fh.write(line)
        fh.write("END\n")
    return out


def run_propka(pdb):
    # propka writes <stem>.pka into the working directory, not beside the input
    work = os.path.dirname(os.path.abspath(pdb))
    subprocess.run([sys.executable, "-m", "propka", os.path.basename(pdb)], cwd=work,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
    pka_file = os.path.join(work, os.path.splitext(os.path.basename(pdb))[0] + ".pka")
    if not os.path.exists(pka_file):
        raise RuntimeError("PROPKA produced no output for %s" % pdb)
    out = {}
    pat = re.compile(r"\s+(%s)\s+(\d+)\s+(\w)\s+([\d.\-]+)\s+([\d.\-]+)" % "|".join(TITRATABLE))
    for line in open(pka_file):
        m = pat.match(line)
        if m:
            out[(m.group(1), int(m.group(2)), m.group(3))] = float(m.group(4))
    return out


def dg_site(pka_free, pka_bound, pH):
    return -RT * math.log((1 + 10 ** (pka_bound - pH)) / (1 + 10 ** (pka_free - pH)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complex", required=True, help="PDB of the binder-target complex")
    ap.add_argument("--binder", required=True, help="binder chain id(s)")
    ap.add_argument("--target", required=True, help="target chain id(s)")
    ap.add_argument("--ph-on", type=float, default=6.5)
    ap.add_argument("--ph-off", type=float, default=7.4)
    ap.add_argument("--min-shift", type=float, default=0.3,
                    help="ignore sites whose pKa moves less than this")
    ap.add_argument("--out", help="write per-site rows to this CSV")
    a = ap.parse_args()

    tmp = tempfile.mkdtemp()
    bound = run_propka(write_chains(a.complex, set(a.binder + a.target), os.path.join(tmp, "cx.pdb")))
    free = {}
    free.update(run_propka(write_chains(a.complex, set(a.binder), os.path.join(tmp, "binder.pdb"))))
    free.update(run_propka(write_chains(a.complex, set(a.target), os.path.join(tmp, "target.pdb"))))

    rows, total = [], 0.0
    for key, pka_free in free.items():
        pka_bound = bound.get(key)
        if pka_bound is None or abs(pka_bound - pka_free) < a.min_shift:
            continue
        ddg = dg_site(pka_free, pka_bound, a.ph_on) - dg_site(pka_free, pka_bound, a.ph_off)
        total += ddg
        res, num, ch = key
        rows.append(dict(residue="%s%d" % (res, num), chain=ch, side="binder" if ch in a.binder else "target",
                         pka_free=round(pka_free, 2), pka_bound=round(pka_bound, 2),
                         d_pka=round(pka_bound - pka_free, 2), ddg_switch=round(ddg, 3)))
    rows.sort(key=lambda r: r["ddg_switch"])

    print("%-10s %-7s %-9s %-9s %-7s %s" % ("residue", "side", "pKa_free", "pKa_bound", "dpKa", "ddG_switch"))
    for r in rows:
        print("%-10s %-7s %-9.2f %-9.2f %-7.2f %+.3f" %
              (r["residue"], r["side"], r["pka_free"], r["pka_bound"], r["d_pka"], r["ddg_switch"]))
    ratio = math.exp(-total / RT)
    print("\nddG_switch total : %+.3f kcal/mol" % total)
    print("affinity ratio   : %.1fx %s at pH %.1f" % (max(ratio, 1 / ratio),
          "stronger" if total < 0 else "weaker", a.ph_on))
    print("coupled sites    : %d (one site caps at %.2f kcal/mol over this window)"
          % (len(rows), 2.303 * RT * abs(a.ph_off - a.ph_on)))
    print("verdict          : %s" % ("PASS" if total <= -1.5 else "FAIL - not a usable switch"))

    if a.out:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else
                               ["residue", "chain", "side", "pka_free", "pka_bound", "d_pka", "ddg_switch"])
            w.writeheader(); w.writerows(rows)
        print("wrote %s" % a.out)


if __name__ == "__main__":
    main()
