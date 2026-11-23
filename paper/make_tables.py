"""Tables III and IV and the numbers quoted in the text of the paper.

Reads only the metrics already computed by the 2026-10-05 campaign
(``STNs-MOCO-MoWFLOP/sharing/metricas_stn_2026-11-17/*_todas.csv``, the
concatenation of the campaign metric CSVs plus a ``paper_no`` column) and the
partition report of the Shannon STNs of instance 4.

Rounding: the values printed in the paper were rounded half-up to three
decimals and then half-up to the precision shown (e.g. 0.0148 -> 0.015 -> 0.02).
``fmt`` reproduces that rule; ``--direct`` rounds each value once instead and
also lists the cells where the two rules differ.

Usage (from any folder):
    python paper/make_tables.py [--direct]
"""
import json
import sys
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
MET = REPO / "STNs-MOCO-MoWFLOP" / "sharing" / "metricas_stn_2026-11-17"
REPORT = (REPO / "STNs-MOCO-MoWFLOP" / "sharing" / "campanhas"
          / "full_ext_entropy_sparse_2026-10-05" / "reports" / "partition_summary_x70norm.json")
P = 100
ALGOS = [("MOEAD", "MOEA/D"), ("NSGA2", "NSGA-II")]
DIRECT = "--direct" in sys.argv
DIFFS = []


def _half_up(value, digits):
    return Decimal(repr(float(value))).quantize(Decimal(1).scaleb(-digits), rounding=ROUND_HALF_UP)


def fmt(value, digits, where=""):
    """Value as printed in the paper (see the module docstring)."""
    two_step = _half_up(_half_up(value, 3), digits) if digits < 3 else _half_up(value, digits)
    direct = _half_up(value, digits)
    if two_step != direct:
        DIFFS.append((where, float(value), str(direct), str(two_step)))
    return str(direct if DIRECT else two_step)


def load():
    stn = pd.read_csv(MET / "stn_metrics_todas.csv", dtype={"run": str})
    part = pd.read_csv(MET / "partition_metrics_todas.csv", dtype={"run": str})
    shared = pd.read_csv(MET / "shared_alg_todas.csv", dtype={"run": str})
    rank = pd.read_csv(MET / "rank_stab_todas.csv")
    agg = lambda d: d[(d["p"] == P) & (d["front"] == "aggregated")]
    stn, part, shared = agg(stn), agg(part), agg(shared)
    keys = ["instance", "tag", "algo", "p", "run"]
    stn = stn.merge(part[keys + ["nodes_raw", "R", "D"]], on=keys, how="left")
    return stn, shared, rank.dropna(subset=["paper_no"])


def table_kappa(stn):
    """Table III: instances 7 and 10, kappa 0.5 / 1 / 2."""
    print("% Table III (tab:kappa)")
    for inst in (7, 10):
        first_inst = True
        for a, (algo, label) in enumerate(ALGOS):
            for k, kappa in enumerate(("0.5", "1.0", "2.0")):
                r = stn[(stn["paper_no"] == inst) & (stn["tag"] == f"g{kappa}norm")
                        & (stn["algo"] == algo)].iloc[0]
                w = f"Tab. III inst {inst} {label} kappa {kappa}"
                c_inst = f"{inst:<2}" if first_inst else "  "
                c_alg = f"{label:<7}" if k == 0 else " " * 7
                c_k = {"0.5": "0.5", "1.0": "1  ", "2.0": "2  "}[kappa]
                print(f"{c_inst} & {c_alg} & {c_k} & {int(r['nodes']):<5} & {fmt(r['R'], 2, w + ' Compression')} & "
                      f"{fmt(r['D'], 2, w + ' Distortion')} & {fmt(r['edge_ratio'], 2, w + ' Edge ratio')} & "
                      f"{fmt(r['par_strength'], 2, w + ' Pareto strength')} & {fmt(r['path_par'], 1, w + ' Path length')} \\\\")
                first_inst = False
            print("\\cline{2-9}" if a == 0 else "\\hline")


def table_behavior(stn, shared):
    """Table IV: kappa = 1, instances 5 to 14."""
    print("% Table IV (tab:behavior)")
    g = stn[stn["tag"] == "g1.0norm"]
    s = shared[shared["tag"] == "g1.0norm"].set_index("paper_no")
    for inst in range(5, 15):
        m = g[(g["paper_no"] == inst) & (g["algo"] == "MOEAD")].iloc[0]
        n = g[(g["paper_no"] == inst) & (g["algo"] == "NSGA2")].iloc[0]
        w = f"Tab. IV inst {inst}"
        cell = {}
        for col, digits, name in (("disc_rate", 3, "Discovery rate"), ("par_strength", 2, "Pareto strength"),
                                  ("path_par", 1, "Path length"), ("shared_vec", 2, "Vector overlap"),
                                  ("n_end", 2, "End fraction")):
            cell[col] = (fmt(m[col], digits, f"{w} {name} M"), fmt(n[col], digits, f"{w} {name} N"))
        ps, pl = cell["par_strength"], cell["path_par"]
        bold = lambda v, best: f"\\textbf{{{v}}}" if best else v
        cells = [f"{inst:<2}", fmt(s.loc[inst, "shared_alg"], 3, f"{w} Overlap"), *cell["disc_rate"],
                 bold(ps[0], float(ps[0]) > float(ps[1])), bold(ps[1], float(ps[1]) > float(ps[0])),
                 bold(pl[0], float(pl[0]) < float(pl[1])), bold(pl[1], float(pl[1]) < float(pl[0])),
                 *cell["shared_vec"], *cell["n_end"]]
        print(" & ".join(cells) + " \\\\")


def text_numbers(stn, rank):
    print("% Numbers quoted in the text")
    rep = next(r for r in json.loads(REPORT.read_text(encoding="utf-8"))
               if r["instance"] == "ns509" and r["config"] == "p100_i50")
    zero = rep["n"] - rep["positions_nonzero"]
    print(f"Inst. 4 (Shannon, x=70%): q = {rep['z']}, zero-entropy positions = {zero:,} "
          f"({fmt(100 * zero / rep['n'], 1)}%), max entropy = {fmt(rep['entropy_max'], 2)}")

    sh = stn[stn["tag"] == "x70norm"]
    oc = stn[stn["tag"] == "g1.0norm"]
    pick = lambda d, inst, algo, col: d[(d["paper_no"] == inst) & (d["algo"] == algo)][col].iloc[0]
    print(f"Inst. 4, MOEA/D edges per node: Shannon {fmt(pick(sh, 4, 'MOEAD', 'edge_ratio'), 1)}, "
          f"occupancy {fmt(pick(oc, 4, 'MOEAD', 'edge_ratio'), 1)}")
    for algo, label in ALGOS:
        f_sh = pick(sh, 4, algo, "nodes") / pick(sh, 1, algo, "nodes")
        f_oc = pick(oc, 4, algo, "nodes") / pick(oc, 2, algo, "nodes")
        print(f"{label}: Shannon STN nodes grow x{fmt(f_sh, 0)} from inst. 1 to 4; "
              f"occupancy STN grows x{fmt(f_oc, 1)} from inst. 2 to 4")
    print("Nodes, inst. 1-4 (unpartitioned / Shannon / occupancy):")
    for algo, label in ALGOS:
        raw = [int(pick(sh, i, algo, "nodes_raw")) for i in range(1, 5)]
        print(f"  {label}: {raw} / {[int(pick(sh, i, algo, 'nodes')) for i in range(1, 5)]} / "
              f"{[int(pick(oc, i, algo, 'nodes')) for i in range(1, 5)]}")

    r = rank[(rank["p"] == P) & (rank["X"] == 70)].set_index("paper_no").sort_index()
    print("Rank stability (p=100, X=70): " + ", ".join(
        f"{int(i)}: {fmt(r.loc[i, 'rank_stab'], 2, f'rank stability inst {int(i)}')}" for i in r.index))
    rest = r.loc[5:14]
    print(f"  instances 5-14: {fmt(rest['rank_stab'].min(), 2)} to {fmt(rest['rank_stab'].max(), 2)}")
    print(f"Shared retained positions: inst. 4 {fmt(100 * r.loc[4, 'retained_overlap'], 0)}%, "
          f"instances 5-14 {fmt(100 * rest['retained_overlap'].min(), 0)}% to "
          f"{fmt(100 * rest['retained_overlap'].max(), 0)}%")

    # trends with kappa on instances 5-14 (Section V-B)
    g = stn[stn["paper_no"].between(5, 14) & stn["tag"].isin(["g0.5norm", "g1.0norm", "g2.0norm"])]
    for algo, label in ALGOS:
        d = g[g["algo"] == algo].pivot(index="paper_no", columns="tag")
        mono = lambda col, sign: int(((sign * d[col]["g1.0norm"] > sign * d[col]["g0.5norm"])
                                      & (sign * d[col]["g2.0norm"] > sign * d[col]["g1.0norm"])).sum())
        weak = lambda col, sign: int(((sign * d[col]["g1.0norm"] >= sign * d[col]["g0.5norm"])
                                      & (sign * d[col]["g2.0norm"] >= sign * d[col]["g1.0norm"])).sum())
        ends = lambda col, sign: int((sign * d[col]["g2.0norm"] > sign * d[col]["g0.5norm"]).sum())
        print(f"{label}, kappa 0.5->1->2 on 10 instances (strictly monotone): nodes decrease "
              f"{mono('nodes', -1)}, compression increases {mono('R', 1)}, distortion increases "
              f"{mono('D', 1)}; Pareto strength increases {weak('par_strength', 1)} (non-decreasing; "
              f"kappa 2 > 0.5 on {ends('par_strength', 1)}), path length decreases "
              f"{weak('path_par', -1)} (non-increasing; kappa 2 < 0.5 on {ends('path_par', -1)})")

    # comparisons in Section V-C (kappa = 1, instances 5-14)
    m = oc[oc["algo"] == "MOEAD"].set_index("paper_no").sort_index().loc[5:14]
    n = oc[oc["algo"] == "NSGA2"].set_index("paper_no").sort_index().loc[5:14]
    print(f"kappa=1, instances 5-14: Discovery rate higher for NSGA-II on "
          f"{int((n['disc_rate'] > m['disc_rate']).sum())}; Pareto strength higher for MOEA/D on "
          f"{int((m['par_strength'] > n['par_strength']).sum())}; Path length shorter for MOEA/D on "
          f"{int((m['path_par'] < n['path_par']).sum())}; Vector overlap higher for NSGA-II on "
          f"{int((n['shared_vec'] > m['shared_vec']).sum())}; End fraction lower for NSGA-II on "
          f"{int((n['n_end'] < m['n_end']).sum())}")


if __name__ == "__main__":
    stn, shared, rank = load()
    table_kappa(stn)
    print()
    table_behavior(stn, shared)
    print()
    text_numbers(stn, rank)
    if DIRECT and DIFFS:
        print("\n% Cells where single (direct) rounding differs from the rule used in the paper")
        for where, value, direct, two_step in DIFFS:
            print(f"{where}: {value} -> direct {direct}, paper {two_step}")
