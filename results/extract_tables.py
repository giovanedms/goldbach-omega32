"""Extract the numerical tables of the paper from the certificate outputs.

Writes one CSV per table into this directory:
  table1_decades.csv   minimal certified margin per decade (Q={3,5}, per-cell K)
  table2_analytic.csv  analytic range: constants and certified M_lo(log Z)
  table3_margins.csv   minimal certified margins of the finite range, all scenarios
  table4_omega31.csv   exploratory float profiles for Omega(ab) <= 31
Run from the repository root:  python results/extract_tables.py
"""
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RIG = ROOT / "certificate" / "rigoroso"
OUT = ROOT / "results"
NUM = r"[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?"


def write(name, header, rows):
    with open(OUT / name, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print(f"{name}: {len(rows)} rows")


def table1():
    rows = []
    for line in (RIG / "out_cover_Q35_interval.txt").read_text().splitlines():
        m = re.match(rf"\s*z in \(({NUM}),({NUM})\]: min cert ({NUM}) on \(({NUM}),({NUM})\] "
                     rf"s=({NUM}) d=({NUM}) K=({NUM})\s+fails=(\d+)", line)
        if m:
            rows.append(m.groups())
    write("table1_decades.csv", ["z_from", "z_to", "min_certified_margin", "cell_from", "cell_to",
                                 "s", "delta", "K", "failing_cells"], rows)


def table2():
    rows, cur = [], None
    for line in (RIG / "out_analytic.txt").read_text().splitlines():
        m = re.match(rf"Q=\(([\d, ]+),?\) Y=({NUM}) K_exact\(Y\)<=({NUM}) K_an\(Y\)<=({NUM}) "
                     rf"K_used=({NUM}) cV>=({NUM})", line)
        if m:
            cur = list(m.groups())
            cur[0] = "{" + cur[0].strip().rstrip(",") + "}"
            continue
        m = re.match(rf"\s*s=({NUM}) alpha=({NUM})", line)
        if m and cur is not None:
            cur += list(m.groups())
            continue
        m = re.match(rf"\s*(CERTIFIED M\(t\) >= |NOT certified: M lower = )({NUM})", line)
        if m and cur is not None:
            rows.append(cur + [m.group(2), "certified" if m.group(1).startswith("CERT") else "not certified"])
            cur = None
    write("table2_analytic.csv", ["Q", "Z", "K_exact_le", "K_an_le", "K_used", "cV_ge", "s", "alpha",
                                  "M_lo", "status"], rows)


def table3():
    files = [("I", "out_cover_Q35_interval.txt"), ("I", "out_cover_Q3_interval.txt"),
             ("E", "out_cover_Q35_envelope.txt"), ("E", "out_cover_Q3_envelope.txt"),
             ("U", "out_cover_Q3_uniform1.146.txt"), ("U", "out_cover_Q35_uniformsup.txt")]
    rows = []
    for scen, fn in files:
        t = (RIG / fn).read_text()
        q = re.search(r"# Q=\(([\d, ]+),?\) XMAX=(\S+)", t)
        cells = re.search(r"cells=(\d+)", t)
        fail = re.search(r"failing\(cert<=0\)=(\d+)", t)
        mn = re.search(rf"min certified margin = ({NUM}) on \(({NUM}),({NUM})\]", t)
        rows.append([scen, "{" + q.group(1).strip().rstrip(",") + "}", q.group(2),
                     cells.group(1) if cells else "", fail.group(1) if fail else "",
                     mn.group(1) if mn else "", mn.group(2) if mn else "", mn.group(3) if mn else "", fn])
    write("table3_margins.csv", ["scenario", "Q", "up_to", "cells", "failing", "min_margin",
                                 "cell_from", "cell_to", "source"], rows)


def table4():
    rows = []
    for line in (ROOT / "certificate" / "profiles" / "out_profiles.txt").read_text().splitlines():
        m = re.match(rf"(odd|even) M=(\d+)(?: \w+)?: Q=\(([\d, ]+),?\) c=({NUM}) r=({NUM}) z0=({NUM})\s+"
                     rf"min logratio=({NUM}) at z=({NUM})", line)
        if m:
            g = list(m.groups())
            g[2] = "{" + g[2].strip().rstrip(",") + "}"
            rows.append(g)
    write("table4_omega31.csv", ["parity", "bound_M", "Q", "c", "r", "z0", "min_log_ratio", "at_z"], rows)


if __name__ == "__main__":
    table1(); table2(); table3(); table4()
