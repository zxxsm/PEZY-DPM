#!/usr/bin/env python3
"""Compare a DPM central-axis dose output with a reference CSV."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Iterable


def parse_number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def read_reference(path: Path, x_cm: float, y_cm: float) -> dict[float, float]:
    rows: dict[float, float] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            x = parse_number(row["x_cm"])
            y = parse_number(row["y_cm"])
            if not (math.isclose(x, x_cm, abs_tol=1.0e-6) and math.isclose(y, y_cm, abs_tol=1.0e-6)):
                continue
            z = round(parse_number(row["z_cm"]), 6)
            dose = parse_number(row["dose_MeV_per_g"])
            rows[z] = dose
    return rows


def read_candidate_csv(path: Path, x_cm: float, y_cm: float) -> dict[float, float]:
    return read_reference(path, x_cm, y_cm)


def numeric_rows(lines: Iterable[str]) -> Iterable[tuple[float, float, float, float]]:
    for line in lines:
        fields = line.split()
        if len(fields) != 5:
            continue
        try:
            x, y, z, dose, _uncertainty = (parse_number(field) for field in fields)
        except ValueError:
            continue
        yield x, y, z, dose


def read_dpm_output(path: Path, x_cm: float, y_cm: float) -> dict[float, float]:
    rows: dict[float, float] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for x, y, z, dose in numeric_rows(handle):
            if not (math.isclose(x, x_cm, abs_tol=1.0e-6) and math.isclose(y, y_cm, abs_tol=1.0e-6)):
                continue
            rows[round(z, 6)] = dose
    return rows


def write_candidate(path: Path, rows: dict[float, float], x_cm: float, y_cm: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["x_cm", "y_cm", "z_cm", "dose_MeV_per_g"])
        for z, dose in sorted(rows.items()):
            writer.writerow([x_cm, y_cm, z, dose])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--mean-threshold", type=float, required=True, help="Allowed mean absolute difference in %%Dmax")
    parser.add_argument("--expected-points", type=int, default=150)
    parser.add_argument("--x", type=float, default=15.25)
    parser.add_argument("--y", type=float, default=15.25)
    parser.add_argument("--write-candidate-csv", type=Path)
    args = parser.parse_args()

    reference = read_reference(args.reference, args.x, args.y)
    if args.candidate.suffix.lower() == ".csv":
        candidate = read_candidate_csv(args.candidate, args.x, args.y)
    else:
        candidate = read_dpm_output(args.candidate, args.x, args.y)

    if args.write_candidate_csv:
        write_candidate(args.write_candidate_csv, candidate, args.x, args.y)

    common = sorted(set(reference) & set(candidate))
    finite = all(math.isfinite(reference[z]) and math.isfinite(candidate[z]) for z in common)
    nonnegative = all(reference[z] >= 0.0 and candidate[z] >= 0.0 for z in common)

    if not reference or not common:
        print("FAIL: no matching central-axis dose points")
        return 1

    dmax = max(reference.values())
    differences = [(candidate[z] - reference[z]) / dmax * 100.0 for z in common]
    mean_abs = sum(abs(value) for value in differences) / len(differences)
    max_abs = max(abs(value) for value in differences)
    rms = math.sqrt(sum(value * value for value in differences) / len(differences))

    print(f"matched_points={len(common)}")
    print(f"mean_abs_diff_percent_Dmax={mean_abs:.6f}")
    print(f"max_abs_diff_percent_Dmax={max_abs:.6f}")
    print(f"rms_diff_percent_Dmax={rms:.6f}")
    print(f"mean_threshold_percent_Dmax={args.mean_threshold:.6f}")

    passed = (
        len(common) == args.expected_points
        and len(reference) == args.expected_points
        and len(candidate) == args.expected_points
        and finite
        and nonnegative
        and mean_abs <= args.mean_threshold
    )
    print("PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
