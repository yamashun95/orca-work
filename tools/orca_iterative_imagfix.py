#!/usr/bin/env python3
"""Iteratively remove imaginary frequencies by displacing along the imaginary normal mode.

This script is designed for ORCA workflows where you have a directory containing:
- orca.inp (uses '* xyzfile ... orca.xyz')
- orca.xyz
- orca.hess
- orca.log

It will:
1) Detect the most negative (imaginary) vibrational frequency from orca.hess.
2) Extract that mode's mass-weighted Cartesian eigenvector from orca.hess.
3) Un-massweight it, scale to a small RMS displacement (in Å), and generate displaced geometries.
4) Run ORCA in new iteration subdirectories until imaginary modes are gone.

Notes
- ORCA's "NORMAL MODES" in the log/hess are mass-weighted by 1/sqrt(m).
  This script converts back to Cartesian displacements via sqrt(m).
- By default it tries BOTH + and - displacement and continues from the better one.

Usage example
  export ORCA_CMD="$HOME/orca_6_1_1/orca"
  python3 tools/orca_iterative_imagfix.py --start-dir afify/b3lyp/wo_correction/Step4_Validate_TS/Step4B_Validate_Forwards_and_Backwards/Backwards \
      --step-rms 0.02 --max-iter 10

"""

from __future__ import annotations

import argparse
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

BOHR_TO_ANG = 0.529177210903


@dataclass(frozen=True)
class ImagResult:
    mode_index: int
    freq_cm1: float


@dataclass(frozen=True)
class CandidateResult:
    path: Path
    imag_modes: List[ImagResult]
    ok: bool
    error: Optional[str] = None


def read_xyz(path: Path) -> Tuple[List[str], List[List[float]], str]:
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"Invalid XYZ (too short): {path}")
    try:
        n = int(lines[0].strip())
    except Exception as e:
        raise ValueError(f"Invalid XYZ first line (natoms): {path}") from e
    comment = lines[1]
    body = lines[2:]
    if len(body) < n:
        raise ValueError(
            f"Invalid XYZ (expected {n} atom lines, got {len(body)}): {path}"
        )
    symbols: List[str] = []
    coords: List[List[float]] = []
    for i in range(n):
        parts = body[i].split()
        if len(parts) < 4:
            raise ValueError(f"Invalid XYZ atom line {i+3} in {path}: {body[i]}")
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return symbols, coords, comment


def write_xyz(
    path: Path, symbols: List[str], coords: List[List[float]], comment: str
) -> None:
    if len(symbols) != len(coords):
        raise ValueError("symbols/coords length mismatch")
    out = [str(len(symbols)), comment]
    for s, (x, y, z) in zip(symbols, coords):
        out.append(f"{s:>2s} {x: 16.10f} {y: 16.10f} {z: 16.10f}")
    path.write_text("\n".join(out) + "\n")


def parse_hess_vibrational_frequencies(hess_path: Path) -> List[float]:
    freqs: List[float] = []
    with hess_path.open() as f:
        for line in f:
            if line.strip() == "$vibrational_frequencies":
                n = int(next(f).strip())
                for _ in range(n):
                    parts = next(f).split()
                    if len(parts) < 2:
                        continue
                    freqs.append(float(parts[1]))
                break
    if not freqs:
        raise ValueError(f"Could not find $vibrational_frequencies in {hess_path}")
    return freqs


def parse_hess_atoms(hess_path: Path) -> Tuple[List[str], List[float]]:
    """Return (labels, masses_amu) from $atoms."""
    labels: List[str] = []
    masses: List[float] = []
    with hess_path.open() as f:
        for line in f:
            if line.strip() == "$atoms":
                n = int(next(f).strip())
                for _ in range(n):
                    parts = next(f).split()
                    if len(parts) < 5:
                        raise ValueError(
                            f"Malformed $atoms line in {hess_path}: {' '.join(parts)}"
                        )
                    labels.append(parts[0])
                    masses.append(float(parts[1]))
                break
    if not labels:
        raise ValueError(f"Could not find $atoms in {hess_path}")
    return labels, masses


def select_imaginary_modes(
    freqs_cm1: List[float], imag_cutoff: float
) -> List[ImagResult]:
    # imag_cutoff is positive, treat freq < -imag_cutoff as imaginary
    out: List[ImagResult] = []
    for i, w in enumerate(freqs_cm1):
        if w < -abs(imag_cutoff):
            out.append(ImagResult(mode_index=i, freq_cm1=w))
    out.sort(key=lambda x: x.freq_cm1)  # most negative first
    return out


def parse_hess_normal_mode_column(hess_path: Path, mode_index: int) -> List[float]:
    """Extract a single column (mode_index) from the $normal_modes matrix.

    Returns a vector of length 3N with mass-weighted Cartesian components.
    """
    target = int(mode_index)
    in_modes = False
    header_cols: List[int] = []
    vec: dict[int, float] = {}

    int_re = re.compile(r"^\s*\d+\s*$")

    with hess_path.open() as f:
        for line in f:
            if not in_modes:
                if line.strip() == "$normal_modes":
                    in_modes = True
                    # Next line: "42 42"
                    dims = next(f).split()
                    if len(dims) < 2:
                        raise ValueError(
                            f"Malformed $normal_modes header in {hess_path}"
                        )
                    nrow = int(dims[0])
                    ncol = int(dims[1])
                    if target < 0 or target >= ncol:
                        raise ValueError(
                            f"mode_index {target} out of range (0..{ncol-1})"
                        )
                    # After this, there are repeated blocks: header row with column indices,
                    # followed by nrow lines each containing row index + values.
                    continue
                else:
                    continue

            # In $normal_modes
            stripped = line.rstrip("\n")
            if not stripped.strip():
                continue
            if (
                stripped.lstrip().startswith("$")
                and stripped.strip() != "$normal_modes"
            ):
                break

            # Header line has only integers (column indices) and spaces.
            # Example: "                    0                  1 ..."
            parts = stripped.split()
            if parts and all(int_re.match(p) for p in parts):
                header_cols = [int(p) for p in parts]
                continue

            # Data line: row_index + len(header_cols) floats
            parts = stripped.split()
            if not parts:
                continue
            try:
                row = int(parts[0])
            except ValueError:
                continue
            if not header_cols:
                continue
            values = parts[1:]
            if len(values) < len(header_cols):
                # Some ORCA versions may wrap; not handled.
                raise ValueError(
                    f"Unexpected wrapped/short $normal_modes line for row {row} in {hess_path}"
                )
            if target in header_cols:
                j = header_cols.index(target)
                vec[row] = float(values[j])

    if not vec:
        raise ValueError(
            f"Failed to extract mode {target} from $normal_modes in {hess_path}"
        )

    # Ensure contiguous 0..3N-1
    n = max(vec.keys()) + 1
    out = [0.0] * n
    for k, v in vec.items():
        out[k] = v
    return out


def mode_vector_to_cart_displacements_ang(
    mode_vec_massweighted: List[float],
    masses_amu: List[float],
    step_rms_ang: float,
) -> List[List[float]]:
    """Convert mass-weighted mode vector to Cartesian displacements in Å.

    Assumes the mode vector is in (bohr / sqrt(amu)) and converts to bohr
    by multiplying sqrt(mass_amu), then to Å.

    The result is scaled so that RMS(|dr_atom|) == step_rms_ang.
    """
    n_atoms = len(masses_amu)
    if len(mode_vec_massweighted) < 3 * n_atoms:
        raise ValueError(
            f"Mode vector too short: got {len(mode_vec_massweighted)}, need {3*n_atoms}"
        )

    # Un-massweight: bohr displacement per coordinate
    disp_bohr: List[float] = []
    for atom_i in range(n_atoms):
        m = masses_amu[atom_i]
        if m <= 0:
            raise ValueError(f"Non-positive mass for atom {atom_i}: {m}")
        sm = math.sqrt(m)
        for k in range(3):
            disp_bohr.append(mode_vec_massweighted[3 * atom_i + k] * sm)

    # Convert to Å
    disp_ang = [d * BOHR_TO_ANG for d in disp_bohr]

    # Compute RMS per atom
    per_atom_sq: List[float] = []
    for atom_i in range(n_atoms):
        dx, dy, dz = disp_ang[3 * atom_i : 3 * atom_i + 3]
        per_atom_sq.append(dx * dx + dy * dy + dz * dz)
    rms = math.sqrt(sum(per_atom_sq) / n_atoms) if n_atoms else 0.0
    if rms == 0.0:
        raise ValueError("Mode displacement RMS is zero; cannot scale")

    scale = step_rms_ang / rms

    out: List[List[float]] = []
    for atom_i in range(n_atoms):
        dx, dy, dz = disp_ang[3 * atom_i : 3 * atom_i + 3]
        out.append([dx * scale, dy * scale, dz * scale])
    return out


def apply_displacement(
    coords: List[List[float]], disp: List[List[float]], sign: float
) -> List[List[float]]:
    if len(coords) != len(disp):
        raise ValueError("coords/disp length mismatch")
    out: List[List[float]] = []
    for (x, y, z), (dx, dy, dz) in zip(coords, disp):
        out.append([x + sign * dx, y + sign * dy, z + sign * dz])
    return out


def ensure_orca_xyzfile_input(inp_path: Path) -> None:
    txt = inp_path.read_text()
    if "* xyzfile" not in txt:
        raise ValueError(
            f"Expected '* xyzfile ...' in {inp_path}. This script assumes XYZFILE input."
        )


def run_orca(
    workdir: Path, orca_cmd: str, inp_name: str, tee_name: str, quiet: bool
) -> None:
    cmd = shlex.split(orca_cmd) + [inp_name]
    log_path = workdir / tee_name

    with log_path.open("w") as logf:
        proc = subprocess.Popen(
            cmd,
            cwd=str(workdir),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            logf.write(line)
            if not quiet:
                sys.stdout.write(line)
        ret = proc.wait()

    if ret != 0:
        raise RuntimeError(f"ORCA failed in {workdir} (exit={ret}). See {log_path}")


def analyze_dir_for_imag(workdir: Path, imag_cutoff: float) -> List[ImagResult]:
    hess_path = workdir / "orca.hess"
    if not hess_path.exists():
        raise FileNotFoundError(f"Missing {hess_path}")
    freqs = parse_hess_vibrational_frequencies(hess_path)
    return select_imaginary_modes(freqs, imag_cutoff=imag_cutoff)


def copy_inputs_to_dir(src_dir: Path, dst_dir: Path, inp_name: str) -> None:
    dst_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src_dir / inp_name, dst_dir / inp_name)


def main() -> int:
    p = argparse.ArgumentParser(
        description="Iteratively displace along imaginary mode and re-optimize until no imaginary frequencies remain."
    )
    p.add_argument(
        "--start-dir",
        type=Path,
        required=True,
        help="Directory containing orca.inp/orca.xyz/orca.hess/orca.log",
    )
    p.add_argument(
        "--out-root",
        type=str,
        default="imagfix_runs",
        help="Output directory (created under start-dir) to store iteration subdirectories (default: imagfix_runs)",
    )
    p.add_argument(
        "--orca-cmd",
        type=str,
        default=os.environ.get("ORCA_CMD", "orca"),
        help="ORCA executable command (default: $ORCA_CMD or 'orca')",
    )
    p.add_argument(
        "--inp",
        type=str,
        default="orca.inp",
        help="Input filename (default: orca.inp)",
    )
    p.add_argument(
        "--step-rms",
        type=float,
        default=0.02,
        help="RMS displacement per atom in Å (default: 0.02)",
    )
    p.add_argument(
        "--imag-cutoff",
        type=float,
        default=1.0,
        help="Treat frequencies < -imag_cutoff (cm^-1) as imaginary (default: 1.0)",
    )
    p.add_argument(
        "--max-iter",
        type=int,
        default=10,
        help="Maximum number of displacement+reopt cycles (default: 10)",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Only create displaced geometries and dirs; do not run ORCA",
    )
    p.add_argument(
        "--quiet",
        action="store_true",
        help="Do not stream ORCA output to stdout (still saved to orca.log)",
    )
    p.add_argument(
        "--one-sign",
        choices=["+", "-"],
        default=None,
        help="Only try one displacement sign ('+' or '-'). Default: try both and keep the better.",
    )

    args = p.parse_args()

    start_dir = args.start_dir.resolve()
    cur_dir = start_dir
    if not cur_dir.exists():
        raise FileNotFoundError(cur_dir)

    run_root = start_dir / args.out_root
    run_root.mkdir(parents=True, exist_ok=True)

    inp_path = cur_dir / args.inp
    if not inp_path.exists():
        raise FileNotFoundError(inp_path)
    ensure_orca_xyzfile_input(inp_path)

    print(f"Start dir: {cur_dir}")
    print(f"Out root : {run_root}")
    print(f"ORCA cmd : {args.orca_cmd}")
    print(f"step_rms : {args.step_rms:.4f} Å")
    print(f"imag_cut : {args.imag_cutoff:.3f} cm^-1")

    # If already ok, exit.
    imag0 = analyze_dir_for_imag(cur_dir, imag_cutoff=args.imag_cutoff)
    if not imag0:
        print("No imaginary modes detected. Nothing to do.")
        return 0
    print(f"Initial imaginary modes: {[(m.mode_index, m.freq_cm1) for m in imag0]}")

    for it in range(1, args.max_iter + 1):
        hess_path = cur_dir / "orca.hess"
        xyz_path = cur_dir / "orca.xyz"
        if not hess_path.exists() or not xyz_path.exists():
            raise FileNotFoundError(f"Missing orca.hess/orca.xyz in {cur_dir}")

        freqs = parse_hess_vibrational_frequencies(hess_path)
        imag = select_imaginary_modes(freqs, imag_cutoff=args.imag_cutoff)
        if not imag:
            print(f"Iteration {it-1}: no imaginary modes. Done.")
            return 0

        target = imag[0]  # most negative
        labels, masses = parse_hess_atoms(hess_path)
        symbols_xyz, coords_xyz, comment_xyz = read_xyz(xyz_path)
        if symbols_xyz != labels:
            # Not necessarily fatal, but usually should match ordering.
            print(
                "WARNING: Atom labels differ between orca.xyz and orca.hess $atoms. Proceeding."
            )

        mode_vec = parse_hess_normal_mode_column(hess_path, target.mode_index)
        disp = mode_vector_to_cart_displacements_ang(
            mode_vec_massweighted=mode_vec,
            masses_amu=masses,
            step_rms_ang=args.step_rms,
        )

        print(
            f"\nIter {it}: mode {target.mode_index} = {target.freq_cm1:.2f} cm^-1; generating displaced geometries"
        )

        signs: List[Tuple[str, float]]
        if args.one_sign == "+":
            signs = [("plus", +1.0)]
        elif args.one_sign == "-":
            signs = [("minus", -1.0)]
        else:
            signs = [("plus", +1.0), ("minus", -1.0)]

        candidates: List[CandidateResult] = []
        for sign_name, sign_val in signs:
            out_dir = run_root / f"iter_{it:03d}_{sign_name}"
            try:
                if out_dir.exists():
                    raise FileExistsError(f"Refusing to overwrite existing {out_dir}")
                copy_inputs_to_dir(cur_dir, out_dir, inp_name=args.inp)

                displaced_coords = apply_displacement(coords_xyz, disp, sign=sign_val)
                write_xyz(
                    out_dir / "orca.xyz",
                    symbols_xyz,
                    displaced_coords,
                    comment=f"imagfix iter {it} {sign_name} from {cur_dir.name} mode {target.mode_index} {target.freq_cm1:.2f} cm^-1",
                )

                if args.dry_run:
                    candidates.append(
                        CandidateResult(path=out_dir, imag_modes=[], ok=True)
                    )
                    print(f"  - created (dry-run): {out_dir}")
                    continue

                print(f"  - running ORCA in: {out_dir}")
                run_orca(
                    workdir=out_dir,
                    orca_cmd=args.orca_cmd,
                    inp_name=args.inp,
                    tee_name="orca.log",
                    quiet=args.quiet,
                )

                im = analyze_dir_for_imag(out_dir, imag_cutoff=args.imag_cutoff)
                candidates.append(CandidateResult(path=out_dir, imag_modes=im, ok=True))
                print(
                    f"    result: n_imag={len(im)}; most_negative={(im[0].freq_cm1 if im else None)}"
                )
            except Exception as e:
                candidates.append(
                    CandidateResult(path=out_dir, imag_modes=[], ok=False, error=str(e))
                )
                print(f"    FAILED: {out_dir}: {e}")

        if args.dry_run:
            print("Dry-run complete. No ORCA jobs executed.")
            return 0

        ok_candidates = [c for c in candidates if c.ok]
        if not ok_candidates:
            print("All candidates failed; aborting.")
            return 2

        # Select best: fewer imaginary modes, then least negative magnitude (closest to 0)
        def score(c: CandidateResult) -> Tuple[int, float]:
            if not c.imag_modes:
                return (0, 0.0)
            most_neg = c.imag_modes[0].freq_cm1  # negative
            return (len(c.imag_modes), abs(most_neg))

        best = min(ok_candidates, key=score)
        cur_dir = best.path

        if not best.imag_modes:
            print(f"\nSuccess: no imaginary modes in {cur_dir}")
            return 0

        print(
            f"\nContinuing from best candidate: {cur_dir} (n_imag={len(best.imag_modes)}, most_negative={best.imag_modes[0].freq_cm1:.2f})"
        )

    print(
        f"Reached max iterations ({args.max_iter}) without eliminating imaginary modes."
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
