#!/usr/bin/env python3
"""
Extract normal modes from ORCA log and displace geometry along a specific mode
"""
import sys
import re
import argparse
import numpy as np


def parse_geometry_from_xyz(xyz_file):
    """Parse XYZ file and return atoms and coordinates in Angstrom"""
    with open(xyz_file, "r") as f:
        lines = f.readlines()

    n_atoms = int(lines[0].strip())
    atoms = []
    coords = []

    for i in range(2, 2 + n_atoms):
        parts = lines[i].split()
        atoms.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    return atoms, np.array(coords)


def parse_normal_modes_from_log(log_file, n_atoms):
    """
    Parse normal modes from ORCA log file
    Returns: numpy array of shape (3*n_atoms, n_modes)
    """
    with open(log_file, "r") as f:
        lines = f.readlines()

    # Find NORMAL MODES section
    start_idx = None
    for i, line in enumerate(lines):
        if "NORMAL MODES" in line and "------------" in lines[i + 1]:
            start_idx = i
            break

    if start_idx is None:
        raise ValueError("NORMAL MODES section not found in log file")

    # Skip header lines
    i = start_idx + 1
    while i < len(lines) and lines[i].strip().startswith(("-", "These modes", "M(i")):
        i += 1

    # Parse mode data
    n_coords = 3 * n_atoms
    modes_data = {}  # mode_idx -> list of values

    while i < len(lines):
        line = lines[i]

        # Check for mode header line (column headers)
        if re.match(r"\s+\d+\s+\d+", line):
            mode_indices = [int(x) for x in line.split()]
            # Initialize modes if not exists
            for mode_idx in mode_indices:
                if mode_idx not in modes_data:
                    modes_data[mode_idx] = []
            i += 1

            # Read data rows
            for row_idx in range(n_coords):
                if i >= len(lines):
                    break
                data_line = lines[i]

                # Parse row: "      0       0.132967  -0.003352   0.006607 ..."
                match = re.match(r"\s+(\d+)\s+(.*)", data_line)
                if match:
                    row_num = int(match.group(1))
                    values = [float(x) for x in match.group(2).split()]

                    # Append values to corresponding modes
                    for j, mode_idx in enumerate(mode_indices):
                        if j < len(values):
                            modes_data[mode_idx].append(values[j])
                    i += 1
                else:
                    break
        else:
            i += 1
            # Stop if we hit another section
            if "IR SPECTRUM" in line or "RAMAN" in line:
                break

    # Convert to numpy array
    n_modes = max(modes_data.keys()) + 1 if modes_data else 0
    modes_array = np.zeros((n_coords, n_modes))

    for mode_idx, values in modes_data.items():
        if len(values) == n_coords:
            modes_array[:, mode_idx] = values

    return modes_array


def displace_geometry(coords, mode_vector, scaling=0.25):
    """
    Displace geometry along mode vector

    Args:
        coords: numpy array of shape (n_atoms, 3) in Angstrom
        mode_vector: numpy array of shape (3*n_atoms,) - normalized displacement
        scaling: displacement scaling factor in Bohr (will convert to Angstrom)

    Returns:
        new_coords: numpy array of shape (n_atoms, 3) in Angstrom
    """
    n_atoms = coords.shape[0]

    # Reshape mode vector to (n_atoms, 3)
    mode_3d = mode_vector.reshape(n_atoms, 3)

    # Convert scaling from Bohr to Angstrom (1 Bohr = 0.529177 Angstrom)
    scaling_ang = scaling * 0.529177

    # Apply displacement
    new_coords = coords + scaling_ang * mode_3d

    return new_coords


def write_xyz(filename, atoms, coords):
    """Write XYZ file"""
    n_atoms = len(atoms)
    with open(filename, "w") as f:
        f.write(f"{n_atoms}\n")
        f.write(f"Displaced geometry\n")
        for atom, coord in zip(atoms, coords):
            f.write(
                f"{atom:2s}  {coord[0]:14.8f}  {coord[1]:14.8f}  {coord[2]:14.8f}\n"
            )


def main():
    parser = argparse.ArgumentParser(description="Displace geometry along normal mode")
    parser.add_argument("--xyz", required=True, help="Input XYZ file")
    parser.add_argument("--log", required=True, help="ORCA log file with normal modes")
    parser.add_argument("--mode", type=int, required=True, help="Mode index to follow")
    parser.add_argument(
        "--scaling", type=float, default=0.25, help="Displacement scaling (Bohr)"
    )
    parser.add_argument("--output", required=True, help="Output XYZ file")

    args = parser.parse_args()

    # Parse geometry
    atoms, coords = parse_geometry_from_xyz(args.xyz)
    n_atoms = len(atoms)
    print(f"Read {n_atoms} atoms from {args.xyz}")

    # Parse normal modes
    modes = parse_normal_modes_from_log(args.log, n_atoms)
    print(f"Parsed normal modes: shape = {modes.shape}")

    # Check if requested mode exists
    if args.mode >= modes.shape[1]:
        print(
            f"Error: Mode {args.mode} not found (max mode index: {modes.shape[1]-1})",
            file=sys.stderr,
        )
        return 1

    # Get mode vector
    mode_vector = modes[:, args.mode]
    print(f"Mode {args.mode}: ||v|| = {np.linalg.norm(mode_vector):.6f}")

    # Displace geometry
    new_coords = displace_geometry(coords, mode_vector, args.scaling)

    # Calculate displacement
    displacement = np.linalg.norm(new_coords - coords)
    print(f"Total displacement: {displacement:.6f} Angstrom")

    # Write output
    write_xyz(args.output, atoms, new_coords)
    print(f"Wrote displaced geometry to {args.output}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
