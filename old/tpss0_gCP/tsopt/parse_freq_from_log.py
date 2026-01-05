#!/usr/bin/env python3
"""
Parse vibrational frequencies from ORCA log file
"""
import sys
import re
import argparse


def parse_frequencies_from_log(log_file, cutoff=-10.0):
    """
    Parse vibrational frequencies from ORCA log file

    Args:
        log_file: Path to ORCA log file
        cutoff: Cutoff for negative frequencies (cm^-1)

    Returns:
        tuple: (frequencies_list, negative_count, negative_indices)
    """
    frequencies = []

    try:
        with open(log_file, "r") as f:
            lines = f.readlines()

        # Find VIBRATIONAL FREQUENCIES section
        in_freq_section = False
        for i, line in enumerate(lines):
            if "VIBRATIONAL FREQUENCIES" in line:
                in_freq_section = True
                continue

            if in_freq_section:
                # Look for frequency lines: "     0:       0.00 cm**-1" or "     6:    -186.89 cm**-1"
                match = re.match(r"\s+(\d+):\s+(-?\d+\.\d+)\s+cm\*\*-1", line)
                if match:
                    mode_idx = int(match.group(1))
                    freq_value = float(match.group(2))
                    frequencies.append(freq_value)
                    continue
                # End of frequency section - look for NORMAL MODES
                if "NORMAL MODES" in line:
                    break

        if not frequencies:
            print(f"Warning: No frequencies found in {log_file}", file=sys.stderr)
            return [], 0, [], 0.0

        # Count negative frequencies below cutoff
        negative_indices = [i for i, freq in enumerate(frequencies) if freq < cutoff]
        neg_count = len(negative_indices)

        # Find minimum frequency (excluding first 6 which should be ~0 for translations/rotations)
        real_freqs = frequencies[6:] if len(frequencies) > 6 else frequencies
        min_freq = min(real_freqs) if real_freqs else 0.0

        return frequencies, neg_count, negative_indices, min_freq

    except FileNotFoundError:
        print(f"Error: Log file not found: {log_file}", file=sys.stderr)
        return [], 0, [], 0.0
    except Exception as e:
        print(f"Error parsing log file: {e}", file=sys.stderr)
        return [], 0, [], 0.0


def main():
    parser = argparse.ArgumentParser(description="Parse frequencies from ORCA log file")
    parser.add_argument("--log", required=True, help="ORCA log file")
    parser.add_argument(
        "--cutoff",
        type=float,
        default=-10.0,
        help="Cutoff for negative frequencies (default: -10.0 cm^-1)",
    )
    parser.add_argument("--output", help="Output file for results")
    parser.add_argument(
        "--target", type=int, default=1, help="Target number of negative frequencies"
    )

    args = parser.parse_args()

    frequencies, neg_count, negative_indices, min_freq = parse_frequencies_from_log(
        args.log, args.cutoff
    )

    if not frequencies:
        result = f"Frequencies (cm^-1): min=0.00, neg_count=0\nNegative indices: []\nSTATUS: ERROR\n"
    else:
        result = f"Frequencies (cm^-1): min={min_freq:.2f}, neg_count={neg_count}\n"
        result += f"Negative indices: {negative_indices}\n"

        # Determine status
        if neg_count > args.target:
            result += "STATUS: TOO_MANY\n"
        elif neg_count < args.target:
            result += "STATUS: TOO_FEW\n"
        else:
            result += "STATUS: CONVERGED\n"

        # Show all frequencies
        result += f"\nAll frequencies:\n"
        for i, freq in enumerate(frequencies):
            if freq < args.cutoff:
                result += f"  {i}: {freq:10.2f} cm^-1  ***imaginary***\n"
            elif i < 6:
                result += f"  {i}: {freq:10.2f} cm^-1  (trans/rot)\n"
            else:
                result += f"  {i}: {freq:10.2f} cm^-1\n"

    # Write to output file if specified
    if args.output:
        with open(args.output, "w") as f:
            f.write(result)

    # Always print to stdout for script capture
    print(result, end="")

    return 0 if neg_count == args.target else 1


if __name__ == "__main__":
    sys.exit(main())
