#!/usr/bin/env python3
"""
IRC エネルギー曲線をプロット
"""
import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def extract_irc_energies_from_xyz(xyz_file):
    """IRC 軌跡 XYZ ファイルからエネルギーを抽出"""
    energies = []
    steps = []

    with open(xyz_file, "r") as f:
        lines = f.readlines()

    step = 0
    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # 原子数行（各構造の最初の行）
        if line.isdigit():
            n_atoms = int(line)
            # 次の行にエネルギー情報がある
            if i + 1 < len(lines):
                comment_line = lines[i + 1]
                # "Coordinates from ORCA-job ... E -399.088864922217" の形式
                energy_match = re.search(r"E\s+(-?\d+\.\d+)", comment_line)
                if energy_match:
                    energy = float(energy_match.group(1))
                    energies.append(energy)
                    steps.append(step)
                    step += 1

            # 次の構造へスキップ (原子数 + コメント行 + 座標行)
            i += n_atoms + 2
        else:
            i += 1

    if len(energies) == 0:
        return None, None

    return np.array(steps), np.array(energies)


def extract_irc_energies(log_file):
    """IRC ログファイルからエネルギーを抽出（後方互換性のため残す）"""
    # 軌跡ファイルから読み込みを試みる
    log_path = Path(log_file)

    if "forward" in log_path.name:
        xyz_file = log_path.parent / "irc_forward_IRC_F_trj.xyz"
    elif "reverse" in log_path.name:
        xyz_file = log_path.parent / "irc_reverse_IRC_B_trj.xyz"
    else:
        return None, None

    if xyz_file.exists():
        return extract_irc_energies_from_xyz(xyz_file)

    return None, None


def main():
    irc_dir = Path("/home/yamamoto/Work/chemistry/orca-work/MEA_water/irc")

    forward_log = irc_dir / "irc_forward.log"
    reverse_log = irc_dir / "irc_reverse.log"

    # Forward IRC
    f_steps, f_energies = extract_irc_energies(forward_log)

    # Reverse IRC
    r_steps, r_energies = extract_irc_energies(reverse_log)

    # プロット
    plt.figure(figsize=(12, 7))

    # Hartree to kcal/mol
    hartree_to_kcal = 627.509

    if f_energies is not None:
        # TS を 0 とした相対エネルギー
        ts_energy = f_energies[0]
        f_rel = (f_energies - ts_energy) * hartree_to_kcal
        plt.plot(
            f_steps,
            f_rel,
            "o-",
            color="blue",
            linewidth=2,
            markersize=6,
            label="Forward IRC",
        )

    if r_energies is not None:
        ts_energy = r_energies[0] if f_energies is None else f_energies[0]
        r_rel = (r_energies - ts_energy) * hartree_to_kcal
        # Reverse は負の座標にプロット
        r_steps_neg = -r_steps if len(r_steps) > 0 else r_steps
        plt.plot(
            r_steps_neg,
            r_rel,
            "o-",
            color="red",
            linewidth=2,
            markersize=6,
            label="Reverse IRC",
        )

    # TS 位置にマーカー
    plt.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
    plt.axvline(x=0, color="gray", linestyle="--", alpha=0.5)
    plt.plot(0, 0, "k*", markersize=20, label="TS")

    plt.xlabel("IRC Step", fontsize=14, fontweight="bold")
    plt.ylabel("Relative Energy (kcal/mol)", fontsize=14, fontweight="bold")
    plt.title(
        "Intrinsic Reaction Coordinate (IRC) Energy Profile",
        fontsize=16,
        fontweight="bold",
    )
    plt.legend(fontsize=12, loc="best")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # 保存
    output_file = irc_dir / "irc_energy_profile.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"プロット保存: {output_file}")

    # データも保存
    data_file = irc_dir / "irc_energies.txt"
    with open(data_file, "w") as f:
        f.write("# IRC Energy Data\n")
        f.write("# Step  Energy(Eh)  Relative(kcal/mol)\n")

        if r_energies is not None:
            ts_e = r_energies[0] if f_energies is None else f_energies[0]
            for step, energy in zip(r_steps_neg, r_energies):
                rel_e = (energy - ts_e) * hartree_to_kcal
                f.write(f"{step:6.0f}  {energy:16.10f}  {rel_e:12.4f}\n")

        if f_energies is not None:
            ts_e = f_energies[0]
            for step, energy in zip(f_steps, f_energies):
                rel_e = (energy - ts_e) * hartree_to_kcal
                f.write(f"{step:6.0f}  {energy:16.10f}  {rel_e:12.4f}\n")

    print(f"データ保存: {data_file}")

    # 統計情報
    print("\n=== IRC エネルギー統計 ===")
    if f_energies is not None:
        print(f"Forward IRC: {len(f_energies)} steps")
        print(f"  開始エネルギー: {f_energies[0]:.8f} Eh")
        print(f"  終了エネルギー: {f_energies[-1]:.8f} Eh")
        print(
            f"  相対エネルギー: {(f_energies[-1] - f_energies[0]) * hartree_to_kcal:.2f} kcal/mol"
        )

    if r_energies is not None:
        print(f"\nReverse IRC: {len(r_energies)} steps")
        print(f"  開始エネルギー: {r_energies[0]:.8f} Eh")
        print(f"  終了エネルギー: {r_energies[-1]:.8f} Eh")
        ts_e = r_energies[0] if f_energies is None else f_energies[0]
        print(
            f"  相対エネルギー: {(r_energies[-1] - ts_e) * hartree_to_kcal:.2f} kcal/mol"
        )


if __name__ == "__main__":
    main()
