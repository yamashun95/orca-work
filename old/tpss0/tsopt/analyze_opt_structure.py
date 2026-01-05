import numpy as np

# 初期構造 (scan.013)
scan013_coords = np.array(
    [
        [1.595097, -0.234332, -0.700429],  # O
        [-2.039963, -0.126794, -0.102604],  # N
        [-0.769206, -0.150447, -0.817740],  # C (MEA)
        [0.428577, -0.224510, 0.109090],  # C (MEA)
        [-0.706053, 0.748128, -1.428487],  # H
        [-0.771404, -1.008449, -1.487156],  # H
        [0.371578, -1.131385, 0.717200],  # H
        [0.430662, 0.638266, 0.780822],  # H
        [-2.164526, -0.958307, 0.461966],  # H
        [-2.100884, 0.669410, 0.520327],  # H
        [2.362115, -0.284664, -0.124638],  # H
        [-4.581938, -0.009519, -0.484937],  # O (CO2)
        [-3.715310, -0.021766, -1.283957],  # C (CO2)
        [-3.285676, -0.007001, -2.380158],  # O
    ]
)

# 最適化後の構造
opt_coords = np.array(
    [
        [1.59521641984219, -0.23436327012506, -0.70008688394837],  # O
        [-2.03997325780028, -0.12676240026975, -0.10313308298181],  # N
        [-0.76895037265932, -0.15047040385740, -0.81804870689051],  # C (MEA)
        [0.42851838803848, -0.22451384496066, 0.10920458992258],  # C (MEA)
        [-0.70559812581428, 0.74806411408681, -1.42880937636032],  # H
        [-0.77094200512502, -1.00850406338704, -1.48740495515346],  # H
        [0.37135789636333, -1.13137032440317, 0.71732505468357],  # H
        [0.43044321549322, 0.63826836330510, 0.78092758021522],  # H
        [-2.16429602113822, -0.95828133702555, 0.46157961677174],  # H
        [-2.10068831365741, 0.66944008404988, 0.51992661572356],  # H
        [2.36210552442730, -0.28456401192413, -0.12415181972501],  # H
        [-4.58194465816717, -0.00967451899293, -0.48436989517430],  # O (CO2)
        [-3.71553876844233, -0.02175907317195, -1.28358748961964],  # C (CO2)
        [-3.28669197166034, -0.00682264082413, -2.38005158836319],  # O
    ]
)


def calc_distance(coords, i, j):
    """2つの原子間の距離を計算"""
    diff = coords[i] - coords[j]
    return np.linalg.norm(diff)


print("=" * 60)
print("構造変化の分析")
print("=" * 60)

# 主要な距離を計算
print("\n主要な原子間距離 (Ångström):")
print("-" * 60)
print(f"{'距離':<30} {'初期値':<15} {'最適化後':<15}")
print("-" * 60)

# N-CO2距離
n_idx = 1
c_co2_idx = 12
d_n_cco2_init = calc_distance(scan013_coords, n_idx, c_co2_idx)
d_n_cco2_opt = calc_distance(opt_coords, n_idx, c_co2_idx)
print(f"{'N - C(CO2)':<30} {d_n_cco2_init:<15.4f} {d_n_cco2_opt:<15.4f}")

# C=O 結合 (CO2)
o1_idx = 11
o2_idx = 13
d_co_init = calc_distance(scan013_coords, c_co2_idx, o1_idx)
d_co_opt = calc_distance(opt_coords, c_co2_idx, o1_idx)
print(f"{'C=O (CO2-O1)':<30} {d_co_init:<15.4f} {d_co_opt:<15.4f}")

d_co2_init = calc_distance(scan013_coords, c_co2_idx, o2_idx)
d_co2_opt = calc_distance(opt_coords, c_co2_idx, o2_idx)
print(f"{'C=O (CO2-O2)':<30} {d_co2_init:<15.4f} {d_co2_opt:<15.4f}")

# 全原子の変位
print("\n" + "-" * 60)
print("各原子の最大変位 (Å):")
print("-" * 60)

atom_names = [
    "O(water)",
    "N",
    "C1(MEA)",
    "C2(MEA)",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H",
    "H(OH)",
    "O(CO2)",
    "C(CO2)",
    "O(CO2)",
]
displacements = []

for i in range(len(scan013_coords)):
    diff = scan013_coords[i] - opt_coords[i]
    disp = np.linalg.norm(diff)
    displacements.append(disp)
    if i < 14:
        print(f"{atom_names[i]:<20} {disp:.6f}")

max_disp = max(displacements)
mean_disp = np.mean(displacements)
print(f"\n最大変位: {max_disp:.6f} Å")
print(f"平均変位: {mean_disp:.6f} Å")

print("\n" + "=" * 60)
print("結論：非常に小さな変位 → scan.013 はすでに最小値付近")
print("=" * 60)
