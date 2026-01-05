#!/bin/bash
# 初回のみ実行：入力ファイルを準備し、スクリプトを実行

cd /home/yamamoto/Work/chemistry/orca-work/MEA_CO2_water_afify/tsopt

# 入力ファイルを確保
cp tsopt.inp tsopt_iter_opt.inp
sed -i 's|./scan.010.xyz|tsopt.xyz|g' tsopt_iter_opt.inp

# 初期XYZを設定
cp scan.010.xyz tsopt.xyz

# 以前の出力を削除（入力ファイルは消さない）
rm -f tsopt_iter.property.txt tsiter_*.* tsopt_iter_opt.xyz tsopt_iter_opt_trj.xyz

# スクリプト実行
./iterative_ts_opt.sh
