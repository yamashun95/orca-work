#!/bin/bash

# 反復的最適化スクリプト (遷移状態探索)
# 効率化版：OptTS と Freq を分離
# - OptTS だけで最適化（Freq は計算しない）
# - 最適化完了後、その最終構造に対してのみ Freq を単発実行
# - Freq の結果から虚振動をチェック
# - 虚振動が目標を超えていれば負モードに沿って構造を変位させて次のイテレーションへ

set -euo pipefail

MAX_ITER=25
CUTOFF=-1.0    # cm^-1 以下を虚振動と判定（-1.0で小さな負振動も検出）
SCALING=0.25   # モード沿い変位のスケール（Bohr単位のまま適用）
TARGET_NEG=1   # TS を想定

ORCA="$HOME/orca_6_1_1/orca"
INPUT_OPT="tsopt_iter_opt.inp"     # OptTS only
INPUT_FREQ="tsopt_iter.inp"         # Freq only
XYZ="tsopt.xyz"
PROP="tsopt_iter.property.txt"
FREQ_PARSER="parse_freq_from_log.py"
DISPLACER="displace_along_mode.py"

for ((iter=1; iter<=MAX_ITER; iter++)); do
  echo "========== Iteration $iter (Optimization) =========="

  # 1. OptTS 実行（Freq なし）
  $ORCA $INPUT_OPT > tsiter_${iter}_opt.log 2>&1 || {
    echo "ERROR: OptTS calculation failed at iteration $iter" >&2
    exit 1
  }

  # 最適化構造は tsopt_iter_opt.xyz に出力される
  if [ ! -f "tsopt_iter_opt.xyz" ]; then
    echo "ERROR: Optimized geometry tsopt_iter_opt.xyz not found" >&2
    exit 1
  fi
  
  # OptTS の出力を tsopt.xyz にコピー（次のイテレーションやFreq用）
  cp tsopt_iter_opt.xyz tsopt.xyz

  echo "========== Iteration $iter (Frequency) =========="

  # 2. Freq 単発実行（OptTS 後の構造に対して）
  $ORCA $INPUT_FREQ > tsiter_${iter}_freq.log 2>&1 || {
    echo "ERROR: Freq calculation failed at iteration $iter" >&2
    exit 1
  }

  if [ ! -f "$PROP" ]; then
    echo "ERROR: Property file $PROP not found" >&2
    exit 1
  fi

  # 3. 振動数を解析
  FREQ_LOG="tsiter_${iter}_freq.log"
  FREQ_REPORT="tsiter_${iter}_freqs.txt"
  
  # set -e を一時的に無効化（終了コード1は正常な動作）
  set +e
  python "$FREQ_PARSER" \
    --log "$FREQ_LOG" \
    --cutoff "$CUTOFF" \
    --target "$TARGET_NEG" \
    --output "$FREQ_REPORT"
  set -e

  # 成功判定（目標の虚振動数に到達）
  if grep -q "STATUS: CONVERGED" "$FREQ_REPORT"; then
    cp "$XYZ" tsopt_final.xyz
    echo "SUCCESS: converged to TS with $TARGET_NEG imaginary frequency after $iter iterations"
    exit 0
  fi
  
  # 虚振動が多すぎる場合、最も負の振動に沿って構造を変位
  if grep -q "STATUS: TOO_MANY" "$FREQ_REPORT"; then
    # 負振動インデックスを取得（最初の要素が最も負）
    NEG_INDICES=$(grep "Negative indices:" "$FREQ_REPORT" | sed 's/Negative indices: \[//; s/\]//' | tr -d ' ')
    FIRST_NEG=$(echo "$NEG_INDICES" | cut -d',' -f1)
    
    if [ -z "$FIRST_NEG" ]; then
      echo "ERROR: No negative mode found despite TOO_MANY status" >&2
      exit 1
    fi
    
    echo "Displacing along mode $FIRST_NEG (most negative)"
    
    python "$DISPLACER" \
      --xyz "$XYZ" \
      --log "$FREQ_LOG" \
      --mode "$FIRST_NEG" \
      --scaling "$SCALING" \
      --output "${XYZ}.tmp" || {
      echo "ERROR: Displacement failed" >&2
      exit 1
    }
    
    mv "${XYZ}.tmp" "$XYZ"
    echo "Structure displaced, continuing to iteration $((iter+1))"
    
  elif grep -q "STATUS: TOO_FEW" "$FREQ_REPORT"; then
    echo "ERROR: TOO_FEW negative frequencies (need $TARGET_NEG for TS)" >&2
    exit 1
  fi
done

echo "ERROR: Maximum iterations ($MAX_ITER) reached without convergence" >&2
exit 1
