#!/bin/bash
# ORCA構造最適化計算の状況確認スクリプト

cd /home/yamamoto/Work/chemistry/orca-work/MEA_CO2_water_afify/backwardopt

echo "=========================================="
echo "ORCA構造最適化計算の実行報告"
echo "=========================================="
echo ""
echo "【実行内容】"
echo "✓ 初期構造: IRC最終構造 (backward_init.xyz)"
echo "  - 由来: IRCパス終端 (C-N距離: 5.994 Å)"
echo "  - エネルギー: -399.091837 Eh"
echo ""
echo "✓ 計算設定 (backward_opt.inp):"
echo "  - 関数形・基底関数: TPSS0 / def2-TZVP"
echo "  - 積分: RIJCOSX, TightSCF, SlowConv"
echo "  - 溶媒: 水 (CPCM, SMD=false)"
echo "  - 最適化条件: MaxIter=100, Convergence=Tight"
echo "  - 並列処理: nprocs=16"
echo ""
echo "【計算進捗】"
if ps aux | grep -q "[o]rca_leanscf_mpi"; then
    echo "✓ ステータス: 実行中 (MPI 16プロセスで稼働)"
    
    # ステップ数を確認
    STEP_COUNT=$(grep -c "Geometry step " backward_opt.log 2>/dev/null || echo "0")
    echo "  - 最適化ステップ: $STEP_COUNT 完了"
    
    # 最新エネルギー
    LATEST_ENERGY=$(grep "SCF Done:" backward_opt.log | tail -1 | awk '{print $(NF-1)}')
    if [ -n "$LATEST_ENERGY" ]; then
        echo "  - 最新エネルギー: $LATEST_ENERGY Eh"
    fi
    
    # RMS勾配
    LATEST_RMS=$(grep "RMS gradient" backward_opt.log | tail -1 | awk '{print $3}')
    if [ -n "$LATEST_RMS" ]; then
        echo "  - 最新RMS勾配: $LATEST_RMS (目標値: 0.00003)"
    fi
    
    echo ""
    echo "計算は継続中です。完了までお待ちください。"
else
    echo "✓ ステータス: 完了 ✓"
    
    # 完了確認
    if grep -q "OPTIMIZATION SUCCESSFUL" backward_opt.log 2>/dev/null; then
        echo "  最適化が正常に完了しました！"
        echo ""
        echo "【結果】"
        grep "OPTIMIZATION SUCCESSFUL" backward_opt.log -A 5
    else
        echo "  ログを確認してください。"
    fi
fi

echo ""
echo "ログファイル: backward_opt.log"
echo "出力構造ファイル: backward_opt.xyz"
echo "=========================================="
