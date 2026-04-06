#!/bin/bash
# ============================================================
#  RunAnalysis.sh
#  pT 分解能解析を Step 0 → Step 1 → Step 2 → Step 3 の順で実行する。
#  alienv 環境内であれば、どのディレクトリからでも実行可能。
#
#  使い方:
#    bash /home/takuma/work/alice/O2Physics/MyUPC01/Resolution/macro/RunAnalysis.sh
#
#  一部のステップだけを再実行したい場合:
#    root -l -b -q '/home/takuma/work/alice/O2Physics/MyUPC01/Resolution/macro/Step0_CheckEventCount.C'
#    root -l -b -q '/home/takuma/work/alice/O2Physics/MyUPC01/Resolution/macro/Step1_MergeAndView2D.C'
#    root -l -b -q '/home/takuma/work/alice/O2Physics/MyUPC01/Resolution/macro/Step2_CrystalBallFit.C'
#    root -l -b -q '/home/takuma/work/alice/O2Physics/MyUPC01/Resolution/macro/Step3_ResolutionCurve.C'
# ============================================================

set -e  # エラーが起きたら即座に止める

MACRO_DIR="/home/takuma/work/alice/O2Physics/MyUPC01/Resolution/macro"

echo "====================================================="
echo " pT Resolution Analysis"
echo "====================================================="

echo ""
echo ">>> [Step 0] Checking event counts per pT bin..."
root -l -b -q "${MACRO_DIR}/Step0_CheckEventCount.C"
echo "    Step 0 done."

echo ""
echo ">>> [Step 1] Merging files and viewing 2D histogram..."
root -l -b -q "${MACRO_DIR}/Step1_MergeAndView2D.C"
echo "    Step 1 done."

echo ""
echo ">>> [Step 2] Crystal Ball fit per pT_true bin..."
root -l -b -q "${MACRO_DIR}/Step2_CrystalBallFit.C"
echo "    Step 2 done."

echo ""
echo ">>> [Step 3] Drawing resolution curve..."
root -l -b -q "${MACRO_DIR}/Step3_ResolutionCurve.C"
echo "    Step 3 done."

echo ""
echo "====================================================="
echo " All steps finished!"
echo " Output files are in: /media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/"
echo "   Step0_EventCount_pT.png"
echo "   Step0_EventCount_pT_each.png"
echo "   Step1_2D_Merged.png"
echo "   Step1_merged.root"
echo "   Step2_CrystalBallFit_bin??.png  (one per pT bin)"
echo "   Step2_fitresults.root"
echo "   Step3_ResolutionCurve.png"
echo "   Step3_BiasCurve.png"
echo "   Step3_Combined.png"
echo "   Step3_ResolutionCurve.root"
echo "====================================================="
