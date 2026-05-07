#!/usr/bin/env python3
import os
import glob
import shutil
import subprocess
import sys

# ============================================================
# ★ 設定ここから ★
# ============================================================

# 入力データの親ディレクトリ
INPUT_BASE = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Output_global0429test/"

# 出力先ディレクトリ
OUTPUT_DIR = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/0429test/pTtest/"

# コンフィグファイル（JSON）が保存されている親ディレクトリ
CONF_BASE_DIR = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/Resolution/conf/" # ★環境に合わせて変更してください

# 解析タスク名
TASK_CMD = "o2-analysis-my-upc-muon-pair-resolutioneta-global-muon"

# 共有メモリサイズ
SHM_SIZE = "10000000000"

# サンプルリスト: (サブディレクトリ名, ファイルプレフィックス, 出力ファイル名, コンフィグファイル名)
SAMPLES = [
    ("jpsi-coh",       "jpsi-coh-*.root",        "jpsi-coh.root",       "jpsi-coh-conf.json"),
    #("jpsi-incoh",     "jpsi-incoh-*.root",       "jpsi-incoh.root",     "jpsi-incoh-conf.json"),
    #("psi2s-coh",      "psi2s-coh-*.root",        "psi2s-coh.root",      "psi2s-coh-conf.json"),
    #("psi2s-incoh",    "psi2s-incoh-*.root",      "psi2s-incoh.root",    "psi2s-incoh-conf.json"),
]

# ============================================================
# ★ 設定ここまで ★
# ============================================================

LIST_FILENAME = "skimmed_aod_list.txt"

os.makedirs(OUTPUT_DIR, exist_ok=True)

results = []  # (sample_name, success, message)

for (subdir, pattern, outname, conf_filename) in SAMPLES:
    target_dir = os.path.join(INPUT_BASE, subdir)
    search_pattern = os.path.join(target_dir, pattern)
    root_files = sorted(glob.glob(search_pattern))
    
    # コンフィグファイルの絶対パス
    conf_filepath = os.path.join(CONF_BASE_DIR, conf_filename)

    print(f"\n{'='*60}")
    print(f"[Sample] {subdir}")
    print(f"{'='*60}")

    if not root_files:
        msg = f"ROOTファイルが見つかりません: {search_pattern}"
        print(f"  [SKIP] {msg}")
        results.append((outname, False, msg))
        continue

    if not os.path.exists(conf_filepath):
        msg = f"コンフィグファイルが見つかりません: {conf_filepath}"
        print(f"  [SKIP] {msg}")
        results.append((outname, False, msg))
        continue

    # AODリストを作成
    list_filepath = os.path.join(target_dir, LIST_FILENAME)
    with open(list_filepath, "w") as lf:
        for rf in root_files:
            lf.write(f"{os.path.abspath(rf)}\n")
    print(f"  ファイル数: {len(root_files)}")
    print(f"  リスト: {list_filepath}")
    print(f"  コンフィグ: {conf_filepath}")

    # 解析実行コマンドの構築
    # ※ O2でJSONを読む場合は json:// プレフィックスを付けるのが標準的です
    cmd = [
        TASK_CMD,
        "--aod-file", f"@{list_filepath}",
        "--configuration", f"json://{conf_filepath}",
        "--shm-segment-size", SHM_SIZE,
        "-b",
    ]
    print(f"  実行: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        msg = f"クラッシュ (exit={e.returncode})"
        print(f"  [ERROR] {msg}")
        results.append((outname, False, msg))
        continue

    # AnalysisResults.root を OUTPUT_DIR/<outname> に移動
    src = "AnalysisResults.root"
    dst = os.path.join(OUTPUT_DIR, outname)
    if not os.path.exists(src):
        msg = "AnalysisResults.root が見つかりません"
        print(f"  [ERROR] {msg}")
        results.append((outname, False, msg))
        continue

    shutil.move(src, dst)
    print(f"  [OK] 保存: {dst}")
    results.append((outname, True, dst))

# ============================================================
# サマリー
# ============================================================
print(f"\n{'='*60}")
print("  実行サマリー")
print(f"{'='*60}")
success_count = 0
for (name, ok, msg) in results:
    status = "✓" if ok else "✗"
    print(f"  [{status}] {name}: {msg}")
    if ok:
        success_count += 1
print(f"\n  完了: {success_count}/{len(results)} サンプル")
print(f"  出力先: {OUTPUT_DIR}")