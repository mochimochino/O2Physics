import os
import glob
import subprocess

target_dir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/psi2s-incoh" #Change

list_filename = "skimmed_aod_list.txt"


search_pattern = os.path.join(target_dir, "psi2s-incoh-*.root") #Change
root_files = glob.glob(search_pattern)

if not root_files:
    print(f"エラー: {search_pattern} に一致するROOTファイルが見つかりません。")
    exit(1)

list_filepath = os.path.join(target_dir, list_filename)
with open(list_filepath, "w") as f:
    for root_file in sorted(root_files):
        f.write(f"{os.path.abspath(root_file)}\n")

print(f"--- リスト作成完了 ---")
print(f"対象ファイル数: {len(root_files)}")
print(f"リストファイル: {list_filepath}\n")

cmd = [
    "o2-analysis-my-upc-muon-resolution",
    "--aod-file", f"@{list_filepath}",
    "-b"
]

print("--- 解析タスク実行開始 ---")
print(f"実行コマンド: {' '.join(cmd)}")

try:
    subprocess.run(cmd, check=True)
    print("\n--- 解析タスクが正常に完了しました ---")
except subprocess.CalledProcessError as e:
    print(f"\nエラー発生: 解析タスクの実行中にクラッシュしました。")
    print(f"終了コード: {e.returncode}")