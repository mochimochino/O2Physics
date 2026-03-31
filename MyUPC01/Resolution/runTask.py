import os
import json
import glob
import subprocess

# --- 設定 ---
# 提示いただいた元のJSONファイルをテンプレートとして使用します
template_json_path = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/conf/mumu-low-conf.json" 

# 実行するディレクトリ (現在のディレクトリ内の .txt を対象とする)
list_files = glob.glob("/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/List/mumu-low/mumu-low-*.txt")

# 失敗したタスクを記録するファイルのパスと、記録用のリスト
failed_log_path = "failed_tasks.log"
failed_tasks = []

# テンプレートとなるJSONの読み込み
try:
    with open(template_json_path, 'r') as f:
        base_conf = json.load(f)
except FileNotFoundError:
    print(f"エラー: テンプレートとなる {template_json_path} が見つかりません。")
    exit(1)

# 各リストファイルに対して処理を実行
for list_file in list_files:
    # ベース名の取得 (例: "mumu-low-544013.txt" -> "mumu-low-544013")
    basename = os.path.splitext(os.path.basename(list_file))[0]
    
    # 1. JSONの設定を動的に書き換える
    base_conf["internal-dpl-aod-reader"]["aod-file-private"] = f"@{list_file}"
    
    # 2. 一時的なJSONファイルとして保存する
    temp_json_path = f"conf_{basename}.json"
    with open(temp_json_path, 'w') as f:
        json.dump(base_conf, f, indent=4)
        
    print(f"\n--- 処理開始: {basename} ---")
    
    # 3. 実行するコマンドの構築
    cmd = [
        "o2-analysis-ud-upc-cand-producer-muon",
        "--configuration", f"json://{temp_json_path}",
        "--aod-writer-keep", "dangling",
        "--aod-writer-resfile", basename,
        "-b",
        "--shm-segment-size", "12000000000"
    ]
    
    # 4. タスクの実行
    try:
        subprocess.run(cmd, check=True)
        print(f"完了: {basename}")
    except subprocess.CalledProcessError as e:
        print(f"エラー発生: {basename} の処理中にタスクがクラッシュしました。")
        print(e)
        # エラーが発生したベース名をリストに追加
        failed_tasks.append(basename)
    
    # (オプション) 実行が終わった一時JSONファイルを削除したい場合は以下のコメントアウトを外す
    os.remove(temp_json_path)

print("\n--- 全タスクの実行が終了しました ---")

# 5. 失敗したタスクの記録と出力
if failed_tasks:
    print(f"警告: {len(failed_tasks)} 件のタスクが失敗しました。")
    with open(failed_log_path, 'w') as f:
        for failed_task in failed_tasks:
            f.write(f"{failed_task}\n")
    print(f"失敗したタスクのリストを '{failed_log_path}' に保存しました。再実行の際にご活用ください。")
else:
    print("すべてのタスクが正常に完了しました。")