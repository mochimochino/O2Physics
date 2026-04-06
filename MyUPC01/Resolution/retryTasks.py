import os
import json
import subprocess

failed_log_path = "failed_tasks.log"
template_json_path = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/conf/psi2s-incoh-fd-conf.json" # Change
list_dir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/List/psi2s-incoh-fd" # Change

max_files_per_retry_list = 10
retry_failed_log_path = "failed_tasks_retry.log"

retry_failed_tasks = []

try:
    with open(template_json_path, 'r') as f:
        base_conf = json.load(f)
except FileNotFoundError:
    print(f"エラー: テンプレート {template_json_path} が見つかりません。")
    exit(1)

if not os.path.exists(failed_log_path):
    print(f"エラー: 失敗ログ {failed_log_path} が見つかりません。再実行するタスクがありません。")
    exit(0)

with open(failed_log_path, 'r') as f:
    failed_tasks = [line.strip() for line in f if line.strip()]

print(f"--- 再実行対象: {len(failed_tasks)} 件のリストをさらに分割して処理します ---")

for basename in failed_tasks:
    original_list_file = os.path.join(list_dir, f"{basename}.txt")
    
    if not os.path.exists(original_list_file):
        print(f"警告: 元のリストファイルが見つかりません -> {original_list_file}")
        continue
        
    with open(original_list_file, 'r') as f:
        paths = [line.strip() for line in f if line.strip()]
        
    total_files = len(paths)
    
    for i in range(0, total_files, max_files_per_retry_list):
        chunk = paths[i:i + max_files_per_retry_list]
        
        sub_index = i // max_files_per_retry_list
        new_basename = f"{basename}_{sub_index}"
        new_list_file = os.path.join(list_dir, f"{new_basename}.txt")
        
        with open(new_list_file, "w") as out_f:
            for p in chunk:
                out_f.write(p + "\n")
                
        base_conf["internal-dpl-aod-reader"]["aod-file-private"] = f"@{new_list_file}"
        temp_json_path = f"conf_{new_basename}.json"
        
        with open(temp_json_path, 'w') as jf:
            json.dump(base_conf, jf, indent=4)
            
        print(f"\n[再実行] 処理開始: {new_basename} ({len(chunk)} ファイル)")
        
        cmd = [
            "o2-analysis-ud-upc-cand-producer-muon",
            "--configuration", f"json://{temp_json_path}",
            "--aod-writer-keep", "dangling",
            "--aod-writer-resfile", new_basename,
            "-b",
            "--shm-segment-size", "12000000000"
        ]

        try:
            subprocess.run(cmd, check=True)
            print(f"完了: {new_basename}")
        except subprocess.CalledProcessError as e:
            print(f"エラー発生: {new_basename} がクラッシュしました。")
            retry_failed_tasks.append(new_basename)
            
        if os.path.exists(temp_json_path):
            os.remove(temp_json_path)

print("\n--- 再実行タスクの処理がすべて終了しました ---")

if retry_failed_tasks:
    print(f"警告: 再実行でも {len(retry_failed_tasks)} 件のタスクが失敗しました。")
    with open(retry_failed_log_path, 'w') as f:
        for failed_task in retry_failed_tasks:
            f.write(f"{failed_task}\n")
    print(f"再実行で失敗したリストを '{retry_failed_log_path}' に保存しました。")
    print("※これらのファイルには破損したROOTファイルが含まれている可能性があります。")
else:
    print("再実行したすべてのタスクが正常に完了しました。")