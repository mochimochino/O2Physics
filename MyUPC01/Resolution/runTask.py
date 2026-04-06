import os
import json
import glob
import subprocess

template_json_path = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/conf/psi2s-incoh-fd-conf.json"  # Change

list_files = glob.glob("/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/List/psi2s-incoh-fd/psi2s-incoh-fd-*.txt") # Change


failed_log_path = "failed_tasks.log"
failed_tasks = []

try:
    with open(template_json_path, 'r') as f:
        base_conf = json.load(f)
except FileNotFoundError:
    print(f"エラー: テンプレートとなる {template_json_path} が見つかりません。")
    exit(1)

for list_file in list_files:
    basename = os.path.splitext(os.path.basename(list_file))[0]
    
    base_conf["internal-dpl-aod-reader"]["aod-file-private"] = f"@{list_file}"
    
    temp_json_path = f"conf_{basename}.json"
    with open(temp_json_path, 'w') as f:
        json.dump(base_conf, f, indent=4)
        
    print(f"\n--- 処理開始: {basename} ---")
    
    cmd = [
        "o2-analysis-ud-upc-cand-producer-muon",
        "--configuration", f"json://{temp_json_path}",
        "--aod-writer-keep", "dangling",
        "--aod-writer-resfile", basename,
        "-b",
        "--shm-segment-size", "12000000000"
    ]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"完了: {basename}")
    except subprocess.CalledProcessError as e:
        print(f"エラー発生: {basename} の処理中にタスクがクラッシュしました。")
        print(e)
        failed_tasks.append(basename)
    
    os.remove(temp_json_path)

print("\n--- 全タスクの実行が終了しました ---")

if failed_tasks:
    print(f"警告: {len(failed_tasks)} 件のタスクが失敗しました。")
    with open(failed_log_path, 'w') as f:
        for failed_task in failed_tasks:
            f.write(f"{failed_task}\n")
    print(f"失敗したタスクのリストを '{failed_log_path}' に保存しました。再実行の際にご活用ください。")
else:
    print("すべてのタスクが正常に完了しました。")