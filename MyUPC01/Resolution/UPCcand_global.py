import os
import json
import glob
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

# 1. 全体設定
BASE_DIR = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/"
LIST_BASE_DIR = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/List/"  # リストのベースディレクトリを変更
OUTPUT_BASE_DIR = os.path.join(BASE_DIR, "Output_global0420")
RETRY_LIST_DIR = os.path.join(BASE_DIR, "List_Retry0420") # 分割リストの保存先
MAX_WORKERS = 40            # 並列実行数
MAX_FILES_PER_RETRY = 5     # 再実行時の1リストあたりの最大ファイル数

DATASET_MAP = {
    "jpsi-coh": "jpsi-coh-conf.json",
    "jpsi-incoh": "jpsi-incoh-conf.json",
    "mumu-high": "mumu-high-conf.json",
    "mumu-low": "mumu-low-conf.json",
    "mumu-mid": "mumu-mid-conf.json",
    "psi2s-coh": "psi2S-coh-conf.json",
    "psi2s-coh-fd": "psi2s-coh-fd-conf.json",
    "psi2s-incoh": "psi2s-incoh-conf.json",
    "psi2s-incoh-fd": "psi2s-incoh-fd-conf.json"
}

def run_o2_task(dataset_name, conf_filename, list_filepath, task_basename):
    """1つのテキストファイルリストに対するO2タスクを実行する関数"""
    template_json_path = os.path.join(BASE_DIR, "conf", conf_filename)
    output_dir = os.path.join(OUTPUT_BASE_DIR, dataset_name)
    os.makedirs(output_dir, exist_ok=True)
    
    temp_json_path = os.path.join(output_dir, f"conf_{task_basename}.json")

    try:
        with open(template_json_path, 'r') as f:
            base_conf = json.load(f)
        
        base_conf["internal-dpl-aod-reader"]["aod-file-private"] = f"@{list_filepath}"
        
        with open(temp_json_path, 'w') as f:
            json.dump(base_conf, f, indent=4)

        cmd = [
            "o2-analysis-ud-upc-cand-producer-global-muon",
            "--configuration", f"json://{temp_json_path}",
            "--aod-writer-keep", "dangling",
            "--aod-writer-resfile", task_basename,
            "-b",
            "--shm-segment-size", "12000000000"
        ]
        
        subprocess.run(
            cmd, 
            check=True, 
            cwd=output_dir, 
            stdout=subprocess.DEVNULL, 
            stderr=subprocess.PIPE
        )
        
        os.remove(temp_json_path)
        return True, task_basename, None
        
    except subprocess.CalledProcessError as e:
        error_msg = e.stderr.decode('utf-8') if e.stderr else str(e)
        if os.path.exists(temp_json_path):
            os.remove(temp_json_path)
        return False, task_basename, error_msg
    except Exception as e:
        if os.path.exists(temp_json_path):
            os.remove(temp_json_path)
        return False, task_basename, str(e)


def main():
    failed_log_path = os.path.join(BASE_DIR, "failed_tasks.log")
    retry_failed_log_path = os.path.join(BASE_DIR, "failed_tasks_retry.log")
    
    tasks_to_run = []
    failed_tasks_info = []

    # ==========================================
    # フェーズ1: 初回並列実行
    # ==========================================
    for dataset, conf_file in DATASET_MAP.items():
        # リストの取得元を LIST_BASE_DIR に変更
        list_dir = os.path.join(LIST_BASE_DIR, dataset)
        list_files = glob.glob(os.path.join(list_dir, "*.txt"))
        
        if not list_files:
            print(f"警告: {list_dir} にテキストファイルが見つかりません。")
            continue
            
        for list_filepath in list_files:
            basename = os.path.splitext(os.path.basename(list_filepath))[0]
            tasks_to_run.append((dataset, conf_file, list_filepath, basename))
            
    print(f"--- フェーズ1: 全 {len(tasks_to_run)} 件のタスクを {MAX_WORKERS} 並列で開始します ---")

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {
            executor.submit(run_o2_task, t[0], t[1], t[2], t[3]): t for t in tasks_to_run
        }
        
        for future in as_completed(futures):
            original_task_info = futures[future] 
            success, task_basename, error_msg = future.result()
            
            if success:
                print(f"[完了] {original_task_info[0]} / {task_basename}")
            else:
                print(f"[エラー] {original_task_info[0]} / {task_basename} がクラッシュしました。")
                failed_tasks_info.append(original_task_info)

    if failed_tasks_info:
        with open(failed_log_path, 'w') as f:
            for task in failed_tasks_info:
                f.write(f"{task[3]}\n")
    else:
        print("すべてのタスクが正常に完了しました。再実行フェーズはスキップします。")
        return

    # ==========================================
    # フェーズ2: 失敗タスクの分割と再実行
    # ==========================================
    print(f"\n--- フェーズ2: 失敗した {len(failed_tasks_info)} 件のタスクを分割して再実行します ---")
    retry_tasks_to_run = []

    for dataset, conf_file, original_list_filepath, basename in failed_tasks_info:
        if not os.path.exists(original_list_filepath):
            print(f"警告: 元のリストが見つかりません -> {original_list_filepath}")
            continue
            
        with open(original_list_filepath, 'r') as f:
            paths = [line.strip() for line in f if line.strip()]
            
        retry_dataset_dir = os.path.join(RETRY_LIST_DIR, dataset)
        os.makedirs(retry_dataset_dir, exist_ok=True)
        
        for i in range(0, len(paths), MAX_FILES_PER_RETRY):
            chunk = paths[i:i + MAX_FILES_PER_RETRY]
            sub_index = i // MAX_FILES_PER_RETRY
            new_basename = f"{basename}_{sub_index}"
            new_list_file = os.path.join(retry_dataset_dir, f"{new_basename}.txt")
            
            with open(new_list_file, "w") as out_f:
                for p in chunk:
                    out_f.write(p + "\n")
                    
            retry_tasks_to_run.append((dataset, conf_file, new_list_file, new_basename))

    print(f"分割後の再実行タスク: 計 {len(retry_tasks_to_run)} 件")
    final_failed_tasks = []

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {
            executor.submit(run_o2_task, t[0], t[1], t[2], t[3]): t for t in retry_tasks_to_run
        }
        
        for future in as_completed(futures):
            original_task_info = futures[future]
            success, task_basename, error_msg = future.result()
            
            if success:
                print(f"[再実行-完了] {original_task_info[0]} / {task_basename}")
            else:
                print(f"[再実行-エラー] {original_task_info[0]} / {task_basename} が再度クラッシュしました。")
                final_failed_tasks.append((task_basename, error_msg))

    # ==========================================
    # 最終結果の出力
    # ==========================================
    print("\n--- 全プロセスの実行が終了しました ---")
    if final_failed_tasks:
        print(f"警告: 再実行でも {len(final_failed_tasks)} 件のサブタスクが失敗しました。")
        with open(retry_failed_log_path, 'w') as f:
            for task_basename, error_msg in final_failed_tasks:
                f.write(f"Task: {task_basename}\nError:\n{error_msg}\n{'-'*40}\n")
        print(f"再実行で失敗したリストとエラー詳細を '{retry_failed_log_path}' に保存しました。")
    else:
        print("再実行したすべてのサブタスクが正常に完了しました。")

if __name__ == "__main__":
    main()