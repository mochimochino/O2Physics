import os
import glob
import subprocess
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed

# ==========================================
# 1. 全体設定
# ==========================================
BASE_DIR = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon/test/"
LIST_BASE_DIR = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/List/"  
OUTPUT_BASE_DIR = os.path.join(BASE_DIR, "Output_QA_0427")
RETRY_LIST_DIR = os.path.join(BASE_DIR, "List_Retry_QA_0427")

MAX_WORKERS = 15            # メモリやCPUに応じて調整してください
MAX_FILES_PER_RETRY = 2     # 再実行時の1リストあたりの最大ファイル数

# 今回はJSON不要なので、処理したいデータセット（ディレクトリ名）のリストだけ定義
DATASETS = [
    "jpsi-coh",
    "jpsi-incoh"
    # "mumu-high",
    # "mumu-low",
]

# QAタスクのパラメータ（意図的にゆるくして全貌を見る）
QA_MAX_CHI2 = "100.0"
QA_MAX_DCAXY = "999.0"

# ==========================================

def run_o2_task(dataset_name, list_filepath, task_basename):
    """1つのテキストファイルリストに対するQAタスクを実行する関数"""
    
    # データセットごとの最終出力先ディレクトリ
    final_output_dir = os.path.join(OUTPUT_BASE_DIR, dataset_name)
    os.makedirs(final_output_dir, exist_ok=True)
    
    # 最終的なヒストグラムのファイル名
    final_root_file = os.path.join(final_output_dir, f"{task_basename}.root")
    
    # 並列実行時の AnalysisResults.root 競合を防ぐため、タスク専用の作業ディレクトリを作る
    work_dir = os.path.join(final_output_dir, f"work_{task_basename}")
    os.makedirs(work_dir, exist_ok=True)

    cmd = [
        "o2-analysis-my-upc-globaltrack-qa",
        "--aod-file", f"@{list_filepath}",
        "-b",
        "--QA_maxChi2", QA_MAX_CHI2,
        "--QA_maxDCAxy", QA_MAX_DCAXY
    ]
    
    try:
        # タスク固有の作業ディレクトリ(work_dir)でコマンドを実行
        subprocess.run(
            cmd, 
            check=True, 
            cwd=work_dir, 
            stdout=subprocess.DEVNULL, 
            stderr=subprocess.PIPE
        )
        
        # 成功したら、AnalysisResults.root を名前を変えて回収
        generated_root = os.path.join(work_dir, "AnalysisResults.root")
        if os.path.exists(generated_root):
            shutil.move(generated_root, final_root_file)
        else:
            raise FileNotFoundError(f"{generated_root} が生成されませんでした。")
        
        # 作業ディレクトリのお掃除
        shutil.rmtree(work_dir)
        return True, task_basename, None
        
    except subprocess.CalledProcessError as e:
        error_msg = e.stderr.decode('utf-8') if e.stderr else str(e)
        return False, task_basename, error_msg
    except Exception as e:
        return False, task_basename, str(e)


def main():
    failed_log_path = os.path.join(BASE_DIR, "failed_qa_tasks.log")
    retry_failed_log_path = os.path.join(BASE_DIR, "failed_qa_tasks_retry.log")
    
    tasks_to_run = []
    failed_tasks_info = []

    # ==========================================
    # フェーズ1: 初回並列実行
    # ==========================================
    for dataset in DATASETS:
        list_dir = os.path.join(LIST_BASE_DIR, dataset)
        list_files = glob.glob(os.path.join(list_dir, "*.txt"))
        
        if not list_files:
            print(f"警告: {list_dir} にテキストファイルが見つかりません。")
            continue
            
        for list_filepath in list_files:
            basename = os.path.splitext(os.path.basename(list_filepath))[0]
            tasks_to_run.append((dataset, list_filepath, basename))
            
    print(f"--- フェーズ1: 全 {len(tasks_to_run)} 件のQAタスクを {MAX_WORKERS} 並列で開始します ---")

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {
            executor.submit(run_o2_task, t[0], t[1], t[2]): t for t in tasks_to_run
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
                f.write(f"{task[2]}\n")
    else:
        print("すべてのタスクが正常に完了しました。再実行フェーズはスキップします。")
        return

    # ==========================================
    # フェーズ2: 失敗タスクの分割と再実行
    # ==========================================
    print(f"\n--- フェーズ2: 失敗した {len(failed_tasks_info)} 件のタスクを分割して再実行します ---")
    retry_tasks_to_run = []

    for dataset, original_list_filepath, basename in failed_tasks_info:
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
            new_basename = f"{basename}_retry_{sub_index}"
            new_list_file = os.path.join(retry_dataset_dir, f"{new_basename}.txt")
            
            with open(new_list_file, "w") as out_f:
                for p in chunk:
                    out_f.write(p + "\n")
                    
            retry_tasks_to_run.append((dataset, new_list_file, new_basename))

    print(f"分割後の再実行タスク: 計 {len(retry_tasks_to_run)} 件")
    final_failed_tasks = []

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {
            executor.submit(run_o2_task, t[0], t[1], t[2]): t for t in retry_tasks_to_run
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
    else:
        print("再実行したすべてのサブタスクが正常に完了しました。")

if __name__ == "__main__":
    main()