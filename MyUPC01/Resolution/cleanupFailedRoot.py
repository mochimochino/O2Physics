import os

# 1. パスの設定（前回の実行スクリプトに完全準拠）
BASE_DIR = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/GlobalMuon"
OUTPUT_BASE_DIR = os.path.join(BASE_DIR, "Output_global0420")

DATASETS = [
    "jpsi-coh",
    "jpsi-incoh",
    "mumu-high",
    "mumu-low",
    "mumu-mid",
    "psi2s-coh",
    "psi2s-coh-fd",
    "psi2s-incoh",
    "psi2s-incoh-fd"
]

# ログファイルはBASE_DIR直下にあるため、パスを結合
LOG_FILES = [
    os.path.join(BASE_DIR, "failed_tasks.log"),
    os.path.join(BASE_DIR, "failed_tasks_retry.log")
]

def parse_failed_tasks(log_path):
    """ログファイルから失敗したタスクのベース名(basename)のみを抽出する"""
    basenames = set()
    with open(log_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # フェーズ2のログ形式: "Task: basename"
            if line.startswith("Task:"):
                basename = line.replace("Task:", "").strip()
                basenames.add(basename)
            # フェーズ1のログ形式: エラー詳細などを除外し、名前だけを抽出
            elif not line.startswith("Error:") and not line.startswith("-") and " " not in line and "Exception" not in line:
                basenames.add(line)
                
    return list(basenames)

def main():
    print(f"--- 一括クリーンアッププロセス開始 ---")
    print(f"解析結果ベースディレクトリ: {OUTPUT_BASE_DIR}")
    
    # 1. すべてのログファイルから削除対象のタスク名をグローバルに収集
    failed_basenames = set()
    for log_path in LOG_FILES:
        if os.path.exists(log_path):
            print(f"ログ読み込み: {os.path.basename(log_path)}")
            parsed_tasks = parse_failed_tasks(log_path)
            failed_basenames.update(parsed_tasks)
        else:
            print(f"情報: ログファイルが見つかりません -> {os.path.basename(log_path)}")
            
    if not failed_basenames:
        print("\n記録されている失敗タスクはありません。クリーンアップを終了します。")
        return

    print(f"\n合計 {len(failed_basenames)} 件のユニークな失敗タスクを抽出しました。削除処理に移行します。")
    total_deleted_count = 0

    # 2. 各データセットディレクトリ内を検索して該当ファイルを削除
    for dataset in DATASETS:
        target_dir = os.path.join(OUTPUT_BASE_DIR, dataset)
        
        if not os.path.exists(target_dir):
            continue
            
        dataset_deleted_count = 0
        
        for basename in failed_basenames:
            # 抽出したbasenameをもとに各ディレクトリ内での存在チェックを行う
            root_file_path = os.path.join(target_dir, f"{basename}.root")
            
            if os.path.exists(root_file_path):
                try:
                    os.remove(root_file_path)
                    print(f"  削除完了: [{dataset}] {basename}.root")
                    dataset_deleted_count += 1
                    total_deleted_count += 1
                except OSError as e:
                    print(f"  エラー: {root_file_path} の削除に失敗しました。({e})")
        
        if dataset_deleted_count > 0:
            print(f"  >> {dataset} フォルダ内で {dataset_deleted_count} 件のファイルを削除しました。\n")

    print(f"--- 一括クリーンアッププロセス終了 ---")
    print(f"合計 {total_deleted_count} 件の不要な .root ファイルを削除しました。")

if __name__ == "__main__":
    main()