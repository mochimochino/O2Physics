import os

target_dir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/mumu-high/" # Change

log_files = ["failed_tasks.log", "failed_tasks_retry.log"]

deleted_count = 0

print(f"--- 削除プロセス開始 ---")
print(f"対象ディレクトリ: {target_dir}")

for log_filename in log_files:
    log_path = os.path.join(target_dir, log_filename)
    
    if not os.path.exists(log_path):
        print(f"\n情報: ログファイルが見つかりません (スキップします) -> {log_filename}")
        continue
        
    print(f"\n--- ログファイル読み込み: {log_filename} ---")
    
    with open(log_path, 'r') as f:
        failed_basenames = [line.strip() for line in f if line.strip()]
        
    if not failed_basenames:
        print("記録されている失敗タスクはありません。")
        continue
        
    for basename in failed_basenames:
        root_file_path = os.path.join(target_dir, f"{basename}.root")
        
        if os.path.exists(root_file_path):
            try:
                os.remove(root_file_path)
                print(f"削除完了: {basename}.root")
                deleted_count += 1
            except OSError as e:
                print(f"エラー: {basename}.root の削除に失敗しました。({e})")
        else:
            print(f"スキップ: ファイルが既に存在しません -> {basename}.root")

print(f"\n--- 削除プロセス終了 ---")
print(f"合計 {deleted_count} 件の不要な .root ファイルを削除しました。")