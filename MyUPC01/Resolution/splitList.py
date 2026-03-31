import os
import re

# --- 設定 ---
input_file = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/mumu_low_list.txt"
output_dir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/List/mumu"
max_files_per_list = 50

# 出力先ディレクトリの作成
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

run_dict = {}

# 1. ファイルの読み込みとラン番号ごとの分類
with open(input_file, "r") as f:
    for line in f:
        path = line.strip()
        if not path:
            continue
        
        # 正規表現でラン番号を抽出
        # "mumu_low/" の直後にある数字の連続をラン番号として取得する
        match = re.search(r'mumu_low/(\d+)/', path)
        if match:
            run_number = match.group(1)
            if run_number not in run_dict:
                run_dict[run_number] = []
            run_dict[run_number].append(path)
        else:
            print(f"Warning: ラン番号が抽出できませんでした -> {path}")

# 2. リストの分割と書き出し
for run_number, paths in run_dict.items():
    total_files = len(paths)
    
    # max_files_per_list (50) ごとにリストを分割して処理
    for i in range(0, total_files, max_files_per_list):
        chunk = paths[i:i + max_files_per_list]
        
        # 50以上のファイルがあり分割される場合はサフィックス(連番)をつける
        suffix = f"_{i // max_files_per_list}" if total_files > max_files_per_list else ""
        output_filename = os.path.join(output_dir, f"mumu-low-{run_number}{suffix}.txt")
        
        with open(output_filename, "w") as out_f:
            for p in chunk:
                out_f.write(p + "\n")
                
print(f"処理が完了しました。分割されたファイルは '{output_dir}' ディレクトリ内にあります。")