import os
import re

input_file = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/psi2s_incoh_fd_list.txt" #Change
output_dir = "/media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/List/psi2s-incoh-fd" #Change
max_files_per_list = 50

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

run_dict = {}

with open(input_file, "r") as f:
    for line in f:
        path = line.strip()
        if not path:
            continue
        
        match = re.search(r'psi2s_incoh_fd/(\d+)/', path) # change
        if match:
            run_number = match.group(1)
            if run_number not in run_dict:
                run_dict[run_number] = []
            run_dict[run_number].append(path)
        else:
            print(f"Warning: ラン番号が抽出できませんでした -> {path}")

for run_number, paths in run_dict.items():
    total_files = len(paths)
    
    for i in range(0, total_files, max_files_per_list):
        chunk = paths[i:i + max_files_per_list]
        
        suffix = f"_{i // max_files_per_list}" if total_files > max_files_per_list else ""
        output_filename = os.path.join(output_dir, f"psi2s-incoh-fd-{run_number}{suffix}.txt") # Change
        
        with open(output_filename, "w") as out_f:
            for p in chunk:
                out_f.write(p + "\n")
                
print(f"処理が完了しました。分割されたファイルは '{output_dir}' ディレクトリ内にあります。")