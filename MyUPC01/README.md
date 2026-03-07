# MyUPC01 - Custom UPC Analysis Module

## 概要 (Overview)
このモジュールは、ALICE O2Physicsフレームワークにおいて、UPC (Ultra-peripheral Collision) イベントにおける前方ラピディティ領域（Forward Rapidity）のミューオン対（di-muon）不変質量分布を求めるためのカスタム解析モジュールです。

ベースとして `UDTutorial_06.cxx` を参考にしつつ、以下の点を拡張・安全化しています。

## 含まれるファイル
- `MyUPCTask.cxx`: 解析のメインコード。ミューオンのトラッキング情報とペアの運動学（pt, eta, rapidity, mass）からヒストグラムを作成します。
- `CMakeLists.txt`: このモジュールを `O2Physics` のビルドプロセスに組み込むためのCMake設定です。

## 技術的・実装上の改善ポイント (Differences from UDTutorial_06)
- **配列外アクセス（Out-of-bounds）の防止**: 元のチュートリアルコードでは、イベント内に前方トラックが「ちょうど2本」あることを前提にして配列アクセスを行っていたため、3本以上、あるいは1本以下のトラックが含まれるAO2Dデータフレームを読み込んだ際にクラッシュ（Segmentation Fault）を引き起こす可能性がありました。
- **任意のトラック数への対応**: イベントごとの前方トラック数が2本以上ある場合でも、すべての可能な「ペア（ペアリング）」を組み合わせて解析するようにループ処理を修正しています。（1本以下の場合は安全にスキップされます）
- **質量カット幅の調整**: 汎用的な $\mu\mu$ 不変質量分布を見れるように、特定のJ/$\psi$ピークだけでなく質量範囲を調整しやすくしています。

## コンパイル方法
`O2Physics` のルートディレクトリにて以下を実行してください：
```bash
aliBuild build O2Physics -j 2
```

## 実行方法
コンパイル後にエイリアス環境に入り、以下のコマンドで実行できます。（AO2D.root は実際のデータファイルのパスに置き換えてください）

```bash
alienv enter O2Physics/latest-MyUPC01-o2
o2-analysis-my-upc-01 --aod-file /media/takuma/ESD-EAWA/Data/UPCcandMuon/AO2D_merged.root -b
```

## 適用している主要なカット（Selection Criteria）

現行の `MyUPCTask.cxx` では以下のALICE UPCミューオン解析のカット条件を適用しています。
（CodiMDを参照）

## 各ミューオンの横運動量（pT）抽出機能
特定の不変質量領域（例：J/$\psi$ピーク付近）にあるミューオン対を構成する「各々のミューオンの$p_T$分布」を出力する機能が備わっています。
既存の `MyUPCTask` 内でペアの質量領域判定を行い、以下のヒストグラムに各ミューオンの$p_T$を個別に詰めています。
- `PtEachMuonUnlikeMassRegion`: Unlike Signペアのうち、指定した不変質量領域に入った各ミューオンの$p_T$成分
- `PtEachMuonLikeMassRegion`: Like Signペアのうち、指定した不変質量領域に入った各ミューオンの$p_T$成分

対象となる質量領域は、Configurableパラメータとして外部から柔軟に変更可能です。デフォルトでは J/$\psi$ 付近（2.8〜3.4 GeV）に設定されています。
実行時に任意の領域（例えば 2.0 〜 4.0 GeV）へ変更する場合は、以下のようにオプションを渡して実行します：

```bash
o2-analysis-my-upc-01 --aod-file /path/to/AO2D.root -b --massMin 2.0 --massMax 4.0
```

## 特定の質量領域におけるダイミューオン（ペア）の横運動量解析
Coherent / Incoherent J/$\psi$ などの物理過程を分離・解析するためには、特定の不変質量領域にある**ミューオン対（ダイミューオン）そのものの横運動量（$p_T$）分布**を評価する必要があります。
本モジュールでは、タスク実行時に 質量 vs $p_T$ の2次元ヒストグラム（`MassVsPtUnlike`, `MassVsPtLike`）を出力します。出力後、専用のマクロを用いて任意の質量領域の $p_T$ 分布を射影（Projection）して抽出することが可能です。

### ZDC情報を用いたトポロジー分類（中性子放出）
ZDC（Zero Degree Calorimeter: ZNA / ZNC）のエネルギー情報を用いて、どちらの原子核が電離によって壊れたか（中性子を放出したか）を判定し、事象を分類します。これによりCoherent（原子核が壊れにくい）とIncoherent（原子核が壊れやすい）の分離精度を高めることができます。

タスク実行時に、以下の4つのトポロジーごとに独立した 質量 vs $p_T$ の2次元ヒストグラム（`MassVsPtUnlike_0n0n` など）が自動的に生成されます。
- `0n0n` (topology 0): 両方の原子核が壊れなかった（ZNA, ZNC 共に閾値以下）
- `Xn0n` (topology 1): A側の原子核のみ壊れた（ZNAのみ閾値より大きい）
- `0nXn` (topology 2): C側の原子核のみ壊れた（ZNCのみ閾値より大きい）
- `XnXn` (topology 3): 両方の原子核が壊れた（ZNA, ZNC 共に閾値より大きい）

ZDCのエネルギー閾値（デフォルトは `1.0` ADC）も実行時のオプションで変更可能です：
```bash
# ZNAとZNCの閾値を 5.0 に設定する場合
o2-analysis-my-upc-01 --aod-file /path/to/AO2D.root -b --cutZNAEnergy 5.0 --cutZNCEnergy 5.0
```

### 解析マクロ (`AnalysisPtReal.C`)
出力された `AnalysisResults.root` を読み込み、特定の質量領域（デフォルトでは 2.8〜3.4 GeV）のダイミューオン横運動量分布を抽出するための ROOT マクロです。

#### 実行方法
1. マクロ内の `massMin` と `massMax` の値を、見たい質量領域に合わせて適宜修正します。
2. `AnalysisResults.root` が存在するディレクトリで以下を実行します：
```bash
root -l -q AnalysisPtReal.C
```
3. 実行結果として以下が出力されます：
   - `RawPt_MassRegion.png`: 指定領域の Unlike Sign と Like Sign の生の $p_T$ 分布比較
   - `SignalPt_MassRegion.png`: 上記から Like Sign バックグラウンドを差し引いた純粋なシグナルの $p_T$ 分布
   - 標準出力に、指定範囲での各 Sign のカウント数（Yield）が表示されます。

## 不変質量分布の解析マクロ (`AnalysisMassReal.C`)
`MyUPCTask.cxx` で作成された `AnalysisResults.root` を用いて、全体の不変質量分布の描画とバックグラウンドの差し引き（Signal Extraction）を行うための ROOT マクロです。

### 処理の概要
1. **Raw Histogramsの描画**: 
   異符号ペア（Unlike Sign: $+-$）と、同符号ペア（Like Sign: $++$ および $--$ の和）の不変質量分布を同一キャンバスに重ねて描画します。
2. **シグナル抽出（Signal Extraction）**: 
   Unlike Sign の分布から、組み合わせバックグラウンド（Combinatorial Background）としての Like Sign の分布を引き去り、純粋な物理シグナルを取り出します。

### バックグラウンド評価の数式
通常、不変質量分布における組み合わせバックグラウンド $N_{bkg}$ は、正の同符号ペア $N_{++}$ と負の同符号ペア $N_{--}$ の幾何平均を用いて以下のように評価されます。

$$ N_{bkg} = 2 \sqrt{N_{++} N_{--}} $$

しかし、今回の `MyUPCTask.cxx` の実装では、Like Sign ($++$ と $--$) を区別せずに1つのヒストグラム（`MMuonLike`）として足し合わせて出力しています。

$$ N_{like} = N_{++} + N_{--} $$

$N_{++} \approx N_{--}$ の場合、算術平均と幾何平均はほぼ等しくなるため、以下のように近似できます。

$$ 2 \sqrt{N_{++} N_{--}} \approx N_{++} + N_{--} = N_{like} $$

したがって、このマクロでは単純に Unlike Sign から Like Sign を差し引くことでシグナル $N_{signal}$ を抽出しています。

$$ N_{signal} = N_{+-} - N_{like} $$

### 実行方法
`AnalysisResults.root` が存在するディレクトリ（例: `/home/takuma/work/alice/O2Physics/MyUPC01/`）で、以下のコマンドを実行します。
```bash
root -l -q AnalysisMassReal.C
```
出力として `RawMass.png`（生のヒストグラム）と `Unlike_bkg.png`（シグナル抽出後のヒストグラム）が生成されます。
