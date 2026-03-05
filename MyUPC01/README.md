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

## 解析マクロ (`AnalysisMassReal.C`)
`MyUPCTask.cxx` で作成された `AnalysisResults.root` を用いて、実際の不変質量分布の描画とバックグラウンドの差し引き（Signal Extraction）を行うための ROOT マクロです。

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
