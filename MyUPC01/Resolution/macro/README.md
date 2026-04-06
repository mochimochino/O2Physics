# pT 分解能解析マクロ

このディレクトリには、UPC µ pair データの横運動量（pT）分解能を
ステップごとに確認する ROOT マクロが含まれています。

## ディレクトリ構成

```
pairpT/
├── jpsi-incoh.root
├── jpsi-coh.root
├── psi2s-*.root
├── mumu-*.root
├── Step1_2D_Merged.png        ← Step1 の出力
├── Step1_merged.root          ← Step1 の出力 (Step2 が読む)
├── Step2_CrystalBallFit_bin??.png ← Step2 の出力（各 pT ビン）
├── Step2_fitresults.root      ← Step2 の出力 (Step3 が読む)
├── Step3_ResolutionCurve.png  ← Step3 の出力
├── Step3_BiasCurve.png        ← Step3 の出力
├── Step3_Combined.png         ← Step3 の出力（σ + μ を1枚に）
├── Step3_ResolutionCurve.root ← Step3 の出力 (グラフ入り)
└── macro/
    ├── Step1_MergeAndView2D.C
    ├── Step2_CrystalBallFit.C
    ├── Step3_ResolutionCurve.C
    ├── RunAnalysis.sh
    └── README.md  ← このファイル
```

## 実行方法

**前提**: `alienv enter O2Physics/latest-MyUPC_mini-o2` で alienv 環境に入っていること。

```bash
cd /media/takuma/ESD-EAWA/Data/UPCcandMuon/MC/0331/pairpT/macro

# 全ステップを一括実行
bash RunAnalysis.sh

# または各ステップを個別実行
root -l -b -q 'Step1_MergeAndView2D.C'
root -l -b -q 'Step2_CrystalBallFit.C'
root -l -b -q 'Step3_ResolutionCurve.C'
```

## 各ステップの説明

### Step 1: マージ → 2D ヒストグラム確認
- **入力**: 各 MC サンプル（9 ファイル）
  - `hPtResoVsPtTrue_PostCut` を読み込む
- **処理**: 全ファイルをマージ → pT_true 軸をリビンニング
- **出力**:
  - `Step1_2D_Merged.png` … COLZ 表示で分解能の全体像を確認
  - `Step1_merged.root` … リビン済み TH2D を保存（Step2 で使用）

### Step 2: Crystal Ball フィット（各 pT ビン）
- **入力**: `Step1_merged.root`
- **処理**: 各 pT_true ビンで 1D 射影 → Crystal Ball 関数でフィット
  - CB 関数パラメータ: N（規格化）, μ（平均/バイアス）, σ（分解能）, α（テール開始点）, n（テール指数）
  - 1st pass: Gauss で大まかな peak/sigma を推定
  - 2nd pass: Crystal Ball でフィット（±3σ 範囲）
- **出力**:
  - `Step2_CrystalBallFit_bin??.png` … 各ビンのフィット結果を画像保存
  - `Step2_fitresults.root` … gr_sigma（σ vs pT）, gr_mu（μ vs pT）

### Step 3: 分解能カーブ描画
- **入力**: `Step2_fitresults.root`
- **出力**:
  - `Step3_ResolutionCurve.png` … σ（分解能）vs pT_true カーブ
  - `Step3_BiasCurve.png` … μ（運動量バイアス）vs pT_true カーブ
  - `Step3_Combined.png` … σ と μ を上下 2 パッドに並べた複合図
  - `Step3_ResolutionCurve.root` … グラフを ROOT ファイルに保存

## 物理的な意味

| 量 | 意味 |
|---|---|
| `σ_CB` | CB フィットのコア幅 = **pT 分解能**（ビン幅決定の基礎） |
| `μ_CB` | CB フィットの平均 = **運動量スケールのバイアス** |
| `α, n` | テール形状パラメータ（放射損失・多重散乱の tail を記述） |

分解能曲線 σ(pT) は、:
- ビン幅の最適化（σ より広いビンを取る）
- 系統誤差評価（smearing 計算）
に使用します。
