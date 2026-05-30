# UPCMuonAnalysisMC — タスク設定変数リファレンス

**ファイル:** `UPCMuonAnalysisMC.cxx`  
**タスク名:** `upc-muon-analysis`  
**O2Physics コマンド:** `o2-analysis-my-upc-muon-analysis-mc`

---

## 概要

UPC 光核反応におけるジミューオンペアの **分解能 (Resolution)** および **再構成効率 (Efficiency)** を
同一ジョブで解析する統合タスク。
MC 真値ループ (`processMCTrue`) と Reco ループ (`processMcReco`) を持つ。

---

## Configurable 変数一覧

### 1. トラック選択

| 変数名 | 型 | デフォルト値 | 説明 |
|--------|------|------------|------|
| `reqMatchMFT` | `int` | `2` | イベントあたりの要求トラック本数（ペア解析では 2 を指定） |
| `reqTrackType` | `int` | `0` (GlobalMuonTrack) | 解析対象のトラックタイプ。`0`: GlobalMuonTrack（MCH + MFT マッチング済み）、`3`: MuonStandaloneTrack（MCH のみ） |

---

### 2. シングルトラック アクセプタンス（MC 真値・Reco 共通）

| 変数名 | 型 | デフォルト値 | 説明 |
|--------|------|------------|------|
| `etaMin` | `float` | `-3.6` | シングルミューオンの擬似ラピディティ下限（厳密な不等号: η > etaMin） |
| `etaMax` | `float` | `-2.5` | シングルミューオンの擬似ラピディティ上限（厳密な不等号: η < etaMax） |
| `pTmin` | `float` | `0.3` GeV/c | シングルミューオンの横運動量下限（pT ≥ pTmin） |

> **注意:** `etaMin`/`etaMax` のカットは MC 真値ループ (`processMCTrue`) にも
> Reco ループ (`processMcReco`) にも同じ値で適用される。

---

### 3. Reco クオリティカット（Reco ループのみ）

| 変数名 | 型 | デフォルト値 | 説明 |
|--------|------|------------|------|
| `rAbsMin` | `float` | `17.6` cm | アブソーバー端での半径方向位置 R の下限 |
| `rAbsMax` | `float` | `89.5` cm | アブソーバー端での半径方向位置 R の上限 |
| `maxChi2` | `float` | `100` | グローバルトラックフィットの χ² 上限 |
| `maxChi2MatchMCHMFT` | `float` | `100` | MCH–MFT マッチングの χ² 上限（`reqTrackType = 0` の場合のみ適用） |

> **pDCA カット:** コードに直接実装されており Configurable ではない。
> R < 26.5 cm のとき pDCA < 350 GeV/c・cm、R ≥ 26.5 cm のとき pDCA < 200 GeV/c·cm。

---

### 4. ジミューオンペア 運動学カット（MC 真値・Reco 共通）

| 変数名 | 型 | デフォルト値 | 説明 |
|--------|------|------------|------|
| `pairPtMax` | `float` | `0.25` GeV/c | ペアの横運動量上限（コヒーレント J/ψ 用）。インコヒーレント解析では大きな値に変更する |
| `pairRapidityMin` | `float` | `-3.6` | ペアのラピディティ下限（厳密な不等号） |
| `pairRapidityMax` | `float` | `-2.5` | ペアのラピディティ上限（厳密な不等号） |
| `pairMassMin` | `float` | `1.0` GeV/c² | ペアの不変質量下限（厳密な不等号） |
| `pairMassMax` | `float` | `10.0` GeV/c² | ペアの不変質量上限（厳密な不等号） |

---

### 5. ヒストグラム ビニング: ペア pT / pT²

| 変数名 | 型 | デフォルト値 | 説明 |
|--------|------|------------|------|
| `nBinsPt` | `int` | `1000` | ペア pT 軸のビン数 |
| `ptMax` | `float` | `0.5` GeV/c | ペア pT ヒストグラムの上限 |
| `nBinsPt2` | `int` | `5000` | ペア pT² 軸のビン数 |
| `pt2Max` | `float` | `0.5` GeV²/c² | ペア pT² ヒストグラムの上限 |

---

### 6. ヒストグラム ビニング: 不変質量

| 変数名 | 型 | デフォルト値 | 説明 |
|--------|------|------------|------|
| `nBinsMass` | `int` | `300` | 質量軸のビン数 |
| `massAxisMin` | `float` | `2.0` GeV/c² | 質量ヒストグラムの下限 |
| `massAxisMax` | `float` | `5.0` GeV/c² | 質量ヒストグラムの上限（J/ψ 解析向け；広い質量範囲が必要なら変更） |

---

### 7. ヒストグラム ビニング: シングルトラック

| 変数名 | 型 | デフォルト値 | 説明 |
|--------|------|------------|------|
| `nBinsPhi` | `int` | `360` | φ 軸のビン数 |
| `nBinsEta` | `int` | `300` | η / ラピディティ軸のビン数 |
| `nBinsP` | `int` | `1000` | 運動量 p 軸（および pT 分解能ヒスト）のビン数 |
| `pTrackMax` | `float` | `50.0` GeV/c | p ヒストグラムの上限 |
| `ptTrackMax` | `float` | `10.0` GeV/c | pT ヒストグラムの上限 |
| `pxpyTrackMax` | `float` | `10.0` GeV/c | px / py ヒストグラムの絶対値上限（±pxpyTrackMax の範囲） |

---

## 設定例（JSON コンフィグ）

```json
{
  "upc-muon-analysis": {
    "reqMatchMFT": 2,
    "reqTrackType": 0,
    "etaMin": -3.6,
    "etaMax": -2.5,
    "pTmin": 0.3,
    "rAbsMin": 17.6,
    "rAbsMax": 89.5,
    "maxChi2": 100.0,
    "maxChi2MatchMCHMFT": 100.0,
    "pairPtMax": 0.25,
    "pairRapidityMin": -3.6,
    "pairRapidityMax": -2.5,
    "pairMassMin": 1.0,
    "pairMassMax": 10.0,
    "nBinsPt": 1000,
    "ptMax": 0.5,
    "nBinsPt2": 5000,
    "pt2Max": 0.5,
    "nBinsMass": 300,
    "massAxisMin": 2.0,
    "massAxisMax": 5.0,
    "nBinsPhi": 360,
    "nBinsEta": 300,
    "nBinsP": 1000,
    "pTrackMax": 50.0,
    "ptTrackMax": 10.0,
    "pxpyTrackMax": 10.0
  }
}
```

---

## カットフロー (CutFlow) ビン説明

`hCutFlowMC`（MC 真値）と `hCutFlowReco`（Reco）で共通のビン定義：

| ビン | ラベル | 説明 |
|------|--------|------|
| 0 | All Cand | 全粒子 / 全コリジョン |
| 1 | Has Requested Track Type / PDG=mu | Reco: 指定トラックタイプが存在 / MC: PDG コードがミューオン |
| 2 | Pass Exact 2 Tracks | candidate レベルで reqMatchMFT 本（Reco のみ） |
| 3 | Pass Exact MatchMFT | quality カット後も reqMatchMFT 本（Reco のみ） |
| 4 | Track All | 各トラックのエントリ（Reco のみ） |
| 6 | Track Pass rAbs | rAbsMin ≤ R ≤ rAbsMax を通過 |
| 7 | Track Pass pDCA | pDCA カットを通過 |
| 8 | Track Pass MatchMFT | χ²（グローバル + MCH–MFT）カットを通過 |
| 9 | Track Pass Eta | etaMin < η < etaMax を通過 |
| 10 | Track Pass pT | pT ≥ pTmin を通過 |
| 11 | Pair All | 全ペアの組み合わせ |
| 12 | Pair Unlike-sign | 反対符号ペア |
| 13 | Pair Both MC Muon | 両トラックが MC ミューオン |
| 14 | Pair Pass Pt | pT < pairPtMax を通過 |
| 15 | Pair Pass Rapidity | pairRapidityMin < y < pairRapidityMax を通過 |
| 16 | Pair Pass Mass | pairMassMin < M < pairMassMax を通過 |
