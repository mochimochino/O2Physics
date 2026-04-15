# Bin_width — ビン幅決定マクロ 説明書

このディレクトリには、unfolding 解析に使うビン幅を決定するための ROOT マクロが含まれる。  
各マクロは `hGenPt2`（Generator レベル）と `hRecoPt2`（Reconstructed レベル）の  
1D ヒストグラムをリビンし、続けて `Make2DHistFromBinnedHists.C` を呼んで  
2D レスポンス行列のリビンと `Cell/Gen [%]` の出力を行う。

---

## ファイル一覧

| ファイル名 | 主な役割 |
|---|---|
| `RebinEqualWidth.C` | 等幅ビン分割 |
| `RebinFlatStats.C` | 等統計ビン分割 |
| `RebinConstantDiff.C` | 等差拡大ビン分割 |
| `RebinTopDown.C` | TopDown 法（Purity/Stability 下限指定）によるビン自動決定 |
| `CheckPurityStability.C` | 任意のビン境界で Purity / Stability を計算 |
| `CompareBeforeAfter.C` | リビン前後の 1D ヒストグラムを比較 |

---

## 1. `RebinEqualWidth.C` — 等幅ビン分割

### ビン計算方法

$$
\text{nBins} = \left\lfloor \frac{x_{\max} - x_{\min}}{\Delta w} \right\rfloor
$$

$$
\text{ビン境界} = x_{\min},\; x_{\min}+\Delta w,\; x_{\min}+2\Delta w,\; \ldots
$$

- すべてのビン幅が **一定** (`binWidth`)
- 最終ビン上端は `xMin + nBins * binWidth`（整数個に揃える）

### 主なパラメータ

| 引数 | デフォルト値 | 説明 |
|---|---|---|
| `binWidth` | `0.01` | 1ビンの幅 |
| `xMin` | `0.0` | 開始値 |
| `xMax` | `0.065` | 終了値（最大） |

### 特徴
- シンプルで扱いやすい
- 統計量の少ない高 $p_T^2$ 領域でビン数が多くなりすぎる場合がある

### 出力ヒスト名
- `hGenPt2_rebin`, `hRecoPt2_rebin`

---

## 2. `RebinFlatStats.C` — 等統計ビン分割

### ビン計算方法

Gen ヒスト (`hGenPt2`) の積分が各ビンで等しくなるようにビン境界を決める。

$$
\text{step} = \frac{\text{totalIntegral}}{\text{nTargetBins}}
$$

左から積算し、累積値が `step` の整数倍を超えるたびにビン境界を置く：

```
cumSum += hGen->GetBinContent(i)
if (cumSum >= targetSum):
    ビン境界 = bin i+1 の左端
    targetSum += step
```

- 最後のビン上端は必ず `xMax` に設定される

### 主なパラメータ

| 引数 | デフォルト値 | 説明 |
|---|---|---|
| `nTargetBins` | `6` | 目標ビン数 |
| `xMin` | `0.0` | 開始値 |
| `xMax` | `0.065` | 終了値 |

### 特徴
- 各ビンの Gen 統計量がほぼ均一 → 統計誤差が均一になる
- 低 $p_T^2$ 側（イベントが多い）は幅が狭く、高 $p_T^2$ 側（イベントが少ない）は幅が広くなる

### 出力ヒスト名
- `hGenPt2_flat`, `hRecoPt2_flat`

---

## 3. `RebinConstantDiff.C` — 等差拡大ビン分割

### ビン計算方法

初期ビン幅 $w_0$ から始め、ビンを追加するたびに幅を $\Delta w$ ずつ大きくする（等差数列）：

$$
w_k = w_0 + k \cdot \Delta w \quad (k = 0, 1, 2, \ldots)
$$

$$
x_{k+1} = x_k + w_k
$$

`respectOverflow = true` の場合、最後のビン上端を強制的に `xMax` に揃える。

### 主なパラメータ

| 引数 | デフォルト値 | 説明 |
|---|---|---|
| `initialWidth` | `0.005` | 最初のビン幅 $w_0$ |
| `deltaW` | `0.005` | ビンごとの幅の増分 $\Delta w$ |
| `xMin` | `0.0` | 開始値 |
| `xMax` | `0.065` | 終了値 |
| `respectOverflow` | `true` | 最後のビンを `xMax` に合わせるか |

### 特徴
- 低 $p_T^2$ 側は細かく、高 $p_T^2$ 側は粗くなる
- Resolution や統計量の $p_T^2$ 依存性に合わせてチューニングしやすい

### 出力ヒスト名
- `hGenPt2_cdiff`, `hRecoPt2_cdiff`

---

## 4. `RebinTopDown.C` — TopDown 法（Purity/Stability 下限指定）

### ビン計算方法

高解像度の 2D レスポンス行列を使い、**Purity と Stability が指定した下限を満たす最大範囲**をビンとして上から順に確定する。

```
1. upperBin = maxPt2 に最も近いファインビン
2. lowerBin を upperBin から 0 に向けて 1 つずつ下げる
3. 区間 [lowerBin, upperBin] の Purity・Stability・Stats を計算
4. 全条件を満たしたら:
     ビン境界 = lowerBin の左端
     upperBin = lowerBin - 1  (次のビンの探索開始点)
5. x=0 まで繰り返す
```

$$
\text{Purity}   = \frac{\text{diag}}{\sum_j M_{j,\text{iy}}} \geq \text{targetPurity}
$$

$$
\text{Stability} = \frac{\text{diag}}{\sum_j M_{\text{ix},j}} \geq \text{targetStability}
$$

$$
\text{sumReco} \geq \text{minStats}, \quad \text{sumTrue} \geq \text{minStats}
$$

### 主なパラメータ

| 引数 | デフォルト値 | 説明 |
|---|---|---|
| `targetPurity` | `0.45` | Purity の下限 |
| `targetStability` | `0.45` | Stability の下限 |
| `minStats` | `1000` | ビン内最低統計数 |
| `maxPt2` | `0.065` | ビンの最大値 |

### 特徴
- データ駆動でビン幅を決めるため、Purity・Stability を保証できる
- 統計量が少ない場合、ビン数が少なくなる（または 0 になることがある）
- `targetPurity` / `targetStability` / `minStats` を緩めるとビン数が増える

### 出力ヒスト名
- `hGenPt2_topdown`, `hRecoPt2_topdown`

---

## 5. `CheckPurityStability.C` — Purity / Stability チェック

Rebin マクロではなく、**ユーザーが直接指定したビン境界**に対して  
Purity と Stability を計算するユーティリティ。

### 計算方法

$$
\text{Purity}_i = \frac{\text{diag}_{ii}}{\sum_j M_{ji}} \quad (\text{Reco ビン} i \text{ に落ちた全イベントのうち、Gen も } i \text{ だった割合})
$$

$$
\text{Stability}_i = \frac{\text{diag}_{ii}}{\sum_j M_{ij}} \quad (\text{Gen ビン} i \text{ のうち、Reco も } i \text{ に落ちた割合})
$$

ここで $M_{ij}$ は高解像度 2D レスポンス行列 (`hResponseMatrixPairPt2`) の $(i, j)$ 成分。

### 主なパラメータ

| 引数 | デフォルト値 | 説明 |
|---|---|---|
| `inputFile` | Step1_merged.root | 高解像度行列を含む ROOT ファイル |
| `outputFile` | Logic1_PurityStability.root | 出力 ROOT ファイル |
| `binEdges` | `{0.0, 0.010, 0.065}` | ビン境界のベクター |

---

## 6. `CompareBeforeAfter.C` — リビン前後の比較

リビン前（細かいビン）とリビン後のヒストグラムを重ねて描画し PNG に保存する。  
どのリビン手法でも共通して使用可能。

### 主なパラメータ

| 引数 | 説明 |
|---|---|
| `originalFile` | 元の（細かいビン）ROOT ファイル |
| `rebinnedFile` | リビン後の ROOT ファイル（各 Rebin マクロの出力） |
| `histName` | 比較するヒスト名（`"hGenPt2"` または `"hRecoPt2"`） |

---

## 各マクロの実行フロー

```
root -l RebinEqualWidth.C
         │
         ▼
  Gen / Reco ヒストをリビン
  → ROOT ファイルに保存 (hGenPt2_rebin, hRecoPt2_rebin)
  → PNG を保存
         │
         ▼
  Make2DHistFromBinnedHists() を呼び出し
  → リビン済みヒストのビン境界で 2D レスポンス行列をリビン
  → Purity, Stability を計算
  → Cell/Gen [%] の 2D ヒストを作成・描画
  → ROOT / PDF に保存
```

同様の流れが `RebinFlatStats.C`, `RebinConstantDiff.C` でも実行される。

---

## ビン分割方法の比較

| 方法 | ビン幅の変化 | 統計の均一性 | チューニングのしやすさ |
|---|---|---|---|
| `RebinEqualWidth` | 一定 | 不均一（低 $p_T^2$ 偏重） | 簡単（幅1つで指定） |
| `RebinFlatStats` | 統計に応じて変化 | 均一 | 中程度（bin 数で指定） |
| `RebinConstantDiff` | 等差数列で増加 | やや不均一 | 柔軟（$w_0$, $\Delta w$ で調整） |
| `RebinTopDown` | データ駆動（自動決定） | Purity/Stability を保証 | 3つのしきい値で制御 |

---

## 関連ファイル

- [`../Make2DHistFromBinnedHists.C`](../Make2DHistFromBinnedHists.C) — 各 Rebin マクロから呼ばれる 2D ヒスト作成マクロ
- [`../Step2_DetermineBinning.C`](../Step2_DetermineBinning.C) — Purity / Stability 基準で自動的にビン境界を決定するマクロ（TopDown 法）
