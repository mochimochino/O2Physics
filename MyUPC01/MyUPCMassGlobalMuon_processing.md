# MyUPCMassGlobalMuon.cxx 処理内容まとめ

`MyUPCMassGlobalMuonTask` は、UPC (Ultra-Peripheral Collision) イベントにおける
GlobalMuon (MFT-MCH-MID) トラックを用いたディミューオン不変質量解析タスク。
データ/MC両対応で、PROCESS_SWITCHにより4種類の処理モードを切り替える。

## 1. 全体構成

- 入力テーブル
  - `CandidatesFwd` = `UDCollisions` + `UDCollisionsSelsFwd` (UPC候補イベント)
  - `ForwardTracks` = `UDFwdTracks` + `UDFwdTracksExtra` (前方トラック)
  - `CompleteFwdTracks` = `ForwardTracks` + `UDMcFwdTrackLabels` (MC reco用)
  - `UDZdcsReduced` (ZDC情報、processDataWithZdcのみ)
  - `UDMcCollisions` / `UDMcParticles` (MC truth用)

- 出力ヒストグラム (HistogramRegistry)
  - 必須: `hMassUnlike`, `hPtUnlike`, `hPt2Unlike` (不変質量・pT・pT^2。良トラックがちょうど2本でmu+/mu-が1つずつのペアのみ、構造上unlike-signのみ)
  - MC truth対応: `hMassTruth`, `hPtTruth`, `hPt2Truth`
  - イベントQA: `eventCounter`, `hNumContrib`, `hNTracksTotal`, `hNGoodTracks`
  - カットフロー: `hCutFlow` (13ビン)
  - トラック種別QA: `hTrackType`
  - シングルトラックQA: `hRAbs`, `hPDCA`, `hChi2`, `hChi2MatchMCHMFT`
  - シングルミューオン運動学QA: `hEtaMuon`, `hPhiMuon`, `hPtMuon`
  - ペア運動学QA: `hRapidityPair`, `hPhiPair`
  - V0A QA: `hV0AAmp`
  - ZDC QA (processDataWithZdcのみ): `hEnergyZNA/ZNC`, `hEnergyZNAvsZNC`, `hTopologyCounter`

## 2. Configurables (主要カット値)

| カテゴリ | パラメータ | デフォルト | 説明 |
|---|---|---|---|
| 検出器スイッチ | `useV0ACut` | true | V0A同一BCアンプリチュード veto |
| | `cutV0AAmp` | 100.0 | V0A最大アンプリチュード |
| トラック種別 | `reqTrackType` | GlobalMuonTrack(0) | 要求するforward track type |
| トラック品質 | `rAbsMin/Max` | 17.6 / 89.5 cm | absorber end でのR範囲 |
| | `maxChi2` | 100 | global track chi2上限 |
| | `maxChi2MatchMCHMFT` | 100 | MCH-MFT match chi2上限 |
| シングルμ運動学 | `etaMin/Max` | -3.5 / -2.5 | シングルミューオンη範囲 |
| | `pTmin` | 0.3 GeV/c | シングルミューオンpT下限 |
| ペア運動学 | `pairPtMax` | 0.25 GeV/c | ペアpT上限 |
| | `pairRapidityMin/Max` | -3.5 / -2.5 | ペアrapidity範囲 |
| | `pairMassMin/Max` | 1.0 / 10.0 GeV/c^2 | ペア質量範囲 |
| ZDC | `cutZNAEnergy/cutZNCEnergy` | 1.0 TeV | 中性子タグのエネルギー閾値 |
| | `targetTopology` | -1 (no cut) | 要求する中性子トポロジー (0n0n/Xn0n/0nXn/XnXn) |

## 3. ヘルパー関数

- **`collectCandIDs`**: 全トラックを `udCollisionId()` ごとにグループ化し、候補イベント単位のトラックIDリストを作成。
- **`passV0ACut`**: 同一BC内のV0Aアンプリチュード最大値をヒストグラム化しつつ、`cutV0AAmp` 未満かどうかを返す（`useV0ACut`がtrueのときのみカットとして使用）。
- **`buildZdcMap`**: `UDZdcsReduced` から候補ID→ZDC情報(`ZDCinfo`)のマップを構築。ZNA/ZNC のtime窓(`kMaxZDCTime`=2ns)とエネルギー閾値から中性子トポロジークラス(0n0n/Xn0n/0nXn/XnXn)を判定。
- **`selectGoodTracks`**: 候補内の各トラックについて、
  1. `trackType` が `reqTrackType` と一致するか
  2. `rAtAbsorberEnd` が `[rAbsMin, rAbsMax]` 内か
  3. `pDca` が閾値 (rAbs<26.5なら350, それ以外200) 以下か
  4. `chi2` が `maxChi2` 以下か
  5. (GlobalMuonのみ) `chi2MatchMCHMFT` が `[0, maxChi2MatchMCHMFT]` 内か

  を順にチェックし、通過したトラックIDのリストを返す。各段階でQA/カットフローを記録。
- **`fillPairHistograms`**: 2トラックからTLorentzVectorを構成し、
  1. シングルμのη/φ/pTをQAに記録
  2. η範囲カット (`etaMin/Max`)
  3. シングルμ pTカット (`pTmin`)
  4. ペアpT/φ をQAに記録
  5. ペアpTカット (`pairPtMax`)
  6. ペアrapidityカット (`pairRapidityMin/Max`)
  7. ペア質量カット (`pairMassMin/Max`)
  8. 全カット通過後、`hMassUnlike`/`hPtUnlike`/`hPt2Unlike` を埋める
     (呼び出し元で既にmu+/mu-が1つずつであることを保証済みのため、unlike-sign判定は不要)

  各段階で `hCutFlow` を更新。
- **`processCandidates`**: 候補イベントごとのメインループ。
  1. トラック数2未満の候補をスキップ
  2. V0Aカット適用 (`useV0ACut`)
  3. (ZDC使用時) ZDC情報取得・トポロジーカウント・`targetTopology`カット
  4. `selectGoodTracks` で良トラック選定。**ちょうど2本でない場合はスキップ**
  5. 残った2本の電荷の和が0でない(=mu+/mu-が1つずつでない)場合はスキップ
  6. 上記を満たす場合のみ `fillPairHistograms` を1回実行

## 4. Process関数 (PROCESS_SWITCH)

| 関数名 | デフォルト | 入力 | 内容 |
|---|---|---|---|
| `processDataNoZdc` | **true** | `CandidatesFwd`, `ForwardTracks` | ZDCテーブル無しで実データ処理。`processCandidates` を `applyZdc=false` で呼ぶ |
| `processDataWithZdc` | false | `CandidatesFwd`, `ForwardTracks`, `UDZdcsReduced` | ZDC付き実データ。`buildZdcMap`でZDC情報マップ作成後、`applyZdc=true`で`processCandidates` |
| `processMcReco` | false | `CandidatesFwd`, `CompleteFwdTracks`, `UDMcCollisions`, `UDMcParticles` | MC reco処理。ZDC未使用、`processCandidates`を`applyZdc=false`で呼ぶ |
| `processMcTruth` | false | `UDMcCollisions`, `UDMcParticles` | MC truthレベル処理（下記参照） |

### processMcTruth の流れ

1. `UDMcParticles` から `|pdgCode| == 13` (muon) かつ η/pT 受容範囲内の粒子を、`udMcCollisionId` ごとにグループ化
2. 各MCコリジョン内で、unlike-sign (PDGコードの積が負) のμペアを全組み合わせで作成
3. ペアpT/rapidity/質量カットを適用
4. 通過したペアを `hMassTruth` / `hPtTruth` / `hPt2Truth` に記録

## 5. カットフロー (`hCutFlow`, 13ビン)

```
0  AllCand           : トラック2本以上の候補
1  V0A_pass          : V0Aカット通過
2  ZDCTopology       : ZDCトポロジーカット通過 (ZDC未使用時は素通り)
3  Track_All         : 全トラック (selectGoodTracks内でループされた回数)
4  Track_Type        : 要求トラックタイプ一致
5  Track_Quality     : 単トラック品質カット (rAbs/pDCA/chi2/chi2MFT) 通過
6  HasTwoGoodTracks  : 良トラックがちょうど2本、かつmu+/mu-が1つずつ
7  Pair_EtaCut       : シングルμ ηカット通過
8  Pair_PtCutMuon    : シングルμ pTカット通過
9  Pair_PtCut        : ペアpTカット通過
10 Pair_RapidityCut  : ペアrapidityカット通過
11 Pair_MassCut      : ペア質量カット通過
12 Filled            : 全カット通過、hMassUnlike等を充填
```

注: `HasTwoGoodTracks` は「良トラック数==2」かつ「電荷の和==0 (mu+とmu-が1つずつ)」の
両方を満たした場合のみカウントされる。これにより同符号(like-sign)ペアは原理上発生しない。

## 6. ZDC中性子トポロジー判定 (`buildZdcMap`)

- `timeZNA`/`timeZNC` の絶対値が `kMaxZDCTime` (2ns) 未満かつ非infのとき、その側に「タイミングあり」と判定
- タイミングありかつ `energyCommonZNA/C` がそれぞれ `cutZNAEnergy`/`cutZNCEnergy` を超える場合、その側を「中性子放出あり (Xn)」と判定
- 結果を `znClass` に格納: 0=0n0n, 1=Xn0n, 2=0nXn, 3=XnXn
- `targetTopology` が -1以外の場合、一致しない候補はスキップ
