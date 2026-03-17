# ALICE O2Physics: Luminosity Calculation for UPC Analysis

## 1. ルミノシティの計算式 (Calculation Formula)

ご自身の用いたデータ (AOD) における積分ルミノシティ（Integrated Luminosity, $\mathcal{L}_{\text{int}}$）は、以下の式で計算されます。

$$ \mathcal{L}_{\text{int}} = \frac{N_{\text{events}}}{\sigma_{\text{trigger}}} $$

- **$N_{\text{events}}$**: 解析に用いたAODの中で、基準となるトリガー（例: 最小バイアストリガーである `sel8` や、前方トリガー `FT0` など）に合格したイベントの総数です。
- **$\sigma_{\text{trigger}}$**: その基準トリガーの可視断面積（Visible Cross-Section）です。これはALICE実験グループが van der Meer (vdM) スキャンという特殊な手法を用いて事前に測定・決定している値です。（例: Run 3 pp 13.6 TeV におけるFT0トリガーの断面積など。具体的な値は解析する年のALICE公式情報やCCDBから取得します）

つまり、「特定のデータセット（Run）に含まれる基準イベント数を数え、それをALICE公式の断面積で割る」ことで求められます。

## 2. O2Physics内の既存のタスクについて

O2Physicsには、イベント数やルミノシティに関する情報を計算・出力する標準のタスクがすでに存在します。
主に以下の2つがよく使われます。

1. **`o2-analysis-event-selection`**
   - イベント選択（Vertex Zカットや `sel8` など）を通過したイベント数をカウントし、`AnalysisResults.root` 内の `EventSelection` ディレクトリ配下に `hCollCounter` などのヒストグラムを出力します。まずはここで $N_{\text{events}}$ を確認できます。
   
2. **`o2-analysis-je-luminosity-calculator`** や **`o2-analysis-lumi-stability`**
   - より詳細なルミノシティ計算や、時間経過（安定性）に伴うルミノシティを監視するためのタスクです。

**ただし、ユーザーが指定したローカルのAODリスト (`input_AODlist.txt`) に対して、「RunNumberごとのルミノシティ」を直接「横軸RunNumber、縦軸Lumi」の形で出力してくれる都合の良いシンプルなタスクは、標準タスクのままではすぐに出力しにくいため（ヒストグラムの軸がRunNumberになっていないことが多い）、自分で簡単なタスクを作るか、取得したイベント数を計算し直すのが一般的です。**

## 3. RunNumberごとのルミノシティを求めるためのアプローチ（推奨）

横軸を「RunNumber」、縦軸を「ルミノシティ」にするためには、各Runごとの基準イベント数 $N_{\text{events}}^{\text{Run}}$ を取得し、それを $\sigma_{\text{trigger}}$ で割る必要があります。

**アプローチ案:**
Runごとのカウントを記録する専用のカスタムタスク (`MyUPC_Luminosity.cxx` など) を新しく作成することをお勧めします。

このタスクでは以下を行います：
1. `aod::BCs` (Bunch Crossings) または `aod::Collisions` をループする。
2. そのイベントの `runNumber()` を取得する。
3. イベントセレクション (例えばFT0や特定のトリガー) をパスしているか判定する。
4. `TH1D` などのヒストグラムに、`Fill(runNumber)` として記録する。（これで各Runのイベント数がわかる）
5. 解析終了後のROOTマクロ（描画マクロ）で、そのヒストグラムの全ビンを $1/\sigma_{\text{trigger}}$ でスケーリングし、縦軸を「ルミノシティ」に変換して描画する。

もしこの「Runごとのイベント数をカウントしてルミノシティグラフを作成する専用のO2タスク (`MyUPC_Luminosity.cxx`)」を作成したい場合は、すぐにご用意できますのでお知らせください。
