#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/AnalysisCut.h"

#include "TFile.h"
#include "TH1F.h"

// 注意: 実プロジェクトのトラック型に合わせて置き換えてください。
// ここでは説明用に簡易 Track 構造体を定義しています。
struct Track {
  float pt;
  float eta;
  float phi;
  // Muon-specific vars used by matchedQualityCuts:
  float rAtAbsorberEnd; // -> VarManager::kMuonRAtAbsorberEnd
  float pDca;           // -> VarManager::kMuonPDca
  float chi2;           // -> VarManager::kMuonChi2
  float chi2MatchMCHMID; // -> VarManager::kMuonChi2MatchMCHMID
  float chi2MatchMCHMFT; // -> VarManager::kMuonChi2MatchMCHMFT
};

// Example: dummy event generator to illustrate usage.
// Replace with your real event reading/loop.
static std::vector<Track> makeDummyTracks()
{
  std::vector<Track> v;
  v.push_back({1.2f, -3.0f, 0.1f, 30.0f, 10.0f, 1.5f, 0.1f, 0.2f});
  v.push_back({0.6f, -3.8f, 2.1f, 25.0f, 50.0f, 2.5f, 0.2f, 0.3f});
  v.push_back({2.5f, -3.0f, -1.1f, 40.0f, 5.0f, 1.0f, 0.05f, 0.05f});
  // ... add more for testing
  return v;
}

int main(int argc, char** argv)
{
  // 出力 ROOT ファイルとヒストグラムを作る
  TFile fout("my_muonPt_output.root", "RECREATE");
  TH1F hPt("hPt", "Selected muon p_{T};p_{T} [GeV/c];counts", 100, 0.0, 20.0);
  TH1F hEta("hEta", "Selected muon #eta;#eta;counts", 100, -5.0, -2.0);
  TH1F hPhi("hPhi", "Selected muon #phi;#phi;counts", 64, -3.1416, 3.1416);

  // 1) CutsLibrary からカットを取得
  AnalysisCut* muonCut = o2::aod::dqcuts::GetAnalysisCut("matchedQualityCuts");
  if (!muonCut) {
    std::cerr << "Error: matchedQualityCuts not found\n";
    return 1;
  }

  // 2) VarManager の名前->index マップから必要な配列サイズを決める
  //    VarManager::fgVarNamesMap は CutsLibrary で使われている変数名マップ
  int maxIndex = -1;
  for (const auto& kv : VarManager::fgVarNamesMap) {
    // kv.second が変数インデックスであることを想定
    if (kv.second > maxIndex) {
      maxIndex = kv.second;
    }
  }
  if (maxIndex < 0) {
    std::cerr << "Error: VarManager::fgVarNamesMap seems empty\n";
    delete muonCut;
    return 1;
  }

  // 値配列をゼロ初期化して確保（IsSelected は float* を受け取る）
  std::vector<float> values(static_cast<size_t>(maxIndex + 1), 0.0f);

  // ここから実際のイベント・トラックループに入る（以下はダミー例）
  std::vector<Track> tracks = makeDummyTracks();
  for (const auto& track : tracks) {
    // --- values[] に VarManager の各インデックスに相当する値をセットする ---
    // 必須: matchedQualityCuts で参照される変数群を正しく埋めること
    // matchedQualityCuts の実装を見ると、以下の変数を参照しています：
    //   VarManager::kEta
    //   VarManager::kMuonRAtAbsorberEnd
    //   VarManager::kMuonPDca
    //   VarManager::kMuonChi2
    //   VarManager::kMuonChi2MatchMCHMID
    //   VarManager::kMuonChi2MatchMCHMFT
    // それ以外に PT を使う可能性があるカットを同時に使うなら VarManager::kPt も埋めてください。
    // ここでは pT/eta/phi もヒスト用に埋めます（ヒストは Track から直接でも可）。
    values[VarManager::kPt] = track.pt;
    values[VarManager::kEta] = track.eta;
    // VarManager に kPhi があるなら設定（IsSelected 内では使われませんが安全）
    if (VarManager::fgVarNamesMap.find("kPhi") != VarManager::fgVarNamesMap.end()) {
      // もし VarManager に kPhi 定義が文字列 "kPhi" で登録されているなら:
      values[VarManager::kPhi] = track.phi; // VarManager::kPhi が存在する前提で
    }
    values[VarManager::kMuonRAtAbsorberEnd] = track.rAtAbsorberEnd;
    values[VarManager::kMuonPDca] = track.pDca;
    values[VarManager::kMuonChi2] = track.chi2;
    values[VarManager::kMuonChi2MatchMCHMID] = track.chi2MatchMCHMID;
    values[VarManager::kMuonChi2MatchMCHMFT] = track.chi2MatchMCHMFT;

    // --- カットを評価 ---
    bool pass = muonCut->IsSelected(values.data());
    if (pass) {
      // 選択されたトラックをヒストグラムに入れる
      hPt.Fill(track.pt);
      hEta.Fill(track.eta);
      hPhi.Fill(track.phi);
    }
  } // end track loop

  // ヒストを書き出し
  fout.cd();
  hPt.Write();
  hEta.Write();
  hPhi.Write();
  fout.Close();

  // Clean up
  delete muonCut;

  std::cout << "Done. Output written to my_muonPt_output.root\n";
  return 0;
}