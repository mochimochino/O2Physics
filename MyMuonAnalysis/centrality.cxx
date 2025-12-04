// myCentralityTask.cxx
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/Centrality.h" // Centrality定義のインクルード

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;

struct myCentralityTask {
  // ヒストグラムの定義
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    // Centrality (0-100%) 用の軸定義
    const AxisSpec axisCent{100, 0, 100, "Centrality FT0M (%)"};
    
    // ヒストグラムの追加
    histos.add("centHistogram", "Centrality FT0M", kTH1F, {axisCent});
  }

  // 衝突テーブル(Collisions)とCentralityテーブル(CentFT0Ms)を結合(Join)して処理
  // Run 2データの場合は CentFT0Ms を CentRun2V0Ms に変更してください
  void process(soa::Join<Collisions, CentFT0Ms>::iterator const& col)
  {
    // Centralityの値を取得
    float cent = col.centFT0M(); // Run 2の場合は col.centRun2V0M()

    // ログに出力（デバッグ用：大量に出るので注意）
    // LOGF(info, "Event Centrality: %f", cent);

    // ヒストグラムに詰める
    histos.fill(HIST("centHistogram"), cent);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myCentralityTask>(cfgc)};
}