#include "larreco/RecoAlg/TCAlg/DataStructs.h"

#include <vector>

namespace tca {

  TCEvent evt;
  TCConfig tcc;
  std::vector<TjForecast> tjfs;
  ShowerTreeVars stv;
  // vector of hits, tjs, etc in each slice
  std::vector<TCSlice> slices;
  std::vector<TrajPoint> seeds;

  const std::vector<std::string> AlgBitNames {
    "MaskHits",
    "MaskBadTPs",
    "Michel",
    "DeltaRay",
    "CTKink",
    "CTStepChk",
    "RvPrp",
    "CHMUH",
    "Split",
    "Comp3DVx",
    "Comp3DVxIG",
    "HED",
    "HamBragg",
    "HamVx",
    "HamVx2",
    "JunkVx",
    "JunkTj",
    "Killed",
    "Merge",
    "LastEndMerge",
    "TEP",
    "CHMEH",
    "FillGap",
    "Ghost",
    "MrgGhost",
    "ChkInTraj",
    "StopBadFits",
    "FixBegin",
    "FTBChg",
    "BeginChg",
    "FixEnd",
    "BraggSplit",
    "UUH",
    "VtxTj",
    "ChkVxTj",
    "MisdVxTj",
    "Photon",
    "HaloTj",
    "NoFitToVx",
    "VxMerge",
    "VxNeutral",
    "NoKinkChk",
    "SoftKink",
    "ChkStop",
    "ChkStopEP",
    "ChkChgAsym",
    "FTBRvProp",
    "StopAtTj",
    "Mat3D",
    "TjHiVx3Score",
    "VtxHitsSwap",
    "SplitHiChgHits",
    "ShowerLike",
    "KillInShowerVx",
    "ShowerTj",
    "ShwrParent",
    "MergeOverlap",
    "MergeSubShowers",
    "MergeSubShowersTj",
    "MergeNrShowers",
    "MergeShChain",
    "CompleteShower",
    "SplitTjCVx",
    "MakePFPTjs",
    "FillGaps3D",
    "TEP3D"
  };
  
  const std::vector<std::string> EndFlagNames {
    "Signal",
    "AtKink",
    "AtVtx",
    "Bragg",
    "AtTj",
    "OutFV",
    "NoFitVx"
  };

  const std::vector<std::string> VtxBitNames {
    "VtxTrjTried",
    "Fixed",
    "OnDeadWire",
    "HiVx3Score",
    "VtxTruMatch",
    "VtxMerged",
    "VtxIndPlnNoChg"
  };

  const std::vector<std::string> StrategyBitNames {
    "Normal",
    "StiffEl",
    "StiffMu",
    "Slowing"
  };

  geo::PlaneID DecodeCTP(CTP_t CTP) {
    auto const cryo = (CTP / Cpad);
    return geo::PlaneID(
      /* Cryostat */ cryo,
           /* TPC */ (CTP - cryo * Cpad) / Tpad,
         /* Plane */ (CTP % 10)
      );
  }

} // namespace tca
