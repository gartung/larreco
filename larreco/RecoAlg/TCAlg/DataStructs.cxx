#include "larreco/RecoAlg/TCAlg/DataStructs.h"

namespace tca {
  const std::vector<std::string> AlgBitNames {
    "MaskHits",
    "CTKink",
    "CTStepChk",
    "TryNextPass",
    "RevProp",
    "CHMH",
    "SplitTraj",
    "Comp3DVx",
    "Comp3DVxIG",
    "HED",
    "HamVx",
    "HamVx2",
    "JunkTj",
    "Killed",
    "EndMerge",
    "TrimEndPts",
    "CHMEH",
    "FillGap",
    "Ghost",
    "ChkInTraj",
    "StopBadFits",
    "FixBegin",
    "FixEnd",
    "UseUnusedHits",
    "VtxTj",
    "RefVtx",
    "MBadTPs",
    "NoKinkChk",
    "SoftKink",
    "ChkStop",
    "ChkAllStop",
    "FTBRevProp",
    "StopAtTj",
    "Match3D",
    "InShower",
    "ShowerParent",
    "ShowerTj"
  };

  const std::vector<std::string> StopFlagNames {
    "Signal",
    "AtKink",
    "AtVtx",
    "Bragg",
    "RvPrp",
    "AtTj",
    "BadFits"
  };
  
  const std::vector<std::string> VtxBitNames {
    "Fixed",
    "VtxTrjTried",
    "OnDeadWire",
    "VtxRefined",
    "NiceVtx"
  } ;
} // namespace tca

