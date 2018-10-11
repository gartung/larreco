////////////////////////////////////////////////////////////////////////
//
//
// TCAlg vertex code
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERALGVERTEX_H
#define TRAJCLUSTERALGVERTEX_H


// C/C++ standard libraries
#include <array>
#include <vector>
#include <bitset>
#include <utility> // std::pair<>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"

namespace tca {

  extern TCEvent evt;
  extern TCConfig tcc;
//  extern HistStuff hist;
  // vector of hits, tjs, etc in each slice
  extern std::vector<TCSlice> slices;

  void MakeJunkVertices(TCSlice& slc, const CTP_t& inCTP);
  void Find2DVertices(TCSlice& slc, const CTP_t& inCTP);
  void FindNeutralVertices(TCSlice& slc);
  void MakeJunkTjVertices(TCSlice& slc, const CTP_t& inCTP);
  void ChkVxTjs(TCSlice& slc, const CTP_t& inCTP, bool prt);
  bool MergeWithVertex(TCSlice& slc, VtxStore& vx2, unsigned short existingVxID);
  void SplitTrajCrossingVertices(TCSlice& slc, CTP_t inCTP);
  void FindHammerVertices(TCSlice& slc, const CTP_t& inCTP);
  void FindHammerVertices2(TCSlice& slc, const CTP_t& inCTP);
  void Find3DVertices(TCSlice& slc);
  void Match3DVtxTjs(TCSlice& slc, bool prt);
  void CompleteIncomplete3DVertices(TCSlice& slc);
  bool RefineVtxPosition(TCSlice& slc, const Trajectory& tj, unsigned short& nearPt, short nPtsToChk, bool prt);
  void CompleteIncomplete3DVerticesInGaps(TCSlice& slc);
  // Improve hit assignments near vertex 
  void VtxHitsSwap(TCSlice& slc, const CTP_t inCTP);

  unsigned short TPNearVertex(TCSlice& slc, const TrajPoint& tp);
  bool AttachPFPToVertex(TCSlice& slc, PFPStruct& pfp, unsigned short end, unsigned short vx3ID, bool prt);
  bool AttachAnyTrajToVertex(TCSlice& slc, unsigned short iv, bool prt);
  bool AttachTrajToVertex(TCSlice& slc, Trajectory& tj, VtxStore& vx, bool prt);
  float TrajPointVertexPull(TCSlice& slc, const TrajPoint& tp, const VtxStore& vx);
  float VertexVertexPull(TCSlice& slc, const Vtx3Store& vx1, const Vtx3Store& vx2);
  float VertexVertexPull(TCSlice& slc, const VtxStore& vx1, const VtxStore& vx2);
  bool FitVertex(TCSlice& slc, VtxStore& vx, bool prt);
  bool FitVertex(TCSlice& slc, VtxStore& vx, std::vector<TrajPoint> vxTp, bool prt);
  bool StoreVertex(TCSlice& slc, VtxStore& vx);
  bool ChkVtxAssociations(TCSlice& slc, const CTP_t& inCTP);
  void ScoreVertices(TCSlice& slc);
  void KillPoorVertices(TCSlice& slc);
  void SetVx2Score(TCSlice& slc);
  void SetVx2Score(TCSlice& slc, VtxStore& vx2);
  void SetVx3Score(TCSlice& slc, Vtx3Store& vx3);
  unsigned short Vx3Topo(TCSlice& slc, Vtx3Store& vx3);
  void SetHighScoreBits(TCSlice& slc, Vtx3Store& vx3);
  bool MakeVertexObsolete(TCSlice& slc, VtxStore& vx2, bool forceKill);
  bool MakeVertexObsolete(TCSlice& slc, Vtx3Store& vx3);
  std::vector<int> GetVtxTjIDs(const TCSlice& slc, const VtxStore& vx2);
  std::vector<int> GetVtxTjIDs(const TCSlice& slc, const Vtx3Store& vx3, float& score);
  std::vector<unsigned short> GetPFPVertices(const TCSlice& slc, const PFPStruct& pfp);
  void PosInPlane(const TCSlice& slc, const Vtx3Store& vx3, unsigned short plane, Point2_t& pos);
  unsigned short IsCloseToVertex(TCSlice& slc, VtxStore& vx);
  unsigned short IsCloseToVertex(TCSlice& slc, Vtx3Store& vx3);
} // namespace

#endif // ifndef TRAJCLUSTERALGVERTEX_H
