////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgDetermineT0
// Author:      L. Whitehead (leigh.howard.whitehead@cern.ch),
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/PMAlgDetermineT0.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

pma::PMAlgDetermineT0::PMAlgDetermineT0(const pma::PMAlgDetermineT0::Config &conf): fDriftWidthMargin(conf.DriftWidthMargin())
{
  fDriftMin.clear();
  fDriftMax.clear();
  fDriftCoord = -999;

  GetTPCDriftWidths();
}

void pma::PMAlgDetermineT0::DetermineT0(pma::TrkCandidateColl& tracks)
{
  
  size_t n = 0;

  // First method is to get T0 for those tracks that cross an entire drift length
  n += GetAPAToCPAT0(tracks);

  mf::LogInfo("pma::PMAlgDetermineT0") << " Calculated T0 for " << n << " tracks." << std::endl;

}

size_t pma::PMAlgDetermineT0::GetAPAToCPAT0(pma::TrkCandidateColl& tracks)
{
  // Remember, tracks actually crossing the CPA or APA have a T0 measured
  // as part of the track stitching process

  size_t n = 0;

  for (auto & t : tracks.tracks())
  {

    // We need to ensure the first and last node of the track are in the same main drift volume
    pma::Node3D* firstNode = t.Track()->FirstElement();
    pma::Node3D* finalNode = t.Track()->LastElement();

    // If one TPC is -1, who knows what to do?
    if(firstNode->TPC() == -1 || finalNode->TPC() == -1) continue;

    // If the minimum drift coordinate is the same in both TPCs should be in same drift volume
    // This ensures the track did not cross the cathode or anode
    if(fDriftMin[firstNode->TPC()] == fDriftMin[finalNode->TPC()])
    {
      // Only need one TPC as we just need to get the drift coordinates for this drift volume
      unsigned int tpc = firstNode->TPC();
    
      // Get the track length and drift length of the TPC
      double driftTrkLen = fabs(firstNode->Point3D()[fDriftCoord] - finalNode->Point3D()[fDriftCoord]);
      double driftTPCLen = fabs(fDriftMax[tpc] - fDriftMin[tpc]);

      if(fabs(driftTrkLen - driftTPCLen) < fDriftWidthMargin)
      {
        // This track should have a T0 assigned, so let's do it
        double minDrift = (firstNode->Point3D()[fDriftCoord] < finalNode->Point3D()[fDriftCoord]) ? 
                          firstNode->Point3D()[fDriftCoord] : finalNode->Point3D()[fDriftCoord];
        double maxDrift = (firstNode->Point3D()[fDriftCoord] >= finalNode->Point3D()[fDriftCoord]) ? 
                          firstNode->Point3D()[fDriftCoord] : finalNode->Point3D()[fDriftCoord];

        // Find which one is clostest to the cathode, then use that to tell us the shift
        double driftShift = 0;
        if(fDriftDir[tpc] > 0){
          // For positive drift, the cathode is always the minimum value
          driftShift = fDriftMin[tpc] - minDrift;
        }
        else{
          // For negative drift, the cathode is always the maximum value
          driftShift = fDriftMax[tpc] - maxDrift;
        }

        // Apply the drift direction shift and the track takes care of T0
        t.Track()->GetRoot()->ApplyDriftShiftInTree(driftShift);

        mf::LogInfo("pma::PMAlgDetermineT0") << "- Applied shift " << driftShift << " to a track in tpc " << tpc << std::endl;
        mf::LogInfo("pma::PMAlgDetermineT0") << " - Old x values = " << minDrift << " and " << maxDrift << std::endl;
        mf::LogInfo("pma::PMAlgDetermineT0") << " - New x values = " << minDrift + driftShift << " and " << maxDrift + driftShift << std::endl;
        mf::LogInfo("pma::PMAlgDetermineT0") << " - TPC edges    = " << fDriftMin[tpc] << " and " << fDriftMax[tpc] << std::endl;

        ++n;
      }

    } 
   
  } // End loop over tracks

  return n;

}

void pma::PMAlgDetermineT0::GetTPCDriftWidths(){

  // Get the geometry
  auto const* geom = lar::providerFrom<geo::Geometry>();

  // Since we can stack TPCs, we can't just use geom::TPCGeom::Height()
  for (geo::TPCID const& tID: geom->IterateTPCIDs()) {
    geo::TPCGeo const& TPC = geom->TPC(tID);

    // We need to make sure we only include the real TPCs
    // We have dummy TPCs in the protoDUNE and DUNE geometries
    // The dummy ones have a drift distance of only ~13 cm.
    if(TPC.DriftDistance() < 25.0){
      continue;
    }

    if(fDriftCoord == -999){
      // Note following functions returns: +/- 1, 2, 3 for +/- x, y, z and 0 for unknown
      short int tempCoord = TPC.DetectDriftDirection();
      if(tempCoord != 0){
        // Convert to array index
        short int sign = 1;
        if(tempCoord < 0) sign = -1;
        fDriftDir.insert(std::make_pair(tID.TPC,sign));
        fDriftCoord = abs(tempCoord) - 1; // x = 0, y = 1, z = 2;
      }
    }
    if(fDriftCoord == -999) continue;

    // get centre in world coordinates
    double origin[3] = {0.};
    double centre[3] = {0.};
    TPC.LocalToWorld(origin, centre);
    double tpcDim[3] = {TPC.HalfWidth(), TPC.HalfHeight(), 0.5*TPC.Length() };
    
    fDriftMin.insert(std::make_pair(tID.TPC,centre[fDriftCoord] - tpcDim[fDriftCoord]));
    fDriftMax.insert(std::make_pair(tID.TPC,centre[fDriftCoord] + tpcDim[fDriftCoord]));
  }

}


