#ifndef FDSELECTIONUTILS_H_SEEN
#define FDSELECTIONUTILS_H_SEEN


///////////////////////////////////////////////
// RecoUtils.h
//
// A few reco utilities like truth matching 
// D Brailsford (adapted from work by D Brailsford and M Wallbank), October 2017
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
//#include "lardataobj/RecoBase/Track.h"
//#include "lardataobj/RecoBase/Shower.h"
//#include "lardataobj/AnalysisBase/MVAPIDResult.h"
//#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

// c++
#include <vector>
#include <map>

// ROOT
#include "TTree.h"


namespace RecoUtils{
  int TrueParticleID(const art::Ptr<recob::Hit>& hit); //Returns the geant4 ID which contributes the most to a single reco hit.  The matching method looks for true particle which deposits the most true energy in the reco hit
  int TrueParticleIDFromTotalTrueEnergy(const std::vector<art::Ptr<recob::Hit> >& hits); //Returns the geant4 ID which contributes the most to the vector of hits.  The matching method looks for which true particle deposits the most true energy in the reco hits
  int TrueParticleIDFromTotalRecoCharge(const std::vector<art::Ptr<recob::Hit> >& hits);  //Returns the geant4 ID which contributes the most to the vector of hits.  The matching method looks for which true particle contributes the most reconstructed charge to the hit selection (the reco charge of each hit is correlated with each maximally contributing true particle and summed)
  int TrueParticleIDFromTotalRecoHits(const std::vector<art::Ptr<recob::Hit> >& hits);  //Returns the geant4 ID which contributes the most to the vector of hits.  The matching method looks for which true particle maximally contributes to the most reco hits
  bool IsInsideTPC(TVector3 position, double distance_buffer); //Checks if a position is within any of the TPCs in the geometry (user can define some distance buffer from the TPC walls)

  float TrueEnergyDepositedFromMCTrack(int TrackID,const std::vector<art::Ptr<sim::SimChannel> >& simchannels); //Returns the total energy deposited from the track id given.
  std::map<geo::PlaneID,int> NumberofMCWiresHit(int TrackID,const std::vector<art::Ptr<sim::SimChannel> > & simchannels); // Returns the number of Wires that saw an energy deposit in Monte Carlo from a track.Might be useful to add an energy cut on this. 

  std::map<geo::PlaneID,int> NumberofHitsThatContainEnergyDepositedByTrack(int TrackID, const std::vector<art::Ptr<recob::Hit> >& hits); //Returns the number of hits in the reconstruction that saw an energy deposition by the a track. Might be useful to add an energy cut on this. 

  float TotalEnergyDepinHits(const std::vector<art::Ptr<recob::Hit> >& hits); //Returns the amount of energy deposited in the detector (before recombination and lifetime effects) in the hits. 

  float TotalEnergyDepinHitsFromTrack(const std::vector<art::Ptr<recob::Hit> >& hits, int TrackID); //Returns the amount of energy deposited in the detector (before recombination and lifetime effects)in the hits from a given particle. 

  std::map<int,float> TrueEnergyDepositedFromMCTracks(const std::vector<art::Ptr<sim::SimChannel> >& simchannels); //Returns a map of the Energies depsoited by each track id.

  std::map<int,std::map<geo::PlaneID,int> > NumberofMCWiresHitMap(const std::vector<art::Ptr<sim::SimChannel> >& simchannels); //Returns of map of the number of wires hit in each plane for each track. 

  std::map<geo::PlaneID,int> NumberofHitsThatContainEnergyDepositedByTracks(std::vector<int> TrackIDs, const std::vector<art::Ptr<recob::Hit> >& hits);//Number of Hits that containt energy from the list of tracks. 

}

#endif
