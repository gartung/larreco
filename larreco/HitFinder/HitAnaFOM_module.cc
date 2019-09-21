////////////////////////////////////////////////////////////////////////
// Class:       HitAnaFOM
// Plugin Type: analyzer (art v3_02_06)
// File:        HitAnaFOM_module.cc
//
// This module calculates a Figure Of Merit (FOM) for a collection of hits
// reconstructed with a hit finder FCL configuration set labeled by
// the FclSet string. FclSet and the FOM are written to the text file
// ResultsFile. The FclSet variable is defined in the Python script fclScan.py
// using a format (for example) <FCL configuration group>:<FCL file name>. The <FCL file name>
// has the format <FCL file template name>_<first variable variant><second variable variant>.fcl
//
// Example: FclSet = "GRPChi:hitopt_42.fcl" refers to a group 
// the FOM produced here is a variant of the template FCL file hitopt.fcl in which
// the user defined a group of two fcl variables with the label GRPChi that are pressumed
// to be correlated. The "4" refers to the fourth step of the first fcl variable and "2" refers
// to the second step of the second variable.
//
// The ResultsFile variable is specified by the user in the template FCL file.
// The FOM calculated and the FclSet variable are written to the file with the
// expectation that the user will select the FCL configuration that has the "best" FOM
// and transfer that configuration to the production FCL file.
//
// Generated at Tue Sep  3 15:51:38 2019 by Bruce Baller (baller@fnal.gov) using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include <fstream>      // std::fstream
#include <ctime>        // std::time

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nug4/ParticleNavigation/ParticleList.h"

#include "TH1.h"


namespace hit {
  class HitAnaFOM;
}

class hit::HitAnaFOM : public art::EDAnalyzer {
public:
  explicit HitAnaFOM(fhicl::ParameterSet const& p);

   HitAnaFOM(HitAnaFOM const&) = delete;
  HitAnaFOM(HitAnaFOM&&) = delete;
  HitAnaFOM& operator=(HitAnaFOM const&) = delete;
  HitAnaFOM& operator=(HitAnaFOM&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  art::InputTag fCalDataModuleLabel;
  art::InputTag fHitsModuleLabel;
  
  float fMinEnergy;
  float fUnresolvedEMFraction;
  
  std::string fResultsFile;
  std::string fFclSet;
  
  std::array<float, 3> nHitsExpected;
  std::array<float, 3> nEMDtr;
  std::array<float, 3> nHitsReconstructed;
  std::array<float, 3> nMatchedHitsReconstructed;
  
  float FOM(unsigned short plane);
  
  TH1F* fPeakAmplitude;
  TH1F* fPeakAmplitude_MCMat;
  TH1F* fMultiplicity;
  TH1F* fRMS;

};

///////////////
hit::HitAnaFOM::HitAnaFOM(fhicl::ParameterSet const& pset) : EDAnalyzer{pset} {
  fCalDataModuleLabel  = pset.get< std::string >("CalDataModuleLabel");
  fHitsModuleLabel  = pset.get<art::InputTag>("HitsModuleLabel", "NA");
  fMinEnergy = pset.get<float>("MinEnergy" , 0.0);
  fUnresolvedEMFraction = pset.get<float>("UnresolvedEMFraction" , 0.2);
  fResultsFile = pset.get<std::string>("ResultsFile");
  fFclSet = pset.get<std::string>("FclSet");
}

///////////////
void hit::HitAnaFOM::analyze(art::Event const& evt) {
  
  auto inputHits = art::Handle<std::vector<recob::Hit>>();
  if(!evt.getByLabel(fHitsModuleLabel, inputHits)) throw cet::exception("HitAnaFOM")<<"Failed to get a handle to hit collection  '"<<fHitsModuleLabel.label()<<"'\n";
  int nInputHits = (*inputHits).size();
  std::cout<<"Found "<<nInputHits<<" input hits\n";
  
  auto wires = art::Handle<std::vector<recob::Wire>>();
  if(!evt.getByLabel(fCalDataModuleLabel, wires)) throw cet::exception("HitAnaFOM")<<"Failed to get a handle to wires from  '"<<fCalDataModuleLabel.label()<<"'\n";

  art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
  art::ServiceHandle<geo::Geometry const> geom;
  
  // find all GEANT track IDEs on all wires in the full readout time. 
  // Ideally we expect one reconstructed hit for all ionizations associated with 
  // one MCParticle on a wire. This assumption will be violated in
  // the following cases:
  //   1) ionization from several MCParticles overlap in a view (unavoidable),
  //   2) delta-rays and other EM daughters overlap in time (corrected),
  //   3) MCParticles travel along the drift direction are reconstructed as many long-pulse hits (rare)
  //   4) below-threshold ioniziation (< 1/3 MIP)
  
  double startTime = 0.;
  // find the max time from the hit collection. This assumes that the
  // hit reconstruction was reasonably efficient. Also count the total number
  // of reconstructed hits
  double endTime = 0.;
  for(unsigned int iht = 0; iht < (*inputHits).size(); ++iht) {
    auto& hit = (*inputHits)[iht];
    double hitMaxTime = hit.PeakTime() + 2.5 * hit.RMS();
    if(hitMaxTime > endTime) endTime = hitMaxTime;
    ++nHitsReconstructed[hit.WireID().Plane];
    fPeakAmplitude->Fill(hit.PeakAmplitude());
    fMultiplicity->Fill(hit.Multiplicity());
    fRMS->Fill(hit.RMS());
  } // iht
  
  for(unsigned int iht = 0; iht < (*inputHits).size(); ++iht) {
    auto& hit = (*inputHits)[iht];
    unsigned int plane = hit.WireID().Plane;
    double hitMinTime = hit.PeakTime() - 2.5 * hit.RMS();
    double hitMaxTime = hit.PeakTime() + 2.5 * hit.RMS();
    auto tides = bt_serv->ChannelToTrackIDEs(hit.Channel(), hitMinTime, hitMaxTime);
    // find a valid GEANT TrackID if the hit is matched. Iterate through all of the
    // IDEs to find the one with the highest energy contribution to the hit
    int trackID = -1;
    float maxE = fMinEnergy;
    for(auto& tide : tides) {
      if(tide.energy < maxE) continue;
      maxE = tide.energy;
      trackID = abs(tide.trackID);
    }
    // ignore this hit if it isn't matched
    if(trackID < 0) continue;
    fPeakAmplitude_MCMat->Fill(hit.PeakAmplitude());
    // get all the energy deposited on this wire by this GEANT track
    tides = bt_serv->ChannelToTrackIDEs(hit.Channel(), startTime, endTime);
    for(auto& tide : tides) {
      if(tide.trackID != trackID) continue;
      if(maxE > tide.energy) maxE = tide.energy;
      nMatchedHitsReconstructed[plane] += maxE / tide.energy;
      break;
    } // itde
  } // iht
  
  // Finally we count the number of track IDEs on each wire within the readout time window
  // with the (unrealized) expectation that there shoule be one reconstructed hit for
  // each IDE
  for(unsigned int iw = 0; iw < (*wires).size(); ++iw) {
    auto& wire = (*wires)[iw];
    raw::ChannelID_t channel = wire.Channel();
    std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
    unsigned short plane = wids[0].Plane;
    auto tides = bt_serv->ChannelToTrackIDEs(channel, startTime, endTime);
    for(auto& tide : tides) {
      // don't expect a hit if the energy deposition is below threshold
      if(tide.energy < fMinEnergy) continue;
      // deal with EM daughters. 
      if(tide.trackID < 0) ++nEMDtr[plane];
      // Here we can just count
      ++nHitsExpected[plane];
    } // tide
  } // iw
  
} // analyze

///////////////
float hit::HitAnaFOM::FOM(unsigned short plane) {
  if(plane > 2) return -1;
  if(nHitsExpected[plane] < 10 || nHitsReconstructed[plane] < 10) return -1;
  float nhExpect = nHitsExpected[plane] - fUnresolvedEMFraction * nEMDtr[plane];
  float eff = nMatchedHitsReconstructed[plane] / nhExpect;
  float pur = nMatchedHitsReconstructed[plane] / nHitsReconstructed[plane];
  return eff * pur;
} // fom


///////////////
void hit::HitAnaFOM::beginJob() {
  for(unsigned short plane = 0; plane < 3; ++plane) {
    nHitsExpected[plane] = 0;
    nEMDtr[plane] = 0;
    nHitsReconstructed[plane] = 0;
    nMatchedHitsReconstructed[plane] = 0;
  } // plane
  // define some histograms
  art::ServiceHandle<art::TFileService const> tfs;
  fPeakAmplitude = tfs->make<TH1F>("fPeakAmplitude", "Hit Peak Amplitude", 1000, 0, 200);
  fPeakAmplitude_MCMat = tfs->make<TH1F>("fPeakAmplitude_MCMat", "Hit Peak Amplitude - MC Matched", 1000, 0, 200);
  fMultiplicity = tfs->make<TH1F>("fMultiplicity", "Hit Multiplicity", 50, 0, 50);
  fRMS = tfs->make<TH1F>("fRMS", "Hit RMS", 50, 0, 50);
  
}  // beginJob

///////////////
void hit::HitAnaFOM::endJob() {
  
  for(unsigned short plane = 0; plane < 3; ++plane) {
    if(nHitsExpected[plane] < 10) continue;
    std::cout<<"Plane "<<plane<<" expected "<<(int)nHitsExpected[plane]<<" nEMDtr "<<(int)nEMDtr[plane];
    std::cout<<" reco "<<(int)nHitsReconstructed[plane]<<" matched "<<(int)nMatchedHitsReconstructed[plane];
    std::cout<<" fom "<<std::setprecision(3)<<FOM(plane);
    std::cout<<"\n";
  } // plane
  
  // see if the file exists by trying to open it
  std::fstream test(fResultsFile);
  bool fileExisted = true;
  if(!test) fileExisted = false;
  std::fstream results;
  // append to the results file
  results.open(fResultsFile, std::fstream::app);
  if(!fileExisted) {
    results<<"Results file format: <fclGroup you selected in fclScan.py>:<name of the fcl file> FOM <Figure of Merit for each plane> <Date/Time when this FOM was calculated>\n";
    results<<"FOM = (number of reconstructed hits in the plane) / (number of GEANT TrackIDEs with cuts/corrections).\n";
    results<<"See HitAnaFOM_Module.cc for details\n";
    results<<"-------------------------------------------------------------------------------------\n";
  }
  results<<fFclSet<<" FOM";
  for(unsigned short plane = 0; plane < 3; ++plane) {
    results<<std::setprecision(3)<<std::setw(6)<<FOM(plane);
  } // plane
  std::time_t stime = std::time(nullptr);
  results<<" "<<std::asctime(std::localtime(&stime));
  results.close();
} // endJob

DEFINE_ART_MODULE(hit::HitAnaFOM)
