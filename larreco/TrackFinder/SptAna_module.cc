
////////////////////////////////////////////////////////////////////////
//
// SptAna class - Check SpacePoint performance using MC-matched hits
//
// Bruce Baller
//
////////////////////////////////////////////////////////////////////////

//#ifndef SPTANA_H
//#define SPTANA_H


#include <iomanip>
#include <vector>
#include <string>
#include <array>
#include <fstream>

//Framework includes
#include "art/Framework/Core/ModuleMacros.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/EDAnalyzer.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

//#endif 


namespace spt {

  class SptAna : public art::EDAnalyzer {
  public:
    
    explicit SptAna(fhicl::ParameterSet const& pset); 
    virtual ~SptAna();
    
    /// read access to event
    void analyze(const art::Event& evt);
    void beginJob();
    void endJob();
  private:
    art::InputTag fHitModuleLabel;
    art::InputTag fHitTruthModuleLabel;
    art::InputTag fSpacePointModuleLabel;
    unsigned int fSptCnt {0};
    // Count of 4 possible SpacePoint -> Hit matches for a 3-plane TPC
    // 0 = not matched to any MCParticle
    // 1, 2, and 3 = matched to 1, 2 or 3 MCParticles. Ideally all hits would be matche
    std::array<unsigned int, 4> fMatchCnt {{0}};

    
  }; // class SpacePoint
  
  DEFINE_ART_MODULE(SptAna)

  //--------------------------------------------------------------------
  SptAna::SptAna(fhicl::ParameterSet const& pset)  
  : EDAnalyzer(pset)
  , fHitModuleLabel        (pset.get< art::InputTag > ("HitModuleLabel"))
  , fHitTruthModuleLabel    (pset.get< art::InputTag > ("HitTruthModuleLabel", "NA"))
  , fSpacePointModuleLabel  (pset.get< art::InputTag > ("SpacePointModuleLabel"))
  {
    
  }
  
  //------------------------------------------------------------------
  SptAna::~SptAna()
  {
    
  }
  
  //------------------------------------------------------------------
  void SptAna::beginJob()
  {
  } // beginJob
  
  //------------------------------------------------------------------
  void SptAna::endJob()
  {
    if(fSptCnt == 0) {
      std::cout<<"No triple-hit SpacePoints found\n";
      return;
    }
    std::cout<<"SptAna: SptCnt triple "<<fSptCnt<<" MatchCnt:";
    for(auto mc : fMatchCnt) std::cout<<" "<<mc;
    float eff = (float)fMatchCnt[1] / (float)fSptCnt;
    std::cout<<std::fixed<<std::setprecision(2)<<" eff "<<eff;
    std::cout<<"\n";
  } // endJob
  
  //------------------------------------------------------------------
  void SptAna::analyze(const art::Event& evt)
  {
    
    if(evt.isRealData()) return;
    // get a handle for the hit collection
    auto inputHits = art::Handle<std::vector<recob::Hit>>();
    if(!evt.getByLabel(fHitModuleLabel, inputHits)) throw cet::exception("SptAnaModule")<<"Failed to get a handle to hit collection '"<<fHitModuleLabel.label()<<"'\n";
    auto mcpHandle = art::Handle<std::vector<simb::MCParticle>>();
    if(!evt.getByLabel("largeant", mcpHandle)) throw cet::exception("SptAnaModule")<<"Failed to get a handle to MCParticles using 'largeant'\n";
    // create an index into the MCParticle vector for each matched hit
    std::vector<unsigned int> mcpIndex((*inputHits).size(), UINT_MAX);
    // BackTrackerHitMatchingData is supposedly available
    art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(inputHits, evt, fHitTruthModuleLabel);
    std::vector<art::Ptr<simb::MCParticle>> particle_vec;
    std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
    bool foundBTHMD = true;
    for(unsigned int iht = 0; iht < (*inputHits).size(); ++iht) {
      particle_vec.clear(); match_vec.clear();
      try{ particles_per_hit.get(iht, particle_vec, match_vec); }
      catch(...) {
        std::cout<<"BackTrackerHitMatchingData not found\n";
        foundBTHMD = false;
        break;
      }
      if(particle_vec.empty()) continue;
      int trackID = 0;
      for(unsigned short im = 0; im < match_vec.size(); ++im) {
        if(match_vec[im]->ideFraction < 0.5) continue;
        trackID = particle_vec[im]->TrackId();
        break;
      } // im
      if(trackID == 0) continue;
      // find the index
      for(unsigned int indx = 0; indx < (*mcpHandle).size(); ++indx) {
        auto& mcp = (*mcpHandle)[indx];
        if(mcp.TrackId() != trackID) continue;
        mcpIndex[iht] = indx;
        break;
      } // indx
    } // iht
    if(!foundBTHMD) {
      std::cout<<"try the old way\n";
    }
    // see how many are matched
    unsigned int nmt = 0;
    for(auto mcpi : mcpIndex) if(mcpi != UINT_MAX) ++nmt;
    std::cout<<"Found "<<nmt<<"/"<<(*inputHits).size()<<" Matched Hits / Total Hits\n";
    auto sptHandle = art::Handle<std::vector<recob::SpacePoint>>();
    if(!evt.getByLabel(fSpacePointModuleLabel, sptHandle)) throw cet::exception("SptAnaModule")<<"Failed to get a handle to SpacePoints\n";
    if((*sptHandle).empty()) return;
    art::FindManyP<recob::Hit> fh(sptHandle, evt, fSpacePointModuleLabel);
    for(unsigned int ispt = 0; ispt < (*sptHandle).size(); ++ispt) {
      //                    mcpIndex       count
      std::vector<std::pair<unsigned int, unsigned int>> mcpCnt;
      // number of hits associated with the Spt
      unsigned short nht = fh.at(ispt).size();
      if(nht != 3) {
//        std::cout<<"Ignoring Spt with "<<nht<<" hits. Ignoring all future occurrences...\n";
        continue;
      }
      ++fSptCnt;
      for(unsigned short ii = 0; ii < nht; ++ii) {
        auto& hit = fh.at(ispt).at(ii);
        unsigned int iht = hit.key();
        // require a MCParticle match
        if(mcpIndex[iht] == UINT_MAX) continue;
        // look for this MCParticle index in mcpCnt
        unsigned int mcpi = mcpIndex[iht];
        unsigned short indx = 0;
        for(indx = 0; indx < mcpCnt.size(); ++indx) if(mcpCnt[indx].first == mcpi) break;
        // not found so add it to the vector
        if(indx == mcpCnt.size()) mcpCnt.push_back(std::make_pair(mcpi, 0));
        ++mcpCnt[indx].second;
      } // ii
      ++fMatchCnt[mcpCnt.size()];
    } // ispt


  } // analyze

} // namespace spt


namespace spt{
  
  
} 
