//############################################################################
//### Name:        ShowerTrackFinderCheater                                ###
//### Author:      Ed Tyley                                                ###
//### Date:        16.07.19                                                ###
//### Description: Cheating tool using truth for shower direction          ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/RecoAlg/SBNShowerAlg.h"
#include "larreco/RecoAlg/SBNShowerCheatingAlg.h"
#include "larreco/RecoAlg/MCRecoUtils/RecoUtils.h"
#include "larreco/RecoAlg/MCRecoUtils/ShowerUtils.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "nusimdata/SimulationBase/MCTruth.h"

//C++ Includes
#include <iostream>
#include <cmath>

//Root Includes
#include "TVector3.h"
#include "TMath.h"
#include "TPrincipal.h"
#include "TVector.h"
#include "TTree.h"

namespace ShowerRecoTools {

  class ShowerTrackFinderCheater:IShowerTool {

    public:

      ShowerTrackFinderCheater(const fhicl::ParameterSet& pset);

      ~ShowerTrackFinderCheater();

      //Generic Direction Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      //Algorithm functions
      shower::SBNShowerAlg         fSBNShowerAlg;
      shower::SBNShowerCheatingAlg fSBNShowerCheatingAlg;


      //fcl
      bool fDebugEVD;
      art::InputTag fPFParticleModuleLabel;
      art::InputTag fHitModuleLabel;
  };


  ShowerTrackFinderCheater::ShowerTrackFinderCheater(const fhicl::ParameterSet& pset)
    : fSBNShowerAlg(pset.get<fhicl::ParameterSet>("SBNShowerAlg")),
    fSBNShowerCheatingAlg(pset.get<fhicl::ParameterSet>("SBNShowerCheatingAlg"))
  {
    fDebugEVD              = pset.get<bool>          ("DebugEVD");
    fHitModuleLabel        = pset.get<art::InputTag> ("HitModuleLabel");
    fPFParticleModuleLabel = pset.get<art::InputTag> ("PFParticleModuleLabel","");
  }

  ShowerTrackFinderCheater::~ShowerTrackFinderCheater()
  {
  }

  int ShowerTrackFinderCheater::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle, art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){

    const simb::MCParticle* trueParticle;

    //Get the hits from the shower:
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("ShowerLinearEnergy") << "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    if (ShowerEleHolder.CheckElement("TrueParticle")){
      ShowerEleHolder.GetElement("TrueParticle",trueParticle);
    } else {

      //Could store these in the shower element holder and just calculate once?
      std::map<int,const simb::MCParticle*> trueParticles = fSBNShowerCheatingAlg.GetTrueParticleMap();
      std::map<int,std::vector<int> > showersMothers = fSBNShowerCheatingAlg.GetTrueChain(trueParticles);


      //Get the clusters
      art::Handle<std::vector<recob::Cluster> > clusHandle;
      if (!Event.getByLabel(fPFParticleModuleLabel, clusHandle)){
        throw cet::exception("ShowerLinearEnergy") << "Could not get the pandora clusters. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
        return 1;
      }
      art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleModuleLabel);
      std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(pfparticle.key());

      //Get the hit association
      art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleModuleLabel);

      std::vector<art::Ptr<recob::Hit> > showerHits;
      for(auto const& cluster: clusters){

        //Get the hits
        std::vector<art::Ptr<recob::Hit> > hits = fmhc.at(cluster.key());
        showerHits.insert(showerHits.end(),hits.begin(),hits.end());
      }

      //Get the true particle from the shower
      std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain(showersMothers,showerHits,2);

      if(ShowerTrackInfo.first==-99999) {
        mf::LogWarning("ShowerStartPosition") << "True Shower Not Found";
        return 1;
      }
      trueParticle = trueParticles[ShowerTrackInfo.first];
    }

    //This is all based on the shower vertex being known. If it is not lets not do the track
    if(!ShowerEleHolder.CheckElement("ShowerStartPosition")){
      mf::LogError("ShowerTrackFinder") << "Start position not set, returning "<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement("ShowerDirection")){
      mf::LogError("ShowerTrackFinder") << "Direction not set, returning "<< std::endl;
      return 1;
    }

    TVector3 ShowerStartPosition = {-999,-999,-999};
    ShowerEleHolder.GetElement("ShowerStartPosition",ShowerStartPosition);

    TVector3 ShowerDirection     = {-999,-999,-999};
    ShowerEleHolder.GetElement("ShowerDirection",ShowerDirection);

    //loop over all of the hits in the event and get the ones from the trueParticle
    //FInd the spacepoints assocsiated to these hits
    //Pass these space points to the FindTrackSpacePoints
    art::Handle<std::vector<recob::Hit> > hitHandle;
    std::vector<art::Ptr<recob::Hit> > hits;
    if(Event.getByLabel(fHitModuleLabel, hitHandle)){
      art::fill_ptr_vector(hits, hitHandle);
    }

    // Get the hits associated with the space points
    art::FindManyP<recob::SpacePoint> fmsph(hitHandle, Event, fPFParticleModuleLabel);
    if(!fmsph.isValid()){
      throw cet::exception("ShowerTrackFinderCheater") << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }

    std::map< art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit> > spacePointHitMap;
    //Get the hits from the true particle
    for (auto hit : hits){
      int trueParticleID = RecoUtils::TrueParticleID(hit);
      if (trueParticleID == trueParticle->TrackId()){
        std::vector<art::Ptr<recob::SpacePoint> > sps = fmsph.at(hit.key());
        if (sps.size() == 1){
          art::Ptr<recob::SpacePoint> sp = sps.front();
          spacePointHitMap[sp] = hit;
        }
      }
    }

    std::vector<art::Ptr<recob::Hit> > trackHits;
    std::vector<art::Ptr<recob::SpacePoint> > trackSpacePoints;

    for (auto spacePoint : trackSpacePoints){
      art::Ptr<recob::Hit> hit = spacePointHitMap[spacePoint];
      trackHits.push_back(hit);
    }

    ShowerEleHolder.SetElement(trackHits, "InitialTrackHits");
    ShowerEleHolder.SetElement(trackSpacePoints,"InitialTrackSpacePoints");

    if (fDebugEVD){
      fSBNShowerCheatingAlg.CheatDebugEVD(trueParticle, Event, ShowerEleHolder, pfparticle);
    }

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackFinderCheater)
