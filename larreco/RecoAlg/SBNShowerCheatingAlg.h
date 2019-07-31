#ifndef SBNShowerCheatingAlg_hxx
#define SBNShowerCheatingAlg_hxx

//Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "larreco/ShowerFinder/ShowerElementHolder.hh"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larreco/RecoAlg/SBNShowerAlg.h"
#include "larreco/RecoAlg/MCRecoUtils/RecoUtils.h"

//C++ Includes
#include <iostream>
#include <vector>
#include <map>

//Root Includes
#include "TVector3.h"
#include "TMath.h"
#include "TPrincipal.h"
#include "TVector.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TString.h"

namespace shower {
  class SBNShowerCheatingAlg;
}

class shower::SBNShowerCheatingAlg {
  public:
    SBNShowerCheatingAlg(const fhicl::ParameterSet& pset);

    std::map<int,const simb::MCParticle*> GetTrueParticleMap();
    std::map<int,std::vector<int> > GetTrueChain(std::map<int,const simb::MCParticle*> &trueParticles);
    void CheatDebugEVD(const simb::MCParticle* trueParticle, art::Event& Event,
        reco::shower::ShowerElementHolder& ShowerEleHolder,
        const art::Ptr<recob::PFParticle>& pfparticle);

  private:

    shower::SBNShowerAlg fSBNShowerAlg;

    art::InputTag                                       fHitModuleLabel;
    art::InputTag                                       fPFParticleModuleLabel;
    art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
    art::ServiceHandle<art::TFileService>   tfs;
};
#endif
