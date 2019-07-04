//############################################################################
//### Name:        ShowerTrackHitDirection                                    ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the shower direction using the         ###
//###              initial track method.                                   ###
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
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/RecoAlg/SBNShowerAlg.h"

//C++ Includes 
#include <iostream>
#include <cmath>

//Root Includes 
#include "TVector3.h"
#include "TMath.h"
#include "TH1.h" 

namespace ShowerRecoTools {

  
  class ShowerTrackHitDirection:IShowerTool {
    
  public:
    
    ShowerTrackHitDirection(const fhicl::ParameterSet& pset);
    
    ~ShowerTrackHitDirection(); 
    
    //Generic Direction Finder
    int CalculateProperty(const art::Ptr<recob::PFParticle>& pfparticle,
			  art::Event& Event,
			  reco::shower::ShowerPropertyHolder& ShowerPropHolder
			  ) override;

  private:
    
    // Define standard art tool interface
    void configure(const fhicl::ParameterSet& pset) override;

    //fcl
    bool fUsePandoraVertex;
    bool fUseStartPosition;
    art::InputTag fHitModuleLabel;
    art::InputTag fPFParticleModuleLabel;

    //services
    detinfo::DetectorProperties const* fDetProp;
  };
  
  
  ShowerTrackHitDirection::ShowerTrackHitDirection(const fhicl::ParameterSet& pset)
    : fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())

  {
    configure(pset);
  }
  
  ShowerTrackHitDirection::~ShowerTrackHitDirection()
  {
  }
  
  void ShowerTrackHitDirection::configure(const fhicl::ParameterSet& pset)
  {
    fUsePandoraVertex       = pset.get<bool>         ("UsePandoraVertex");
    fUseStartPosition       = pset.get<bool>         ("UseStartPosition");
    fHitModuleLabel         = pset.get<art::InputTag>("HitModuleLabel");     
    fPFParticleModuleLabel  = pset.get<art::InputTag>("PFParticleModuleLabel","");
  }

  int ShowerTrackHitDirection::CalculateProperty(const art::Ptr<recob::PFParticle>& pfparticle,
					    art::Event& Event,
					    reco::shower::ShowerPropertyHolder& ShowerPropHolder){

    std::cout << "on shower track hit direction" << std::endl;

    //Check the Track Hits has been defined
    if(!ShowerPropHolder.CheckInitialTrackHits()){
      mf::LogError("ShowerDirection") << "Initial track hits not set"<< std::endl;
      return 0;
    }
    
    //Check the Track Hits has been defined
    if(!ShowerPropHolder.CheckInitialTrack()){
      mf::LogError("ShowerDirection") << "Initial track not set"<< std::endl;
      return 0;
    }
   
    //Check the start position is set.
    if(fUsePandoraVertex && !ShowerPropHolder.CheckShowerStartPosition()){
      mf::LogError("ShowerTrackFinder") << "Start position not set, returning "<< std::endl;
      return 0;
    }
    
    //Get the spacepoints associated to hits
    art::Handle<std::vector<recob::Hit> > hitHandle;
    if (!Event.getByLabel(fHitModuleLabel, hitHandle)){
      throw cet::exception("ShowerStartPosition") << "Could not configure the spacepoint handle. Something is configured incorrectly. Stopping"; 
      return 1;
    } 

    art::FindManyP<recob::SpacePoint> fmsp(hitHandle, Event, fPFParticleModuleLabel);
    if(!fmsp.isValid()){
      throw cet::exception("ShowerStartPosition") << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }
    std::cout << "fmsp.size(): " << fmsp.size() << std::endl;
     

    std::vector<art::Ptr<recob::Hit> > InitialTrackHits = ShowerPropHolder.GetInitialTrackHits();

    recob::Track InitialTrack = ShowerPropHolder.GetInitialTrack();
    
    TH1D* XHist = new TH1D("XHist","XHist",30,-1,1);
    TH1D* YHist = new TH1D("YHist","YHist",30,-1,1);
    TH1D* ZHist = new TH1D("ZHist","ZHist",30,-1,1);

    if(fUseStartPosition){ 
      geo::Point_t StartPosition;
      if(fUsePandoraVertex){
	TVector3 StartPosition_vec = ShowerPropHolder.GetShowerStartPosition();
	StartPosition.SetCoordinates(StartPosition_vec.X(),StartPosition_vec.Y(),StartPosition_vec.Z()); 
      }
      else{ 
	StartPosition = InitialTrack.Start();
      }
      
      //Loop over the hits
      for (unsigned int hitIt=0; hitIt<InitialTrackHits.size(); ++hitIt) {
	art::Ptr<recob::Hit> hit = InitialTrackHits.at(hitIt);
	std::cout << "peaktime: " << hit->PeakTime() << std::endl;
      }

      geo::Vector_t Mean = {XHist->GetMean(), YHist->GetMean(), ZHist->GetMean()};
      geo::Vector_t RMS  = {XHist->GetRMS(), YHist->GetRMS(), ZHist->GetRMS()};

      //Remove trajectory points from the mean that are not with one sigma.
      float N = 0.;
      TVector3 Direction_Mean = {0,0,0};
      for(unsigned int traj=0; traj<InitialTrack.NumberTrajectoryPoints(); ++traj){
	geo::Vector_t Direction = InitialTrack.DirectionAtPoint(traj);
	if((TMath::Abs((Direction-Mean).X()) < 2*RMS.X()) && 
	   (TMath::Abs((Direction-Mean).Y())< 2*RMS.Y()) && 
	   (TMath::Abs((Direction-Mean).Z()) < 2*RMS.Z())){
	  TVector3 Direction_vec = {Direction.X(),Direction.Y(),Direction.Z()};
	  Direction_Mean += Direction_vec;
	  ++N;
	}
      }
      
      //Take the mean value
      TVector3 Direction = Direction_Mean*(1/N);
      ShowerPropHolder.SetShowerDirection(Direction);
      return 0;
    }

    return 0;
  }
}

  
DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackHitDirection)
  
