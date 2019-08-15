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
#include "larreco/RecoAlg/TRACSAlg.h"

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
			  reco::shower::ShowerElementHolder& ShowerEleHolder
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
					    reco::shower::ShowerElementHolder& ShowerEleHolder){

    std::cout << "on shower track hit direction" << std::endl;

    //Check the Track Hits has been defined
    if(!ShowerEleHolder.CheckElement("InitialTrackHits")){
      mf::LogError("ShowerTrackHitDirection") << "Initial track hits not set"<< std::endl;
      return 0;
    }
    
    //Check the Track Hits has been defined
    if(!ShowerEleHolder.CheckElement("InitialTrack")){
      mf::LogError("ShowerTrackHitDirection") << "Initial track not set"<< std::endl;
      return 0;
    }
   
    //Check the start position is set.
    if(fUsePandoraVertex && !ShowerEleHolder.CheckElement("ShowerStartPosition")){
      mf::LogError("ShowerTrackHitDirection") << "Start position not set, returning "<< std::endl;
      return 0;
    }
    
    //Get the spacepoints associated to hits
    art::Handle<std::vector<recob::Hit> > hitHandle;
    if (!Event.getByLabel(fHitModuleLabel, hitHandle)){
      throw cet::exception("ShowerTrackHitDirection") << "Could not configure the spacepoint handle. Something is configured incorrectly. Stopping"; 
      return 1;
    } 

    //Get the spacepoint handle. We need to do this in 3D.
    art::FindManyP<recob::SpacePoint> fmsp(hitHandle, Event, fPFParticleModuleLabel);
    if(!fmsp.isValid()){
      throw cet::exception("ShowerTrackHitDirection") << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }

    //Get the initial track hits.
    std::vector<art::Ptr<recob::Hit> > InitialTrackHits;
    ShowerEleHolder.GetElement("InitialTrackHits",InitialTrackHits);

    //Get the initial track?
    recob::Track InitialTrack;
    ShowerEleHolder.GetElement("InitialTrack");

    //Get the start poistion 
    TVector3 StartPosition;
    if(fUsePandoraVertex){
      StartPosition = ShowerEleHolder.GetShowerStartPosition();
    }
    else{ 
      geo::Point_t Start_point = InitialTrack.Start();
      StartPosition = {Start_point.X(),Start_point.Y(),Start_point.Z()};
    }
    
    //Calculate the mean direction and the the standard deviation
    float sumX=0, sumX2=0;
    float sumY=0, sumY2=0;
    float sumZ=0, sumZ2=0;

    //Get the spacepoints associated to the track hit 
    std::vector<art::Ptr<recob::SpacePoint > > intitaltrack_sp;
    for(auto const hit: InitialTrackHits){  
      std::vector<art::Ptr<recob::SpacePoint > > sps = fmsp.at(hit.key());
      for(auto const sp: sps){
	intitaltrack_sp.push_back(sp);
	
	//Get the direction relative to the start positon 
	fTRACSAlg.SpacePointPosition(sp) - StartPosition
      }
    }

    //Calculate the mean direction and the the standard deviation
    float sumX=0, sumX2=0;
    float sumY=0, sumY2=0;
    float sumZ=0, sumZ2=0;
    for(unsigned int traj=0; traj< InitialTrack.NumberTrajectoryPoints(); ++traj){
      geo::Vector_t TrajPosition = (InitialTrack.LocationAtPoint(traj) - StartPosition).Unit();
      sumX += TrajPosition.X(); sumX2 += TrajPosition.X()*TrajPosition.X();
      sumY += TrajPosition.Y(); sumY2 += TrajPosition.Y()*TrajPosition.Y();
      sumZ += TrajPosition.Z(); sumZ2 += TrajPosition.Z()*TrajPosition.Z();
    }

    float NumTraj =  InitialTrack.NumberTrajectoryPoints();
    geo::Vector_t Mean = {sumX/NumTraj,sumY/NumTraj,sumZ/NumTraj};
    Mean = Mean.Unit();
      
    float RMSX = 999;
    float RMSY = 999;
    float RMSZ = 999;
    if(sumX2/NumTraj - ((sumX/NumTraj)*((sumX/NumTraj))) > 0){
      RMSX = TMath::Sqrt(sumX2/NumTraj - ((sumX/NumTraj)*((sumX/NumTraj))));
    }
    if(sumY2/NumTraj - ((sumY/NumTraj)*((sumY/NumTraj))) > 0){
      RMSY = TMath::Sqrt(sumY2/NumTraj - ((sumY/NumTraj)*((sumY/NumTraj))));
    }
    if(sumZ2/NumTraj - ((sumZ/NumTraj)*((sumZ/NumTraj))) > 0){
      RMSZ = TMath::Sqrt(sumZ2/NumTraj - ((sumZ/NumTraj)*((sumZ/NumTraj))));
    }

    
    

      
    //Loop over the hits
    for (unsigned int hitIt=0; hitIt<InitialTrackHits.size(); ++hitIt) {
      art::Ptr<recob::Hit> hit = InitialTrackHits.at(hitIt);
      
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
      ShowerEleHolder.SetShowerDirection(Direction);
      return 0;
    }

    return 0;
  }
}

  
DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackHitDirection)
  
