//############################################################################
//### Name:        ShowerSmartTrackTrajectoryPointDirection                     ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the shower direction using the         ###
//###              first trajectory of the initial track                   ###
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

#include "TFile.h"

namespace ShowerRecoTools {

  
  class ShowerSmartTrackTrajectoryPointDirection:IShowerTool {
    
  public:
    
    ShowerSmartTrackTrajectoryPointDirection(const fhicl::ParameterSet& pset);
    
    ~ShowerSmartTrackTrajectoryPointDirection(); 
    
    //Generic Direction Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
			 art::Event& Event,
			 reco::shower::ShowerElementHolder& ShowerEleHolder
			 ) override;
    
  private:
    
    // Define standard art tool interface
    void configure(const fhicl::ParameterSet& pset) override;

    //fcl
    bool  fUsePandoraVertex;
    bool  fDebugEVD; 
    bool  fAllowDynamicSliding;
    bool  fUsePositionInfo;
    bool  fUseStartPos;
    float fAngleCut;
  };
  
  
  ShowerSmartTrackTrajectoryPointDirection::ShowerSmartTrackTrajectoryPointDirection(const fhicl::ParameterSet& pset)
  {
    configure(pset);
  }
  
  ShowerSmartTrackTrajectoryPointDirection::~ShowerSmartTrackTrajectoryPointDirection()
  {
  }
  
  void ShowerSmartTrackTrajectoryPointDirection::configure(const fhicl::ParameterSet& pset)
  {
    fUsePandoraVertex         = pset.get<bool>("UsePandoraVertex");
    fDebugEVD                 = pset.get<bool>("DebugEVD");  
    fAllowDynamicSliding      = pset.get<bool>("AllowDynamicSliding"); 
    fUseStartPos              = pset.get<bool>("UseStartPos");
    fUsePositionInfo          = pset.get<bool>("UsePositionInfo");
    fAngleCut                 = pset.get<float>("AngleCut");
  }

  int ShowerSmartTrackTrajectoryPointDirection::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
					     art::Event& Event,
					     reco::shower::ShowerElementHolder& ShowerEleHolder){

    std::cout << "ShowerSmartTrackTrajectoryPointDirection_tool" << std::endl;

    //Check the Track has been defined
    if(!ShowerEleHolder.CheckElement("InitialTrack")){
      mf::LogError("ShowerDirection") << "Initial track not set"<< std::endl;
      return 1;
    }
    recob::Track InitialTrack;
    ShowerEleHolder.GetElement("InitialTrack",InitialTrack);

    //Smartly choose the which trajectory point to look at by ignoring the smush of hits at the vertex. 
    if(InitialTrack.NumberTrajectoryPoints() == 1){
      mf::LogError("ShowerDirection") << "Not Enough trajectory points."<< std::endl;
      return 1;
    }

    //Trajectory point which the direction is calcualted for.
    int trajpoint = 0;
    geo::Vector_t Direction_vec;
    
    if(fUsePositionInfo){
      //Get the start position.
      geo::Point_t StartPosition;
        
      if(fUsePandoraVertex){
	//Check the Track has been defined
	if(!ShowerEleHolder.CheckElement("ShowerStartPosition")){
	  mf::LogError("ShowerDirection") << "Shower start position not set"<< std::endl;
	  return 1;
	}
	TVector3 StartPosition_vec = {-999,-999,-999};
	ShowerEleHolder.GetElement("ShowerStartPosition",StartPosition_vec);
	StartPosition.SetCoordinates(StartPosition_vec.X(),StartPosition_vec.Y(),StartPosition_vec.Z()); 
      }
      else{ 
	StartPosition = InitialTrack.Start();
      }

      //Loop over the trajectory points and find two corresponding trajectory points where the angle between themselves (or themsleves and the start position) is less the fMinAngle.
      for(unsigned int traj=0; traj<InitialTrack.NumberTrajectoryPoints()-2; ++traj){
	++trajpoint;
      
	geo::Vector_t TrajPosition_vec      = InitialTrack.LocationAtPoint(traj) - StartPosition ;
	geo::Vector_t NextTrajPosition_vec;
	geo::Vector_t NextNextTrajPosition_vec;
	if(fUseStartPos){
	  NextTrajPosition_vec  = InitialTrack.LocationAtPoint(traj+1) - StartPosition;
	  NextNextTrajPosition_vec  = InitialTrack.LocationAtPoint(traj+2) - StartPosition;
	}
	else{
	  NextTrajPosition_vec  = InitialTrack.LocationAtPoint(traj+1)  - InitialTrack.LocationAtPoint(traj); 
	  NextNextTrajPosition_vec  = InitialTrack.LocationAtPoint(traj+2) - InitialTrack.LocationAtPoint(traj+1);
	}

	TVector3 TrajPosition = {TrajPosition_vec.X(), TrajPosition_vec.Y(),TrajPosition_vec.Z()};
	if(TrajPosition.Mag() == 0){continue;}

	TVector3 NextTrajPosition = {NextTrajPosition_vec.X(), NextTrajPosition_vec.Y(),NextTrajPosition_vec.Z()};
	TVector3 NextNextTrajPosition = {NextNextTrajPosition_vec.X(), NextNextTrajPosition_vec.Y(),NextNextTrajPosition_vec.Z()};
      
	if(TrajPosition.Angle(NextTrajPosition) < fAngleCut && TrajPosition.Angle(NextNextTrajPosition) <fAngleCut){std::cout << "broken" << std::endl;break;} 

	if(fAllowDynamicSliding){
	  StartPosition = {InitialTrack.LocationAtPoint(traj).X(), InitialTrack.LocationAtPoint(traj).Y(), InitialTrack.LocationAtPoint(traj).Z()};
	}
      }

      geo::Point_t  TrajPosition   = InitialTrack.LocationAtPoint(trajpoint);
      Direction_vec  = (TrajPosition - StartPosition).Unit();
    }
    else{ 
      //Loop over the trajectory points and find two corresponding trajectory points where the angle between themselves (or themsleves and the start position) is less the fMinAngle.
      for(unsigned int traj=0; traj<InitialTrack.NumberTrajectoryPoints()-2; ++traj){
	++trajpoint;

	geo::Vector_t TrajDirection_vec;

	if(fUseStartPos){
	  TrajDirection_vec         = InitialTrack.DirectionAtPoint(0);
	}
	else if(fAllowDynamicSliding && traj!=0){
	  TrajDirection_vec         = InitialTrack.DirectionAtPoint(traj-1);
	}
	else{
	  TrajDirection_vec         = InitialTrack.DirectionAtPoint(traj);
	}
	
	geo::Vector_t NextTrajDirection_vec     = InitialTrack.DirectionAtPoint(traj+1);
	geo::Vector_t NextNextTrajDirection_vec = InitialTrack.DirectionAtPoint(traj+2);

	TVector3 TrajDirection = {TrajDirection_vec.X(), TrajDirection_vec.Y(),TrajDirection_vec.Z()};
	if(TrajDirection.Mag() == 0){continue;}

	TVector3 NextTrajDirection = {NextTrajDirection_vec.X(), NextTrajDirection_vec.Y(),NextTrajDirection_vec.Z()};
	TVector3 NextNextTrajDirection = {NextNextTrajDirection_vec.X(), NextNextTrajDirection_vec.Y(),NextNextTrajDirection_vec.Z()};
      
	if(TrajDirection.Angle(NextTrajDirection) < fAngleCut && TrajDirection.Angle(NextNextTrajDirection) <fAngleCut){std::cout << "broken" << std::endl;break;} 
      }
      Direction_vec  = InitialTrack.DirectionAtPoint(trajpoint).Unit();
    }

    if(trajpoint == (int) InitialTrack.NumberTrajectoryPoints() -3){
      mf::LogError("ShowerDirection") << "Trajectory point not set."<< std::endl;
      return 1;
    }

    TVector3 Direction = {Direction_vec.X(), Direction_vec.Y(),Direction_vec.Z()};
    TVector3 DirectionErr = {-999,-999,-999};
    ShowerEleHolder.SetElement(Direction,DirectionErr,"ShowerDirection");

    return 0;
  }
}

  
DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerSmartTrackTrajectoryPointDirection)
  
