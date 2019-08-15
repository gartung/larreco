//############################################################################
//### Name:        ShowerTrackTrajectoryPointDirection                     ###
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


  class ShowerTrackTrajectoryPointDirection:IShowerTool {

    public:

      ShowerTrackTrajectoryPointDirection(const fhicl::ParameterSet& pset);

      ~ShowerTrackTrajectoryPointDirection();

      //Generic Direction Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      //fcl
      bool fUsePandoraVertex;
      bool fDebugEVD;
      bool fUsePositonInfo;
      int  fTrajPoint;
  };


  ShowerTrackTrajectoryPointDirection::ShowerTrackTrajectoryPointDirection(const fhicl::ParameterSet& pset)
  {
    fUsePandoraVertex         = pset.get<bool>("UsePandoraVertex");
    fDebugEVD                 = pset.get<bool>("DebugEVD");
    fUsePositonInfo           = pset.get<bool>("UsePositonInfo");
    fTrajPoint                = pset.get<int>("TrajPoint");
  }

  ShowerTrackTrajectoryPointDirection::~ShowerTrackTrajectoryPointDirection()
  {
  }


  int ShowerTrackTrajectoryPointDirection::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle, art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){

    //Check the Track has been defined
    if(!ShowerEleHolder.CheckElement("InitialTrack")){
      mf::LogError("ShowerTrackTrajectoryPointDirection") << "Initial track not set"<< std::endl;
      return 1;
    }
    recob::Track InitialTrack;
    ShowerEleHolder.GetElement("InitialTrack",InitialTrack);

    if((int)InitialTrack.NumberTrajectoryPoints()-1 < fTrajPoint){
      mf::LogError("ShowerTrackTrajectoryPointDirection") << "Less that fTrajPoint trajectory points, bailing."<< std::endl;
      fTrajPoint = InitialTrack.NumberTrajectoryPoints()-1;
    }

    //ignore bogus info.
    auto flags = InitialTrack.FlagsAtPoint(fTrajPoint);
    if(flags.isSet(recob::TrajectoryPointFlags::InvalidHitIndex)
        || flags.isSet(recob::TrajectoryPointFlagTraits::NoPoint))
    {
      mf::LogError("ShowerTrackTrajectoryPointDirection") << "Bogus trajectory poitn bailing."<< std::endl;
      return 1;
    }


    geo::Vector_t Direction_vec;
    //Get the difference between the point and the start position.
    if(fUsePositonInfo){
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
      //Get the specific trajectory point and look and and the direction from the start position
      geo::Point_t  TrajPosition = InitialTrack.LocationAtPoint(fTrajPoint);
      Direction_vec  = (TrajPosition - StartPosition).Unit();
    }
    else{
      //Use the direction of the trajection at tat point;
      Direction_vec = InitialTrack.DirectionAtPoint(fTrajPoint);
    }

    TVector3 Direction = {Direction_vec.X(), Direction_vec.Y(),Direction_vec.Z()};
    TVector3 DirectionErr = {-999,-999,-999};
    ShowerEleHolder.SetElement(Direction,DirectionErr,"ShowerDirection");
    return 0;
  }
}


DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackTrajectoryPointDirection)

