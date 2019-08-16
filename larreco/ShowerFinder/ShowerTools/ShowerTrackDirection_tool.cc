//############################################################################
//### Name:        ShowerTrackDirection                                    ###
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

#include "TFile.h"

namespace ShowerRecoTools {


  class ShowerTrackDirection:IShowerTool {

    public:

      ShowerTrackDirection(const fhicl::ParameterSet& pset);

      ~ShowerTrackDirection();

      //Generic Direction Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      //Algoritma
      shower::TRACSAlg       fTRACSAlg;

      //fcl
      bool fUsePandoraVertex;
      bool fUsePositionInfo;
      bool fDebugEVD;

  };


  ShowerTrackDirection::ShowerTrackDirection(const fhicl::ParameterSet& pset)
    :fTRACSAlg(pset.get<fhicl::ParameterSet>("TRACSAlg"))
  {
    fUsePandoraVertex         = pset.get<bool>("UsePandoraVertex");
    fUsePositionInfo          = pset.get<bool>("UsePositionInfo");
    fDebugEVD                 = pset.get<bool>("DebugEVD");
  }

  ShowerTrackDirection::~ShowerTrackDirection()
  {
  }

  int ShowerTrackDirection::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){

    //Check the Track has been defined
    if(!ShowerEleHolder.CheckElement("InitialTrack")){
      mf::LogError("ShowerTrackDirection") << "Initial track not set"<< std::endl;
      return 0;
    }

    //Check the start position is set.
    if(fUsePandoraVertex && !ShowerEleHolder.CheckElement("ShowerStartPosition")){
      mf::LogError("ShowerTrackDirection") << "Start position not set, returning "<< std::endl;
      return 0;
    }

    //Get the track
    recob::Track InitialTrack;
    ShowerEleHolder.GetElement("InitialTrack",InitialTrack);

    if(fUsePositionInfo){
      geo::Point_t StartPosition;
      if(fUsePandoraVertex){
        TVector3 StartPosition_vec = {-999,-999,-999};
        ShowerEleHolder.GetElement("ShowerStartPosition",StartPosition_vec);
        StartPosition.SetCoordinates(StartPosition_vec.X(),StartPosition_vec.Y(),StartPosition_vec.Z());
      }
      else{
        StartPosition = InitialTrack.Start();
      }

      //Calculate the mean direction and the the standard deviation
      float sumX=0, sumX2=0;
      float sumY=0, sumY2=0;
      float sumZ=0, sumZ2=0;
      for(unsigned int traj=0; traj< InitialTrack.NumberTrajectoryPoints(); ++traj){

        auto flags = InitialTrack.FlagsAtPoint(traj);
        if(flags.isSet(recob::TrajectoryPointFlags::InvalidHitIndex)
            || flags.isSet(recob::TrajectoryPointFlagTraits::NoPoint))
        {continue;}


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

      TVector3 Direction_Mean = {0,0,0};
      int N = 0;
      //Remove trajectory points from the mean that are not with one sigma.
      for(unsigned int traj=0; traj< InitialTrack.NumberTrajectoryPoints(); ++traj){

        auto flags = InitialTrack.FlagsAtPoint(traj);
        if(flags.isSet(recob::TrajectoryPointFlags::InvalidHitIndex)
            || flags.isSet(recob::TrajectoryPointFlagTraits::NoPoint))
        {continue;}


        geo::Point_t TrajPosition = InitialTrack.LocationAtPoint(traj);
        geo::Vector_t Direction   = (TrajPosition - StartPosition).Unit();

        if((TMath::Abs((Direction-Mean).X()) < 1*RMSX) &&
            (TMath::Abs((Direction-Mean).Y()) < 1*RMSY) &&
            (TMath::Abs((Direction-Mean).Z()) < 1*RMSZ)){
          TVector3 Direction_vec = {Direction.X(),Direction.Y(),Direction.Z()};
          if(Direction_vec.Mag() == 0){continue;}
          Direction_Mean += Direction_vec;
          ++N;
        }
      }


      //Take the mean value
      if(N > 0){
        TVector3 Direction    = Direction_Mean.Unit();
        TVector3 DirectionErr = {RMSX,RMSY,RMSZ};
        ShowerEleHolder.SetElement(Direction,DirectionErr,"ShowerDirection");
      }
      else{
        mf::LogError("ShowerTrackDirection") << "None of the points are within 1 sigma"<< std::endl;
        return 1;
      }

      if (fDebugEVD){
        fTRACSAlg.DebugEVD(pfparticle,Event,ShowerEleHolder);
      }
      return 0;

    }else{ // if(fUsePositionInfo)
      float sumX=0, sumX2=0;
      float sumY=0, sumY2=0;
      float sumZ=0, sumZ2=0;
      for(unsigned int traj=0; traj< InitialTrack.NumberTrajectoryPoints(); ++traj){

        auto flags = InitialTrack.FlagsAtPoint(traj);
        if(flags.isSet(recob::TrajectoryPointFlags::InvalidHitIndex)
            || flags.isSet(recob::TrajectoryPointFlagTraits::NoPoint))
        {continue;}


        geo::Vector_t  Direction = InitialTrack.DirectionAtPoint(traj);
        sumX += Direction.X(); sumX2 += Direction.X()*Direction.X();
        sumY += Direction.Y(); sumY2 += Direction.Y()*Direction.Y();
        sumZ += Direction.Z(); sumZ2 += Direction.Z()*Direction.Z();
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

      //Remove trajectory points from the mean that are not with one sigma.
      float N = 0.;
      TVector3 Direction_Mean = {0,0,0};
      for(unsigned int traj=0; traj<InitialTrack.NumberTrajectoryPoints(); ++traj){

        auto flags = InitialTrack.FlagsAtPoint(traj);
        if(flags.isSet(recob::TrajectoryPointFlags::InvalidHitIndex)
            || flags.isSet(recob::TrajectoryPointFlagTraits::NoPoint))
        {continue;}


        geo::Vector_t Direction = InitialTrack.DirectionAtPoint(traj).Unit();
        if((TMath::Abs((Direction-Mean).X()) < 1*RMSX) &&
            (TMath::Abs((Direction-Mean).Y()) < 1*RMSY) &&
            (TMath::Abs((Direction-Mean).Z()) < 1*RMSZ)){
          TVector3 Direction_vec = {Direction.X(),Direction.Y(),Direction.Z()};
          if(Direction_vec.Mag() == 0){continue;}
          Direction_Mean += Direction_vec;
          ++N;
        }
      }

      //Take the mean value
      if(N>0){
        TVector3 Direction = Direction_Mean.Unit();
        TVector3 DirectionErr = {RMSX,RMSY,RMSZ};
        ShowerEleHolder.SetElement(Direction,DirectionErr,"ShowerDirection");
      }
      else{
        mf::LogError("ShowerTrackDirection") << "None of the points are within 1 sigma"<< std::endl;
        return 1;
      }

      if (fDebugEVD){
        fTRACSAlg.DebugEVD(pfparticle,Event,ShowerEleHolder);
      }


    }
    return 0;
  }
}


DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackDirection)



