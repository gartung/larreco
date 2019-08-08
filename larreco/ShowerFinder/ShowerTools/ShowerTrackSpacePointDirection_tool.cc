//############################################################################
//### Name:        ShowerTrackSpacePointDirection                          ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the shower direction using the         ###
//###              the average direction of theinitial track spacepoints.  ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"

//LArSoft Includes 
#include "larcore/Geometry/Geometry.h"
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

  
  class ShowerTrackSpacePointDirection:IShowerTool {
    
  public:
    
    ShowerTrackSpacePointDirection(const fhicl::ParameterSet& pset);
    
    ~ShowerTrackSpacePointDirection(); 
    
    //Generic Direction Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
			 art::Event& Event,
			 reco::shower::ShowerElementHolder& ShowerEleHolder
			 ) override;
    
  private:
    
    // Define standard art tool interface
    void configure(const fhicl::ParameterSet& pset) override;

    //fcl
    bool fUsePandoraVertex;
    
    //Algorithm function
    shower::SBNShowerAlg fSBNShowerAlg;
  };
  
  
  ShowerTrackSpacePointDirection::ShowerTrackSpacePointDirection(const fhicl::ParameterSet& pset)
    : fSBNShowerAlg(pset.get<fhicl::ParameterSet>("SBNShowerAlg"))
  {
    configure(pset);
  }

  ShowerTrackSpacePointDirection::~ShowerTrackSpacePointDirection()
  {
  }
  
  void ShowerTrackSpacePointDirection::configure(const fhicl::ParameterSet& pset)
  {
    fUsePandoraVertex       = pset.get<bool>         ("UsePandoraVertex");
  }

  int ShowerTrackSpacePointDirection::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
						art::Event& Event,
						reco::shower::ShowerElementHolder& ShowerEleHolder){

    std::cout << "on shower track hit direction" << std::endl;

    //Check the Track Hits has been defined
    if(!ShowerEleHolder.CheckElement("InitialTrackSpacePoints")){
      mf::LogError("ShowerTrackSpacePointDirection") << "Initial track spacepoints not set"<< std::endl;
      return 0;
    }
    
    //Check the start position is set.
    if(fUsePandoraVertex && !ShowerEleHolder.CheckElement("ShowerStartPosition")){
      mf::LogError("ShowerTrackSpacePointDirection") << "Start position not set, returning "<< std::endl;
      return 0;
    }

    //Get the start poistion 
    TVector3 StartPosition = {-999,-999,-999};
    if(fUsePandoraVertex){
      ShowerEleHolder.GetElement("ShowerStartPosition",StartPosition);
    }
    else{ 
      //Check the Tracks has been defined
      if(!ShowerEleHolder.CheckElement("InitialTrack")){
	mf::LogError("ShowerTrackSpacePointDirection") << "Initial track not set"<< std::endl;
	return 0;
      }
      recob::Track InitialTrack;
      ShowerEleHolder.GetElement("InitialTrack",InitialTrack);
      geo::Point_t Start_point = InitialTrack.Start();
      StartPosition = {Start_point.X(),Start_point.Y(),Start_point.Z()};
    }

    //Get the initial track hits.
    std::vector<art::Ptr<recob::SpacePoint> > intitaltrack_sp;
    ShowerEleHolder.GetElement("InitialTrackSpacePoints",intitaltrack_sp);
    
    //Calculate the mean direction and the the standard deviation
    float sumX=0, sumX2=0;
    float sumY=0, sumY2=0;
    float sumZ=0, sumZ2=0;

    //Get the spacepoints associated to the track hit 
    for(auto const& sp: intitaltrack_sp){
      
      //Get the direction relative to the start positon 
      TVector3 pos = fSBNShowerAlg.SpacePointPosition(sp) - StartPosition;
      if(pos.Mag() == 0){continue;}
      
      sumX = pos.X(); sumX2 += pos.X()*pos.X();
      sumY = pos.Y(); sumY2 += pos.Y()*pos.Y();
      sumZ = pos.Z(); sumZ2 += pos.Z()*pos.Z();
    }

    float NumSps = intitaltrack_sp.size();
    TVector3 Mean = {sumX/NumSps,sumY/NumSps,sumZ/NumSps};
    Mean = Mean.Unit();
      
    float RMSX = 999;
    float RMSY = 999;
    float RMSZ = 999;
    if(sumX2/NumSps - ((sumX/NumSps)*((sumX/NumSps))) > 0){
      RMSX = TMath::Sqrt(sumX2/NumSps - ((sumX/NumSps)*((sumX/NumSps))));
    }
    if(sumY2/NumSps - ((sumY/NumSps)*((sumY/NumSps))) > 0){
      RMSY = TMath::Sqrt(sumY2/NumSps - ((sumY/NumSps)*((sumY/NumSps))));
    }
    if(sumZ2/NumSps - ((sumZ/NumSps)*((sumZ/NumSps))) > 0){
      RMSZ = TMath::Sqrt(sumZ2/NumSps - ((sumZ/NumSps)*((sumZ/NumSps))));
    }

    
    //Loop over the spacepoints and remove ones the relative direction is not within one sigma.
    TVector3 Direction_Mean = {0,0,0};
    int N = 0; 
    for(auto const sp: intitaltrack_sp){
      TVector3 Direction = fSBNShowerAlg.SpacePointPosition(sp) - StartPosition;
      if((TMath::Abs((Direction-Mean).X()) < 1*RMSX) && 
	 (TMath::Abs((Direction-Mean).Y())< 1*RMSY) && 
	 (TMath::Abs((Direction-Mean).Z()) < 1*RMSZ)){
	if(Direction.Mag() == 0){continue;}
	++N;
	Direction_Mean += Direction;
      }
    }

    if(N>0){
      //Take the mean value
      TVector3 Direction = Direction_Mean.Unit();
      ShowerEleHolder.SetElement(Direction,"ShowerDirection");
    }  
    else{
      mf::LogError("ShowerTrackSpacePointDirection") << "None of the points are within 1 sigma"<< std::endl;
      return 1;
    }
    return 0;
  }
}

  
DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackSpacePointDirection)
  
