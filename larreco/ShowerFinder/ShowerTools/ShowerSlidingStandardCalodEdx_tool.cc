//############################################################################
//### Name:        ShowerSlidingStandardCalodEdx                           ###
//### Author:      Dominic Barker (dominic.barker@sheffield.ac.uk)         ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the dEdx of the start track of the     ###
//###              shower using the standard calomitry module. This        ###
//###              takes the sliding fit trajectory to make a 3D dEdx.     ###
//###              This module is best used with the sliding linear fit    ###
//###              and ShowerTrackTrajToSpacepoint                         ###
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
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardataobj/RecoBase/Hit.h" 
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larreco/RecoAlg/SBNShowerAlg.h"

//C++ Includes 
#include <iostream>
#include <vector> 

//Root Includes 
#include "TVector3.h"

namespace ShowerRecoTools{

  class ShowerSlidingStandardCalodEdx:IShowerTool {

  public:

    ShowerSlidingStandardCalodEdx(const fhicl::ParameterSet& pset);

    ~ShowerSlidingStandardCalodEdx(); 

    //Generic Direction Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
			 art::Event& Event,
			 reco::shower::ShowerElementHolder& ShowerEleHolder
			 ) override;

  private:

    // Define standard art tool interface
    void configure(const fhicl::ParameterSet& pset) override;

    art::ServiceHandle<geo::Geometry> fGeom;
    calo::CalorimetryAlg fCalorimetryAlg;
    shower::SBNShowerAlg fSBNShowerAlg;
    detinfo::DetectorProperties const* fDetProp; 

    float fMinAngleToWire;
    float fShapingTime;
    float fMinDistCutOff;
    float fMaxDist;
    float fdEdxTrackLength;
    bool  fUseMedian;
    art::InputTag fPFParticleModuleLabel;

    bool fCutStartPosition;
    bool fTrajDirection;
  };


  ShowerSlidingStandardCalodEdx::ShowerSlidingStandardCalodEdx(const fhicl::ParameterSet& pset):
    fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
    fSBNShowerAlg(pset.get<fhicl::ParameterSet>("SBNShowerAlg")),
    fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())
  {
    configure(pset);
  }

  ShowerSlidingStandardCalodEdx::~ShowerSlidingStandardCalodEdx()
  {
  }

  void ShowerSlidingStandardCalodEdx::configure(const fhicl::ParameterSet& pset)
  {
    fMinDistCutOff         = pset.get<float>("MinDistCutOff");
    fMaxDist               = pset.get<float>("MaxDist");
    fMinAngleToWire        = pset.get<float>("MinAngleToWire");
    fShapingTime           = pset.get<float>("ShapingTime");
    fdEdxTrackLength       = pset.get<float>("dEdxTrackLength");
    fUseMedian             = pset.get<bool> ("UseMedian");
    fPFParticleModuleLabel = pset.get<art::InputTag>("PFParticleModuleLabel");
    
    fCutStartPosition = pset.get<bool> ("CutStartPosition");
    fTrajDirection    = pset.get<bool> ("TrajDirection");
  }

  int ShowerSlidingStandardCalodEdx::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
					       art::Event& Event,
					       reco::shower::ShowerElementHolder& ShowerEleHolder
					       ){


    // Shower dEdx calculation
    std::cout << "hello world dEdx" << std::endl;

    if(!ShowerEleHolder.CheckElement("ShowerStartPosition")){
      mf::LogError("ShowerSlidingStandardCalodEdx") << "Start position not set, returning "<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement("InitialTrackSpacePoints")){
      mf::LogError("ShowerSlidingStandardCalodEdx") << "Initial Track Spacepoints is not set returning"<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement("InitialTrack")){
      mf::LogError("ShowerSlidingStandardCalodEdx") << "Initial Track is not set"<< std::endl;
      return 1;
    }

    //Get the initial track hits
    std::vector<art::Ptr<recob::SpacePoint> > tracksps;
    ShowerEleHolder.GetElement("InitialTrackSpacePoints",tracksps);

    if(tracksps.size() == 0){
      mf::LogWarning("ShowerSlidingStandardCalodEdx") << "no spacepointsin the initial track" << std::endl;
      return 0;
    }

    // Get the spacepoints
    art::Handle<std::vector<recob::SpacePoint> > spHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, spHandle)){
      throw cet::exception("ShowerTrackTrajToSpacepoint") << "Could not configure the spacepoint handle. Something is configured incorrectly. Stopping";
      return 1;
    }
    
    // Get the hits associated with the space points
    art::FindManyP<recob::Hit> fmsp(spHandle, Event, fPFParticleModuleLabel);
    if(!fmsp.isValid()){
      throw cet::exception("ShowerTrackTrajToSpacepoint") << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }



    //Only consider hits in the same tpcs as the vertex.
    TVector3 ShowerStartPosition = {-999,-999,-999};
    ShowerEleHolder.GetElement("ShowerStartPosition",ShowerStartPosition);
    geo::TPCID vtxTPC = fGeom->FindTPCAtPosition(ShowerStartPosition); 

    //Get the initial track
    recob::Track InitialTrack;
    ShowerEleHolder.GetElement("InitialTrack",InitialTrack);
    
    //Don't care that I could use a vector. 
    std::map<int,std::vector<double > > dEdx_vec; 
    std::map<int,std::vector<double> >  dEdx_vecErr;
    std::map<int,int> num_hits;

    for(geo::PlaneID plane_id: fGeom->IteratePlaneIDs()){
      dEdx_vec[plane_id.Plane] = {};
      dEdx_vecErr[plane_id.Plane] = {}; 
      num_hits[plane_id.Plane] = 0;
    }

    //Loop over the spacepoints
    for(auto const sp: tracksps){

      //Get the associated hit
      std::vector<art::Ptr<recob::Hit> > hits = fmsp.at(sp.key()); 
      if(hits.size() == 0){
	mf::LogWarning("ShowerSlidingStandardCalodEdx") << "no hit for the spacepoint. This suggest the find many is wrong."<< std::endl;
	continue;
      }
      const art::Ptr<recob::Hit> hit = hits[0];
      double wirepitch = fGeom->WirePitch((geo::PlaneID)hit->WireID()); 

      //Only consider hits in the same tpc
      geo::PlaneID planeid = hit->WireID();
      geo::TPCID TPC = planeid.asTPCID();
      if (TPC !=vtxTPC){continue;}

      //Ignore spacepoints within a few wires of the vertex.
      double dist_from_start = (fSBNShowerAlg.SpacePointPosition(sp) - ShowerStartPosition).Mag();

      std::cout << "hit on plane: " << planeid.Plane << " dist_from_start: " << dist_from_start << std::endl;

      if(fCutStartPosition){
	if(dist_from_start < fMinDistCutOff*wirepitch){std::cout << "too close" << std::endl; continue;}
	
	if(dist_from_start > fdEdxTrackLength){std::cout << " too far " << std::endl;continue;}
      }

      //Find the closest trajectory point of the track. These should be in order if the user has used ShowerTrackTrajToSpacepoint_tool but the sake of gernicness I'll get the cloest sp.
      unsigned int index = 999;
      double MinDist = 999;  
      for(unsigned int traj=0; traj< InitialTrack.NumberTrajectoryPoints(); ++traj){

	geo::Point_t TrajPositionPoint = InitialTrack.LocationAtPoint(traj);
	TVector3 TrajPosition = {TrajPositionPoint.X(),TrajPositionPoint.Y(),TrajPositionPoint.Z()};


	//ignore bogus info.
	auto flags = InitialTrack.FlagsAtPoint(traj);
	if(flags.isSet(recob::TrajectoryPointFlagTraits::NoPoint))
	  {continue;}

	TVector3 pos = fSBNShowerAlg.SpacePointPosition(sp) - TrajPosition;
	
	if(pos.Mag() < MinDist && pos.Mag()< fMaxDist*wirepitch){
	  MinDist = pos.Mag();
	  index = traj;
	}
      }

      //If there is no matching trajectory point then bail. 
      if(index == 999){std::cout << "traj not found" << std::endl; continue;}

      std::cout << "index: " << index << std::endl;

      geo::Point_t TrajPositionPoint = InitialTrack.LocationAtPoint(index);
      TVector3 TrajPosition = {TrajPositionPoint.X(),TrajPositionPoint.Y(),TrajPositionPoint.Z()};
      
      geo::Point_t TrajPositionStartPoint = InitialTrack.LocationAtPoint(0);
      TVector3 TrajPositionStart = {TrajPositionStartPoint.X(),TrajPositionStartPoint.Y(),TrajPositionStartPoint.Z()};
      
      //Ignore values with 0 mag from the start position 
      if((TrajPosition - TrajPositionStart).Mag() == 0){std::cout << "mag is 0" << std::endl;continue;}
      if((TrajPosition - ShowerStartPosition).Mag() == 0){std::cout << "mag1 is 0" << std::endl;continue;}

   
      std::cout << "(TrajPosition-TrajPositionStart).Mag(): " << (TrajPosition-TrajPositionStart).Mag() << std::endl;
      
      if((TrajPosition-TrajPositionStart).Mag() < fMinDistCutOff*wirepitch){std::cout << "minmidistance cut off" << std::endl; continue;}    
      
      if((TrajPosition-TrajPositionStart).Mag() > fdEdxTrackLength){std::cout << " too far " << std::endl;continue;}




      //Get the direction of the trajectory point 
      geo::Vector_t TrajDirection_vec = InitialTrack.DirectionAtPoint(index); 
      TVector3 TrajDirection = {TrajDirection_vec.X(),TrajDirection_vec.Y(),TrajDirection_vec.Z()};

      //If the direction is in the same direction as the wires within some tolerance the hit finding struggles. Let remove these. 
      TVector3 PlaneDirection = fGeom->Plane(planeid).GetIncreasingWireDirection(); 

      if(TrajDirection.Angle(PlaneDirection) < fMinAngleToWire){ std::cout << "remove from angle cut" << std::endl;continue;}

      //If the direction is too much into the wire plane then the shaping amplifer cuts the charge. Lets remove these events.
      double velocity = fDetProp->DriftVelocity(fDetProp->Efield(), fDetProp->Temperature());
      double distance_in_x = TrajDirection.X()*(wirepitch/TrajDirection.Dot(PlaneDirection));
      double time_taken = TMath::Abs(distance_in_x/velocity); 

      //Shaping time doesn't seem to exist in a global place so add it as a fcl.
      if(fShapingTime < time_taken){std::cout << "move for shaping time" << std::endl; continue;}

      if(fTrajDirection){
	std::cout << "using track direction" << std::endl;
	ShowerEleHolder.GetElement("ShowerDirection",TrajDirection);
      }

	//If we still exist then we can be used in the calculation. Calculate the 3D pitch
      double trackpitch = (TrajDirection*(wirepitch/TrajDirection.Dot(PlaneDirection))).Mag();

      double angleToVert = fGeom->WireAngleToVertical(fGeom->Plane(planeid.Plane).View(),planeid) - 0.5*TMath::Pi();
      double cosgamma = std::abs(sin(angleToVert)*TrajDirection.Y()+cos(angleToVert)*TrajDirection.Z());

      std::cout << "trackpitch: " << trackpitch << " angle: " << TMath::ATan(TrajDirection.X()/TrajDirection.Z()) << " x pitch: " << wirepitch*TMath::Tan(TMath::ATan(TrajDirection.X()/TrajDirection.Z())) << " wirepitch: " << wirepitch << " plane: " << planeid.Plane << " pitch w: " << wirepitch/cosgamma << " pitch again: " << TMath::Sqrt(((wirepitch/cosgamma)*(wirepitch/cosgamma)) + (wirepitch*TMath::Tan(TMath::ATan(TrajDirection.X()/TrajDirection.Z()))*wirepitch*TMath::Tan(TMath::ATan(TrajDirection.X()/TrajDirection.Z()))))  << std::endl;

      if(trackpitch == 0){
	mf::LogWarning("ShowerSlidingStandardCalodEdx") << "pitch is zero so we are not using the hit "<< std::endl;
	continue;
      }

      //Calculate the dQdx
      double dQdx = hit->Integral()/trackpitch;
      
      //Calculate the dEdx
      double dEdx = fCalorimetryAlg.dEdx_AREA(dQdx, hit->PeakTime(), planeid.Plane);
      std::cout << "TrajDirection.X(): " << TrajDirection.X() << " Y: " << TrajDirection.Y() << TrajDirection.Z() << std::endl;
      std::cout << "dQdx: " <<dQdx << " dEdx: " << dEdx << " trackpitch: " << trackpitch << std::endl;
      

      //Add the value to the dEdx 
      dEdx_vec[planeid.Plane].push_back(dEdx);

      //Iterate the number of hits on the plane
      ++num_hits[planeid.Plane];
    }


    //Never have the stats to do a landau fit and get the most probable value. User decides if they want the median value or the mean.
    std::vector<double> dEdx_val;
    std::vector<double> dEdx_valErr; 
    for(auto const& dEdx_plane: dEdx_vec){

      if((dEdx_plane.second).size() == 0){
	dEdx_val.push_back(-999);
	dEdx_valErr.push_back(-999);
	continue;
      }

      if(fUseMedian){dEdx_val.push_back(TMath::Median((dEdx_plane.second).size(), &(dEdx_plane.second)[0]));}
      else{
	//Else calculate the mean value.
	double dEdx_mean = 0;
	for(auto const& dEdx: dEdx_plane.second){
	  if(dEdx > 10 || dEdx < 0){continue;}
	  dEdx_mean += dEdx;
	}
	dEdx_val.push_back(dEdx_mean/(float)(dEdx_plane.second).size());
      }
    }

    //Work out which is the best plane from the most hits.
    int max_hits   = -999;
    int best_plane = -999;
    for(auto const& num_hits_plane: num_hits){
      if(num_hits_plane.second > max_hits){
	best_plane = num_hits_plane.first;
	max_hits = num_hits_plane.second;
      }
    }

    //To do
    ShowerEleHolder.SetElement(dEdx_val,dEdx_valErr,"ShowerdEdx");
    ShowerEleHolder.SetElement(best_plane,"ShowerBestPlane");

    std::cout<<"Best Plane: "<<best_plane<<" and plane with most hits: "<< max_hits<<std::endl;
    std::cout<<"main dEdx: "<<dEdx_val[best_plane]<<std::endl;

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerSlidingStandardCalodEdx)

