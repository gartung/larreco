//############################################################################
//### Name:        ShowerPCADirection                                      ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the shower direction using PCA         ###
//###              methods. Derviced from PandoraShowers Method            ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
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

//C++ Includes 
#include <iostream>

//Root Includes 
#include "TVector3.h"
#include "TMath.h"
#include "TPrincipal.h"
#include "TVector.h"

namespace ShowerRecoTools {

  
  class ShowerPCADirection:IShowerTool {
    
  public:
    
    ShowerPCADirection(const fhicl::ParameterSet& pset);
    
    ~ShowerPCADirection(); 
    
    //Generic Direction Finder
    int findMetric(const art::Ptr<recob::PFParticle>& pfparticle,
		   art::Event& Event,
		    reco::shower::ShowerPropertyHolder& ShowerPropHolder
		    ) override;

  private:
    
    // Define standard art tool interface
    void configure(const fhicl::ParameterSet& pset) override;

    TVector3 ShowerPCAVector(std::vector<art::Ptr<recob::SpacePoint> >& spacePoints_pfp, art::FindManyP<recob::Hit>& fmh, TVector3& ShowerCentre);
    double  RMSShowerGradient(std::vector<art::Ptr<recob::SpacePoint> >& sps, TVector3& ShowerCentre, TVector3& Direction);

    //Algorithm functions
    TVector3 ShowerCentre(std::vector<art::Ptr<recob::SpacePoint> >& showerspcs, art::FindManyP<recob::Hit>& fmh, float& TotalCharge);
    TVector3 SpacePointPosition(const art::Ptr<recob::SpacePoint>& sp);
    double SpacePointCharge(art::Ptr<recob::SpacePoint> sp,
			    art::FindManyP<recob::Hit>& fmh
			    );

    art::InputTag fPFParticleModuleLabel;
    int fNSegments; 
    bool fUseCollectionOnly;

  };
  
  
  ShowerPCADirection::ShowerPCADirection(const fhicl::ParameterSet& pset)
  {
    configure(pset);
  }
  
  ShowerPCADirection::~ShowerPCADirection()
  {
  }
  
  void ShowerPCADirection::configure(const fhicl::ParameterSet& pset)
  {
    fPFParticleModuleLabel  = pset.get<art::InputTag>("PFParticleModuleLabel","");
    fNSegments              = pset.get<float>        ("NSegments"); 
    fUseCollectionOnly      = pset.get<bool>         ("UseCollectionOnly");
  }

  int ShowerPCADirection::findMetric(const art::Ptr<recob::PFParticle>& pfparticle,
				      art::Event& Event,
				      reco::shower::ShowerPropertyHolder& ShowerPropHolder){

    // Get the assocated pfParicle vertex PFParticles
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("ShowerStartPosition") << "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }
    art::FindManyP<recob::SpacePoint> fmspp(pfpHandle, Event, fPFParticleModuleLabel);
    
    if (!fmspp.isValid()){
      throw cet::exception("ShowerStartPosition") << "Trying to get the spacepoint and failed. Something is not configured correctly. Stopping ";
      return 1;
    }

    //Get the spacepoints handle and the hit assoication
    art::Handle<std::vector<recob::SpacePoint> > spHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, spHandle)){
      throw cet::exception("ShowerStartPosition") << "Coquld not configure the spacepoint handle. Something is configured incorrectly. Stopping"; 
      return 1;
    } 
    art::FindManyP<recob::Hit> fmh(spHandle, Event, fPFParticleModuleLabel);
    if(!fmh.isValid()){
      throw cet::exception("ShowerStartPosition") << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }
    

    //Spacepoints
    std::vector<art::Ptr<recob::SpacePoint> > spacePoints_pfp = fmspp.at(pfparticle.key());

    //Find the PCA vector
    TVector3 ShowerCentre;
    TVector3 Eigenvector = ShowerPCAVector(spacePoints_pfp,fmh,ShowerCentre);
    
    TVector3 Direction        = { Eigenvector[0], Eigenvector[1], Eigenvector[2] }; 

    //Check if we are pointing the correct direction or not, First try the start position
    if(ShowerPropHolder.CheckShowerStartPosition()){
	
      //Get the General direction as the vector between the start position and the centre
      TVector3 StartPositionVec = ShowerPropHolder.GetShowerStartPosition();
      TVector3 GenralDir        = ShowerCentre - StartPositionVec;
          
      //Dot product
      double DotProduct = Direction.Dot(GenralDir);
      
      //If the dotproduct is negative the Direction needs Flipping
      if(DotProduct < 0){
	Eigenvector[0] = - Eigenvector[0];
	Eigenvector[1] = - Eigenvector[1];
	Eigenvector[2] = - Eigenvector[2];
      }
      
      ShowerPropHolder.SetShowerDirection(Eigenvector);
      return 0; 
    }
    
    //Otherwise Check against the RMS of the shower. Method adapated from EMShower Thanks Mike.
    double RMSGradient = RMSShowerGradient(spacePoints_pfp,ShowerCentre,Direction);
    if(RMSGradient < 0){
      	Eigenvector[0] = - Eigenvector[0];
	Eigenvector[1] = - Eigenvector[1];
	Eigenvector[2] = - Eigenvector[2];
      }
    
    ShowerPropHolder.SetShowerDirection(Direction);
    return 0;

	}

  //Function to calculate the RMS at segements of the shower and calculate the gradient of this. If negative then the direction is pointing the opposite way to the correct one
  double ShowerPCADirection::RMSShowerGradient(std::vector<art::Ptr<recob::SpacePoint> >& sps, TVector3& ShowerCentre, TVector3& Direction){
    
    //Get the length of the shower.
    TVector3 firstpoint = SpacePointPosition(sps[0]);
    TVector3 lastpoint  = SpacePointPosition(sps[sps.size()-1]);
    double length = (firstpoint-lastpoint).Mag();
    double segmentsize = length/fNSegments;
    
    std::map<float, std::vector<float> > len_segment_map;
    
    //Split the the spacepoints into segments.
    for(auto const& sp: sps){
      
      //Get the position of the spacepoint
      TVector3 pos = SpacePointPosition(sp) - ShowerCentre;

      //Get the the projected length
      double len = pos.Dot(Direction);

      //Get the length to the projection
      TVector3 perp = pos - len*Direction;
      double len_perp = perp.Mag();
      
      len_segment_map[segmentsize/len].push_back(len_perp);
    }
    
    float sumx  = 0;
    float sumy  = 0;
    float sumx2 = 0;
    float sumy2 = 0;
    float sumxy = 0;

    //Get the rms of the segments and caclulate the gradient. 
    for(auto const& segment: len_segment_map){

      float RMS = TMath::RMS((segment.second).begin(),(segment.second).end());

      //Calculate the gradient using regression
      sumx  += segment.first;
      sumy  += RMS;
      sumx2 += segment.first * segment.first;
      sumy2 += RMS*RMS;
      sumxy += RMS * segment.first;
    }

    double RMSgradient = (sumxy - sumx*sumy)/(sumx2 - sumx*sumx);
    return RMSgradient;
  
  } 

  //Function to calculate the shower direction using a charge weight 3D PCA calculation.
  TVector3 ShowerPCADirection::ShowerPCAVector(std::vector<art::Ptr<recob::SpacePoint> >& sps, art::FindManyP<recob::Hit>& fmh, TVector3& ShowerCentre){

    //Initialise the the PCA.
    TPrincipal *pca = new TPrincipal(3,"");

    float TotalCharge = 0;

    //Get the Shower Centre
    ShowerCentre = ShowerPCADirection::ShowerCentre(sps, fmh, TotalCharge);

    //Normalise the spacepoints, charge weight and add to the PCA.
    for(auto& sp: sps){
      
      //Normalise the spacepoint position.
      TVector3 sp_position = SpacePointPosition(sp);
      sp_position = sp_position - ShowerCentre;

      //Get the charge.
      float Charge = SpacePointCharge(sp,fmh);

      //Charge Weight
      float wht = TMath::Sqrt(Charge/TotalCharge);

      double sp_coord[3];
      sp_coord[0] = sp_position.X()*wht;
      sp_coord[1] = sp_position.Y()*wht;
      sp_coord[2] = sp_position.Z()*wht;

      //Add to the PCA
      pca->AddRow(sp_coord);
    }

    //Evaluate the PCA 
    pca->MakePrincipals();
   
    //Get the Eigenvectors.
    const TMatrixD * Eigenvectors = pca->GetEigenVectors();

    //Get the primary vector 
    const TVectorD EigenvectorD = (*Eigenvectors)[0];
    TVector3 Eigenvector = { EigenvectorD[0], EigenvectorD[1], EigenvectorD[2] };

    return Eigenvector;
   }

  //####################################################################################
  //Algorithm Functions
  TVector3 ShowerPCADirection::ShowerCentre(const std::vector<art::Ptr<recob::SpacePoint> >& showersps, art::FindManyP<recob::Hit>& fmh, float& totalCharge) {
    
    TVector3 pos, chargePoint = TVector3(0,0,0);
    
    //Loop over the spacepoints and get the charge weighted center.
    for(auto const& sp: showersps){
      
      //Get the position of the spacepoint 
      pos = SpacePointPosition(sp);
      
      //Get the associated hits 
      std::vector<art::Ptr<recob::Hit> > hits = fmh.at(sp.key());
      
      //Average the charge unless sepcified.
      float charge  = 0;
      float charge2 = 0;
      for(auto const& hit: hits){
	
	if(fUseCollectionOnly){
	  if(hit->SignalType() == geo::kCollection){ 
	    charge = hit->Integral();
	    break;
	  }
	}
	
	//Check if any of the points are not withing 2 sigma.
	if(!fUseCollectionOnly){
	  charge += hit->Integral();
	  charge2 += hit->Integral();
	}
      }
      
      if(!fUseCollectionOnly){
	//Calculate the unbiased standard deviation and mean. 
	float mean = charge/((float) hits.size()); 
	float rms  = TMath::Sqrt((charge2 - charge*charge)/((float)(hits.size()-1)));
	
	charge = 0;
	for(auto const& hit: hits){
	  if(hit->Integral() > (mean - 2*rms) && hit->Integral() < (mean + 2*rms))  
	    charge += hit->Integral();
	}
      }
      
      chargePoint += charge * pos;
      totalCharge += charge;
      
      if(charge == 0){
	mf::LogWarning("ShowerStartPosition") << "Averaged charge, within 2 sigma, for a spacepoint is zero, Maybe this not a good method";
      }
    }
    
    double intotalcharge = 1/totalCharge;
    TVector3 centre = chargePoint *  intotalcharge;
    return centre;
    
  }
  
  TVector3 ShowerPCADirection::SpacePointPosition(const art::Ptr<recob::SpacePoint>& sp){
    
    const Double32_t* sp_xyz = sp->XYZ();
    TVector3 sp_postiion = {sp_xyz[0], sp_xyz[1], sp_xyz[2]};
    return sp_postiion;
  }

  
  double ShowerPCADirection::SpacePointCharge(art::Ptr<recob::SpacePoint> sp,
					      art::FindManyP<recob::Hit>& fmh
					      ){
    return 1;
  }
  


}
  
DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerPCADirection)
  
