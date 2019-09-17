//############################################################################
//### Name:        ShowerResidualTrackHitFinder                            ###
//### Author:      You                                                     ###
//### Date:        13.05.19                                                ###
//### Description: Generic form of the shower tools                        ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"

//Root Includes
#include "TVector3.h"
#include "TMath.h"
#include "TPrincipal.h"


namespace ShowerRecoTools {


  class ShowerResidualTrackHitFinder: public IShowerTool {
    
  public:
    
    ShowerResidualTrackHitFinder(const fhicl::ParameterSet& pset);
    
    ~ShowerResidualTrackHitFinder();
    
    //Generic Direction Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
			   art::Event& Event,
			   reco::shower::ShowerElementHolder& ShowerEleHolder
			   ) override;
    
  private:
    
    double CalculateResidual(std::vector<art::Ptr<recob::SpacePoint> >& sps, 
			     TVector3& PCAEigenvector, 
			     TVector3& ShowerStartPosition);

    TVector3 ShowerPCAVector(std::vector<art::Ptr<recob::SpacePoint> >& sps, 
			     art::FindManyP<recob::Hit>& fmh,
			     TVector3& ShowerStartPosition);


    //Services
    detinfo::DetectorProperties const* fDetProp;
    
    art::InputTag fPFParticleModuleLabel;
    bool          fUseShowerDirection;
    bool          fChargeWeighted;
    bool          fForwardHitsOnly;
    float         fMaxResidualDiff;
    
  };


  ShowerResidualTrackHitFinder::ShowerResidualTrackHitFinder(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>()),
    fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel","")),
    fUseShowerDirection(pset.get<bool>("UseShowerDirection")),
    fChargeWeighted(pset.get<bool>("ChargeWeighted")),
    fForwardHitsOnly(pset.get<bool>("ForwardHitsOnly")),
    fMaxResidualDiff(pset.get<float>("MaxResidualDiff"))
  {
  }

  ShowerResidualTrackHitFinder::~ShowerResidualTrackHitFinder()
  {
  }


  int ShowerResidualTrackHitFinder::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){

    //This is all based on the shower vertex being known. If it is not lets not do the track
    if(!ShowerEleHolder.CheckElement("ShowerStartPosition")){
      mf::LogError("ShowerResidualTrackHitFinder") << "Start position not set, returning "<< std::endl;
      return 1;
    }
    
    // Get the assocated pfParicle Handle
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("ShowerResidualTrackHitFinder") << "Could not get the pandora pf particles. Something is not cofingured correctly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    // Get the spacepoint - PFParticle assn
    art::FindManyP<recob::SpacePoint> fmspp(pfpHandle, Event, fPFParticleModuleLabel);
    if (!fmspp.isValid()){
      throw cet::exception("ShowerResidualTrackHitFinder") << "Trying to get the spacepoint and failed. Something is not configured correctly. Stopping ";
      return 1;
    }

    // Get the spacepoints
    art::Handle<std::vector<recob::SpacePoint> > spHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, spHandle)){
      throw cet::exception("ShowerResidualTrackHitFinder") << "Could not configure the spacepoint handle. Something is configured incorrectly. Stopping";
      return 1;
    }

    // Get the hits associated with the space points
    art::FindManyP<recob::Hit> fmh(spHandle, Event, fPFParticleModuleLabel);
    if(!fmh.isValid()){
      throw cet::exception("ShowerResidualTrackHitFinder") << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }

    //Holder for the track stub spacepoints.
    std::vector<art::Ptr<recob::SpacePoint> > track_sps;

    // Get the SpacePoints
    std::vector<art::Ptr<recob::SpacePoint> > spacePoints = fmspp.at(pfparticle.key());

    //We cannot progress with no spacepoints.
    if(spacePoints.size() == 0){
      mf::LogError("ShowerResidualTrackHitFinder") << "No space points, returning "<< std::endl;
      return 1;
    }

    TVector3 ShowerStartPosition = {-999,-999,-999};
    ShowerEleHolder.GetElement("ShowerStartPosition",ShowerStartPosition);


    //Decide if the you want to use the direction of the shower or make one.
    if(fUseShowerDirection){ 

      if(!ShowerEleHolder.CheckElement("ShowerDirection")){
      mf::LogError("ShowerResidualTrackHitFinder") << "Direction not set, returning "<< std::endl;
      return 1;
      }

      TVector3 ShowerDirection     = {-999,-999,-999};
      ShowerEleHolder.GetElement("ShowerDirection",ShowerDirection);

      //Order the spacepoints
      IShowerTool::GetTRACSAlg().OrderShowerSpacePoints(spacePoints,ShowerStartPosition,ShowerDirection);
      
      //Remove the back hits if requird.
      if (fForwardHitsOnly){
	int back_sps=0;
	for (auto spacePoint : spacePoints){
	  double proj = IShowerTool::GetTRACSAlg().SpacePointProjection(spacePoint,ShowerStartPosition, ShowerDirection);
	  if(proj<0){
	    ++back_sps;
	  }
	  if(proj>0){
	    break;
	  }
	}
	spacePoints.erase(spacePoints.begin(), spacePoints.begin() + back_sps);
      }


    }
    else{
      //Order the spacepoint using the magnitude away from the vertex
      IShowerTool::GetTRACSAlg().OrderShowerSpacePoints(spacePoints,ShowerStartPosition);
    }
    
    if(spacePoints.size() < 3){
      mf::LogError("ShowerResidualTrackHitFinder") << "Not enough spacepoints bailing"<< std::endl;
      return 1;
    }

    double old_residual = 0;


    for(unsigned int sp=0; sp<spacePoints.size(); ++sp){

      //No more spacepoints to do the analysis with bail.
      if(sp+3 > spacePoints.size()){std::cout << "breaking at the start" << std::endl;break;}


      //Add the first 3 spacepoints
      std::vector<art::Ptr<recob::SpacePoint> > spacePoints_fit;
      spacePoints_fit.push_back(spacePoints.at(sp));
      spacePoints_fit.push_back(spacePoints.at(sp+1));
      spacePoints_fit.push_back(spacePoints.at(sp+2));
      sp += 2;

      //Calculate PCA
      TVector3 Eigenvector = ShowerPCAVector(spacePoints_fit,fmh,ShowerStartPosition);

      //Calculate  residual
      double residual = CalculateResidual(spacePoints_fit,Eigenvector,ShowerStartPosition);

      //We are starting at a bad fit
      if((residual-old_residual) > fMaxResidualDiff){

	if(sp+1 == spacePoints.size()){
	  std::cout << "too few sps" << std::endl;
	  //we are somehow at the last point anyways so finish up
	  break;
	}
	
	//Check to see if its an odd one out.
	spacePoints_fit.pop_back(); 
	spacePoints_fit.push_back(spacePoints.at(sp+1)); 
	TVector3 Eigenvector = ShowerPCAVector(spacePoints_fit,fmh,ShowerStartPosition);
	double residual = CalculateResidual(spacePoints_fit,Eigenvector,ShowerStartPosition);
	std::cout << "next residual: " << residual << std::endl;
	if(residual- old_residual > fMaxResidualDiff){
	  //Move on we have finished the vast traul of spacepoint exists to the end of trackstub (hopefully).
	  std::cout << "breaking" << std::endl;
	  break;
	}
      }

      //Add a point recalculate
      while(sp<spacePoints.size()){
      
	//Add the next point to consider 
	spacePoints_fit.push_back(spacePoints.at(sp)); 

	//Calculate PCA
	TVector3 Eigenvector = ShowerPCAVector(spacePoints_fit,fmh,ShowerStartPosition);
      
	//Calculate  residual
	double residual = CalculateResidual(spacePoints_fit,Eigenvector,ShowerStartPosition);

	std::cout << "residual: " << residual << std::endl;

	//Check to see if the residual has blown up 
	if((residual-old_residual) > fMaxResidualDiff){
	
	  if(sp+1 == spacePoints.size()){
	    //we are somehow at the last point anyways so finish up
	    spacePoints_fit.pop_back();
	    break;
	  }

	  //Check to see if its an odd one out.
	  spacePoints_fit.pop_back(); 
	  spacePoints_fit.push_back(spacePoints.at(sp+1)); 
	  TVector3 Eigenvector = ShowerPCAVector(spacePoints_fit,fmh,ShowerStartPosition);
	  double residual = CalculateResidual(spacePoints_fit,Eigenvector,ShowerStartPosition);
	  if(residual-old_residual > fMaxResidualDiff){
	    //Move on to a new greater faster, more productive fit.
	    spacePoints_fit.pop_back();
	    break;
	  }
	  else{
	    //Remove the element and contine
	    spacePoints_fit.pop_back();
	  }
	}

	old_residual = residual;
	++sp;
      }
      //Add to the track spacepoints
      std::cout << "adding hit" << std::endl;
      track_sps.insert(track_sps.end(), spacePoints_fit.begin(), spacePoints_fit.end());
    }
    
    // Get the hits associated to the space points and seperate them by planes
    std::vector<art::Ptr<recob::Hit> > trackHits;
    for(auto const& spacePoint: track_sps){
      std::vector<art::Ptr<recob::Hit> > hits = fmh.at(spacePoint.key());
      for(auto const& hit: hits){
	trackHits.push_back(hit);
      }
    }

    //Add to the holderer
    ShowerEleHolder.SetElement(trackHits, "InitialTrackHits");
    ShowerEleHolder.SetElement(track_sps,"InitialTrackSpacePoints");
    return 0;
  }

  //Function to calculate the shower direction using a charge weight 3D PCA calculation.
    TVector3 ShowerResidualTrackHitFinder::ShowerPCAVector(std::vector<art::Ptr<recob::SpacePoint> >& sps, art::FindManyP<recob::Hit>& fmh, TVector3& ShowerStartPosition){

    //Initialise the the PCA.
    TPrincipal *pca = new TPrincipal(3,"");

    float TotalCharge = 0;

    //Normalise the spacepoints, charge weight and add to the PCA.
    for(auto& sp: sps){

      TVector3 sp_position = IShowerTool::GetTRACSAlg().SpacePointPosition(sp);

      float wht = 1;

      //Normalise the spacepoint position.
      sp_position = sp_position - ShowerStartPosition;

      if(fChargeWeighted){

        //Get the charge.
        float Charge = IShowerTool::GetTRACSAlg().SpacePointCharge(sp,fmh);
	std::cout << "Charge: " << Charge << std::endl;

        //Get the time of the spacepoint
        float Time = IShowerTool::GetTRACSAlg().SpacePointTime(sp,fmh);

        //Correct for the lifetime at the moment.
        Charge *= TMath::Exp((fDetProp->SamplingRate() * Time ) / (fDetProp->ElectronLifetime()*1e3));
	std::cout << "Charge: "<< Charge << std::endl;

        //Charge Weight
        wht *= TMath::Sqrt(Charge/TotalCharge);
      }
      
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
    const TMatrixD* Eigenvectors = pca->GetEigenVectors();

    TVector3 Eigenvector = { (*Eigenvectors)[0][0], (*Eigenvectors)[1][0], (*Eigenvectors)[2][0] };

    return Eigenvector;
  }
  
  
  double ShowerResidualTrackHitFinder::CalculateResidual(std::vector<art::Ptr<recob::SpacePoint> >& sps, TVector3& PCAEigenvector, TVector3& ShowerStartPosition){
    
    double Residual = 0;
    
    for(auto const& sp: sps){
      
      //Get the relative position of the spacepoint
      TVector3 pos= IShowerTool::GetTRACSAlg().SpacePointPosition(sp) - ShowerStartPosition; 
      
      //Gen the perpendicular distance 
      double len  = pos.Dot(PCAEigenvector);
      double perp = (pos - len*PCAEigenvector).Mag();
      
      Residual += perp;
      
    }
    return Residual;
  } 
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerResidualTrackHitFinder)

