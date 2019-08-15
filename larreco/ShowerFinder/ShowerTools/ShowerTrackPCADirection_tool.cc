//############################################################################
//### Name:        ShowerTrackPCADirection                                 ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the shower direction using PCA         ###
//###              methods using the initial track hits                    ###
//############################################################################

//Warning! Currently as pandora gives each hit a spacepoint, rather than
//         matching up some energy depositions are double counted.
//         This could lead to a bais in the PCA analysis.

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
#include "lardataobj/RecoBase/Wire.h"
//C++ Includes
#include <iostream>
#include <cmath>

//Root Includes
#include "TVector3.h"
#include "TMath.h"
#include "TPrincipal.h"
#include "TVector.h"

namespace ShowerRecoTools {


  class ShowerTrackPCADirection:IShowerTool {

    public:

      ShowerTrackPCADirection(const fhicl::ParameterSet& pset);

      ~ShowerTrackPCADirection();

      //Generic Direction Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      TVector3 ShowerPCAVector(std::vector<art::Ptr<recob::SpacePoint> >& spacePoints_pfp, art::FindManyP<recob::Hit>& fmh, TVector3& ShowerCentre);

      //Algorithm functions
      shower::SBNShowerAlg fSBNShowerAlg;

      //Services
      detinfo::DetectorProperties const* fDetProp;

      //fcl
      art::InputTag fPFParticleModuleLabel;
      art::InputTag fHitModuleLabel;
      bool fChargeWeighted;
      bool fDebugEVD;
      unsigned int  fMinPCAPoints;
  };


  ShowerTrackPCADirection::ShowerTrackPCADirection(const fhicl::ParameterSet& pset)
    : fSBNShowerAlg(pset.get<fhicl::ParameterSet>("SBNShowerAlg")),
    fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())
  {
    fPFParticleModuleLabel  = pset.get<art::InputTag>("PFParticleModuleLabel","");
    fHitModuleLabel         = pset.get<art::InputTag>("HitModuleLabel");
    fChargeWeighted         = pset.get<bool>         ("ChargeWeighted");
    fDebugEVD               = pset.get<bool>         ("DebugEVD");
    fMinPCAPoints           = pset.get<unsigned int>          ("MinPCAPoints");
  }

  ShowerTrackPCADirection::~ShowerTrackPCADirection()
  {
  }

  int ShowerTrackPCADirection::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event,
      reco::shower::ShowerElementHolder& ShowerEleHolder){

    std::cout <<"#########################################\n"<<
      "hello world Track PCA direction\n" <<"#########################################\n"<< std::endl;
    if(!ShowerEleHolder.CheckElement("ShowerStartPosition")){
      mf::LogError("ShowerTrackPCA") << "Start Position not set. Stopping" << std::endl;;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement("InitialTrackSpacePoints")){
      mf::LogError("ShowerTrackPCA") << "TrackSpacePoints not set, returning "<< std::endl;
      return 1;
    }

    //Get the spacepoints handle and the hit assoication
    art::Handle<std::vector<recob::SpacePoint> > spHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, spHandle)){
      throw cet::exception("ShowerTrackPCA") << "Could not configure the spacepoint handle. Something is configured incorrectly. Stopping";
      return 1;
    }
    art::FindManyP<recob::Hit> fmh(spHandle, Event, fPFParticleModuleLabel);
    if(!fmh.isValid()){
      throw cet::exception("ShowerTrackPCA") << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }

    std::vector<art::Ptr<recob::SpacePoint> > trackSpacePoints;
    ShowerEleHolder.GetElement("InitialTrackSpacePoints",trackSpacePoints);

    std::cout<<"Num initial track spacepoints: "<<trackSpacePoints.size()<<std::endl;

    //We cannot progress with no spacepoints.
    if(trackSpacePoints.size() < fMinPCAPoints){
      mf::LogError("ShowerTrackPCA") << "Not enough spacepoints for PCA, returning "<< std::endl;
      return 1;
    }

    //Find the PCA vector
    TVector3 trackCentre;
    TVector3 Eigenvector = ShowerPCAVector(trackSpacePoints, fmh, trackCentre);


    //Get the General direction as the vector between the start position and the centre
    TVector3 StartPositionVec = {-999, -999, -999};
    ShowerEleHolder.GetElement("ShowerStartPosition",StartPositionVec);
    TVector3 GeneralDir       = (trackCentre - StartPositionVec).Unit();

    //Dot product
    double DotProduct = Eigenvector.Dot(GeneralDir);

    //If the dotproduct is negative the Direction needs Flipping
    if(DotProduct < 0){
      Eigenvector[0] = - Eigenvector[0];
      Eigenvector[1] = - Eigenvector[1];
      Eigenvector[2] = - Eigenvector[2];
    }

    TVector3 EigenvectorErr = {-999,-999,-999};

    ShowerEleHolder.SetElement(Eigenvector,EigenvectorErr,"ShowerDirection");

    if (fDebugEVD){
      std::cout<<"Do DebugEVD"<<std::endl;
      fSBNShowerAlg.DebugEVD(pfparticle,Event,ShowerEleHolder);
    }

    std::cout <<"#########################################\n"<<
      "Track PCA direction Done\n" <<"#########################################\n"<< std::endl;
    return 0;
  }


  //Function to calculate the shower direction using a charge weight 3D PCA calculation.
  TVector3 ShowerTrackPCADirection::ShowerPCAVector(std::vector<art::Ptr<recob::SpacePoint> >& sps, art::FindManyP<recob::Hit>& fmh, TVector3& ShowerCentre){

    //Initialise the the PCA.
    TPrincipal *pca = new TPrincipal(3,"");

    float TotalCharge = 0;

    //Get the Shower Centre
    ShowerCentre = fSBNShowerAlg.ShowerCentre(sps, fmh, TotalCharge);


    //Normalise the spacepoints, charge weight and add to the PCA.
    for(auto& sp: sps){

      TVector3 sp_position = fSBNShowerAlg.SpacePointPosition(sp);

      float wht = 1;

      //Normalise the spacepoint position.
      sp_position = sp_position - ShowerCentre;

      if(fChargeWeighted){

        //Get the charge.
        float Charge = fSBNShowerAlg.SpacePointCharge(sp,fmh);

        //Get the time of the spacepoint
        float Time = fSBNShowerAlg.SpacePointTime(sp,fmh);

        //Correct for the lifetime at the moment.
        Charge *= TMath::Exp((fDetProp->SamplingRate() * Time ) / (fDetProp->ElectronLifetime()*1e3));

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

    // TODO: should be a neater way to get the column
    //TVectorD Principal = TMatrixDColumn(Eigenvectors,0);
    //std::cout<<Principal[0];
    //Get the primary vector
    //const TVectorD EigenvectorD = (*Eigenvectors)[0]; // Gets row not column

    TVector3 Eigenvector = { (*Eigenvectors)[0][0], (*Eigenvectors)[1][0], (*Eigenvectors)[2][0] };

    return Eigenvector;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackPCADirection)

