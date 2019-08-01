//############################################################################
//### Name:        ShowerLinearEnergyCheater                                ###
//### Author:      Ed Tyley                                                ###
//### Date:        16.07.19                                                ###
//### Description: Cheating tool using truth for shower direction          ###
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
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/RecoAlg/SBNShowerAlg.h"
#include "larreco/RecoAlg/SBNShowerCheatingAlg.h"
#include "larreco/RecoAlg/MCRecoUtils/RecoUtils.h"
#include "larreco/RecoAlg/MCRecoUtils/ShowerUtils.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "nusimdata/SimulationBase/MCTruth.h"

//C++ Includes
#include <iostream>
#include <cmath>

//Root Includes
#include "TVector3.h"
#include "TMath.h"
#include "TPrincipal.h"
#include "TVector.h"
#include "TTree.h"

namespace ShowerRecoTools {

  class ShowerLinearEnergyCheater:IShowerTool {

    public:

      ShowerLinearEnergyCheater(const fhicl::ParameterSet& pset);

      ~ShowerLinearEnergyCheater();

      //Generic Direction Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      // Define standard art tool interface
      void configure(const fhicl::ParameterSet& pset) override;

      double CalculateEnergy(std::vector<art::Ptr<recob::Hit> >& hits, geo::View_t& view);

      //SBN I think should only use U,V and W
      double fUGradient;
      double fUIntercept;
      double fVGradient;
      double fVIntercept;
      double fZGradient;
      double fZIntercept;
      double fXGradient;
      double fXIntercept;
      double fYGradient;
      double fYIntercept;
      double f3DGradient;
      double f3DIntercept;


      art::InputTag fPFParticleModuleLabel;
      art::InputTag fHitModuleLabel;
      detinfo::DetectorProperties const* detprop = nullptr;
      art::ServiceHandle<geo::Geometry> fGeom;

      //Algorithm functions
      shower::SBNShowerAlg         fSBNShowerAlg;
      shower::SBNShowerCheatingAlg fSBNShowerCheatingAlg;
  };


  ShowerLinearEnergyCheater::ShowerLinearEnergyCheater(const fhicl::ParameterSet& pset):
    detprop(lar::providerFrom<detinfo::DetectorPropertiesService>()),
    fSBNShowerAlg(pset.get<fhicl::ParameterSet>("SBNShowerAlg")),
    fSBNShowerCheatingAlg(pset.get<fhicl::ParameterSet>("SBNShowerCheatingAlg"))
  {
    configure(pset);
  }

  ShowerLinearEnergyCheater::~ShowerLinearEnergyCheater()
  {
  }

  void ShowerLinearEnergyCheater::configure(const fhicl::ParameterSet& pset)
  {
    fUGradient   = pset.get<double>("UGradient");
    fUIntercept  = pset.get<double>("UIntercept");
    fVGradient   = pset.get<double>("VGradient");
    fVIntercept  = pset.get<double>("VIntercept");
    fZGradient   = pset.get<double>("ZGradient");
    fZIntercept  = pset.get<double>("ZIntercept");
    fXGradient   = pset.get<double>("XGradient");
    fXIntercept  = pset.get<double>("XIntercept");
    fYGradient   = pset.get<double>("YGradient");
    fYIntercept  = pset.get<double>("YIntercept");
    f3DGradient  = pset.get<double>("ThreeDGradient");
    f3DIntercept = pset.get<double>("ThreeDIntercept");

    fHitModuleLabel        = pset.get<art::InputTag> ("HitModuleLabel");
    fPFParticleModuleLabel = pset.get<art::InputTag> ("PFParticleModuleLabel","");
  }

  int ShowerLinearEnergyCheater::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle, art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){

    std::cout <<"#########################################\n"<<
      "hello world linear energy cheater\n" <<"#########################################\n"<< std::endl;

    const simb::MCParticle* trueParticle;

    //Get the hits from the shower:
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("ShowerLinearEnergy") << "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    if (ShowerEleHolder.CheckElement("TrueParticle")){
      ShowerEleHolder.GetElement("TrueParticle",trueParticle);
    } else {

      //Could store these in the shower element holder and just calculate once?
      std::map<int,const simb::MCParticle*> trueParticles = fSBNShowerCheatingAlg.GetTrueParticleMap();
      std::map<int,std::vector<int> > showersMothers = fSBNShowerCheatingAlg.GetTrueChain(trueParticles);


      //Get the clusters
      art::Handle<std::vector<recob::Cluster> > clusHandle;
      if (!Event.getByLabel(fPFParticleModuleLabel, clusHandle)){
        throw cet::exception("ShowerLinearEnergy") << "Could not get the pandora clusters. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
        return 1;
      }
      art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleModuleLabel);
      std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(pfparticle.key());

      //Get the hit association
      art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleModuleLabel);

      std::vector<art::Ptr<recob::Hit> > showerHits;
      for(auto const& cluster: clusters){

        //Get the hits
        std::vector<art::Ptr<recob::Hit> > hits = fmhc.at(cluster.key());       showerHits.insert(showerHits.end(),hits.begin(),hits.end());
        showerHits.insert(showerHits.end(),hits.begin(),hits.end());
      }

      //Get the true particle from the shower
      std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain(showersMothers,showerHits,2);

      if(ShowerTrackInfo.first==-99999) {
        std::cout<<"True Shower Not Found"<<std::endl;
        return 1;
      }
      trueParticle = trueParticles[ShowerTrackInfo.first];
    }

    //This is all based on the shower vertex being known. If it is not lets not do the track
    if(!ShowerEleHolder.CheckElement("ShowerStartPosition")){
      mf::LogError("ShowerLinearEnergy") << "Start position not set, returning "<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement("ShowerDirection")){
      mf::LogError("ShowerLinearEnergy") << "Direction not set, returning "<< std::endl;
      return 1;
    }
    art::Handle<std::vector<recob::Hit> > hitHandle;
    std::vector<art::Ptr<recob::Hit> > hits;
    if(Event.getByLabel(fHitModuleLabel, hitHandle)){
      art::fill_ptr_vector(hits, hitHandle);
    }



    //Holder for the final product
    std::vector<double> ShowerLinearEnergy;
    unsigned int numPlanes = fGeom->Nplanes();

    std::map<geo::View_t, std::vector<art::Ptr<recob::Hit> > > viewHitMap;
    for (auto hit : hits){
      int trueParticleID = TMath::Abs(RecoUtils::TrueParticleID(hit));
      if (trueParticleID == trueParticle->TrackId()){
        geo::View_t view = hit->View();
        viewHitMap[view].push_back(hit);
      }
    }

    std::map<unsigned int, double > view_energies;
    //Accounting for events crossing the cathode.
    for(auto const& viewHit: viewHitMap){

      geo::View_t view = viewHit.first;
      std::vector<art::Ptr<recob::Hit> > hits = viewHit.second;

      //Calculate the Energy for
      double Energy = CalculateEnergy(hits,view);

      std::cout<<"view: "<<view<<" and energy: "<<Energy<<std::endl;

      unsigned int viewNum = view;
      view_energies[viewNum] = Energy;
    }

    //TODO think of a better way of doing this
    for (unsigned int plane=0; plane<numPlanes; ++plane) {
      //std::cout<<"Plane: "<<plane<<std::endl;
      double Energy;
      try{
        Energy = view_energies.at(plane);
        if (Energy<0){
          mf::LogWarning("ShowerLinearEnergy") << "Negative shower energy: "<<Energy;
          Energy=-999;
        }

      } catch(...){
        mf::LogError("ShowerLinearEnergy") <<"No energy calculation for plane "<<plane<<std::endl;
        Energy = -999;
      }
      ShowerLinearEnergy.push_back(Energy);
    }

    if(ShowerLinearEnergy.size() == 0){
      throw cet::exception("ShowerLinearEnergy") << "Energy Vector is empty";
      return 1;
    }

    //Tod do
    std::vector<double> EnergyError = {-999,-999,-999};

    ShowerEleHolder.SetElement(ShowerLinearEnergy,EnergyError,"ShowerEnergy");
    return 0;

    std::cout <<"#########################################\n"<<
      "linear energy cheater Done\n" <<"#########################################\n"<< std::endl;
    return 0;
  }

  double ShowerLinearEnergyCheater::CalculateEnergy(std::vector<art::Ptr<recob::Hit> >& hits, geo::View_t& view){

    double totalCharge = 0, totalEnergy = 0;

    for (art::PtrVector<recob::Hit>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit){
      totalCharge += ( (*hit)->Integral() * TMath::Exp( (detprop->SamplingRate() * (*hit)->PeakTime()) / (detprop->ElectronLifetime()*1e3) ) );
    }

    switch (view) {
      case geo::kU:
        totalEnergy = (totalCharge * fUGradient) + fUIntercept;
        break;
      case geo::kV:
        totalEnergy = (totalCharge * fVGradient) + fVIntercept;
        break;
      case geo::kW: //same as geo::kZ
        totalEnergy = (totalCharge * fZGradient) + fZIntercept;
        break;
      case geo::kX:
        totalEnergy = (totalCharge * fXGradient) + fXIntercept;
        break;
      case geo::kY:
        totalEnergy = (totalCharge * fYGradient) + fYIntercept;
        break;
      case geo::k3D:
        totalEnergy = (totalCharge * f3DGradient) + f3DIntercept;
        break;
      default:
        throw cet::exception("ShowerLinearEnergy") << "View Not configured";
        return 1;

    }
    return totalEnergy;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerLinearEnergyCheater)
