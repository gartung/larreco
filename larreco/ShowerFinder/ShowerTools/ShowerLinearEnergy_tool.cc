//############################################################################
//### Name:        ShowerLinearEnergy                                      ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the Energy of the shower. Derived      ###
//###              from the linear energy algorithm, written for           ###
//###              the EMShower_module.cc                                  ###
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
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

//C++ Includes 
#include <iostream>
#include <vector> 

//Root Includes 
#include "TVector3.h"

namespace ShowerRecoTools {

  class ShowerLinearEnergy:IShowerTool {
    
  public:

    ShowerLinearEnergy(const fhicl::ParameterSet& pset);
    
    ~ShowerLinearEnergy(); 
    
    //Generic Direction Finder
    int findMetric(const art::Ptr<recob::PFParticle>& pfparticle,
		    art::Event& Event,
		    reco::shower::ShowerPropertyHolder& ShowerPropHolder
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
    detinfo::DetectorProperties const* detprop = nullptr;
  };
  
  
  ShowerLinearEnergy::ShowerLinearEnergy(const fhicl::ParameterSet& pset):
    detprop(lar::providerFrom<detinfo::DetectorPropertiesService>())
  {
    configure(pset);
  }
  
  ShowerLinearEnergy::~ShowerLinearEnergy()
  {    
  }
  
  void ShowerLinearEnergy::configure(const fhicl::ParameterSet& pset)
  {
    
    fPFParticleModuleLabel  = pset.get<art::InputTag>("PFParticleModuleLabel","");
  
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

  }
  
  int ShowerLinearEnergy::findMetric(const art::Ptr<recob::PFParticle>& pfparticle,
				     art::Event& Event,
				     reco::shower::ShowerPropertyHolder& ShowerPropHolder
				     ){

    //Holder for the final product
    std::vector<double> ShowerLinearEnergy;

    // Get the assocated pfParicle vertex PFParticles
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("ShowerLinearEnergy") << "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }
    
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
    
    //Match up views
    std::map<geo::View_t,std::vector<art::Ptr<recob::Hit> > > view_clusters;

    //Loop over the clusters in the plane and get the hits 
    for(auto const& cluster: clusters){
       
      //Get the hits 
      std::vector<art::Ptr<recob::Hit> > hits = fmhc.at(cluster.key());

      //Get the view. 
      geo::View_t view = cluster->View();
 
      view_clusters[view].insert(view_clusters[view].end(),hits.begin(),hits.end());
    }

    //Accounting for events crossing the cathode.
    for(auto const& cluster: view_clusters){

      std::vector<art::Ptr<recob::Hit> > hits = cluster.second;
      geo::View_t view = cluster.first;

      //Calculate the Energy for 
      double Energy = CalculateEnergy(hits,view);
      
      ShowerLinearEnergy.push_back(Energy);
    }
    
    if(ShowerLinearEnergy.size() == 0){
      throw cet::exception("ShowerLinearEnergy") << "Energy Vector is empty";
      return 1;
    }

    ShowerPropHolder.SetShowerEnergy(ShowerLinearEnergy);
    return 0;
  }

  //Function to calculate the energy of a shower in a plane. Using a linear map between charge and Energy.
  //Exactly the same method as the ShowerEnergyAlg.cxx. Thanks Mike. 
  double ShowerLinearEnergy::CalculateEnergy(std::vector<art::Ptr<recob::Hit> >& hits, geo::View_t& view){

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

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerLinearEnergy)
  

