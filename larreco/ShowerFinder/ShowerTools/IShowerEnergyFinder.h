//############################################################################
//### Name:        IShowerEnergyFinder                                     ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Generic Tool for finding the shower energy. Used in     ###
//###              SBNShower_Module.cc                                     ###
//############################################################################

#ifndef IShowerEnergyFinder_H
#define IShowerEnergyFinder_H

//Framwork Includes
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Event.h"

//LArSoft Includes 
#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/ShowerFinder/ShowerPropertyHolder.h"

//C++ Includes
#include <vector>

namespace ShowerRecoTools{
  class IShowerEnergyFinder{

  public:
    
    virtual ~IShowerEnergyFinder() noexcept = default;
      
    //Generic Energy Finder
    virtual std::vector<double> findEnergy(const art::Ptr<recob::PFParticle>& pfparticle) = 0;

    void InitialiseEvent(art::Event& evt, reco::shower::ShowerPropertyHolder sprop_holder){
      Event            = &evt;
      ShowerPropHolder = &sprop_holder;
    }; 

  private:
    
    // Define standard art tool interface
    virtual void configure(const fhicl::ParameterSet& pset) = 0;

    art::Event* Event;
    reco::shower::ShowerPropertyHolder* ShowerPropHolder;  

  };
}

#endif 
