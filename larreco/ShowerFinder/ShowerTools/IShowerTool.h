//############################################################################
//### Name:        IShowerTool                                             ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Generic Tool for finding the shower energy. Used in     ###
//###              SBNShower_Module.cc                                     ###
//############################################################################

#ifndef IShowerTool_H
#define IShowerTool_H

//Framwork Includes
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Event.h"

//LArSoft Includes 
#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/ShowerFinder/ShowerPropertyHolder.h"

#include <string>

namespace ShowerRecoTools{
  class IShowerTool{

  public:
    
    virtual ~IShowerTool() noexcept = default;
      
    //Generic Energy Finder
    virtual int CalculateProperty(const art::Ptr<recob::PFParticle>& pfparticle,
				  art::Event& Event,
				  reco::shower::ShowerPropertyHolder& ShowerPropHolder
				  ) = 0;

  private:
    
    // Define standard art tool interface
    virtual void configure(const fhicl::ParameterSet& pset) = 0;

  };
}

#endif 
