//############################################################################
//### Name:        IShowerInitialTrackFinder                               ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Generic Tool for finding the initial track of the       ###
//###              shower. Used inSBNShower_Module.cc                      ###
//############################################################################

#ifndef IShowerInitialTrackFinder_H
#define IShowerInitialTrackFinder_H

//Framwork Includes
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Event.h"

//LArSoft Includes 
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "larreco/ShowerFinder/ShowerPropertyHolder.h"

//Root Includes 
#include "TVector3.h"

namespace ShowerRecoTools{
  class IShowerInitialTrackFinder{

  public:

    virtual ~IShowerInitialTrackFinder() noexcept = default; 
      
    //Generic Direction Finder
    virtual recob::Track findInitialTrack(const art::Ptr<recob::PFParticle>& pfparticle) = 0 ;

    void InitialiseEvent(art::Event& evt,reco::shower::ShowerPropertyHolder sprop_holder){
      Event                = &evt;
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
