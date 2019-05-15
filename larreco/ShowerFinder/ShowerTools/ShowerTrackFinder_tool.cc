//############################################################################
//### Name:        ShowerTrackFinder                                       ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the initial shower track using a rms   ###
//###              based method to define when the shower starts to        ###
//###              shower. This methd is derived from the EMShower_module  ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerInitialTrackFinder.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"

//LArSoft Includes 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"

//C++ Includes 
#include <iostream>

namespace ShowerRecoTools{

  class ShowerTrackFinder:IShowerInitialTrackFinder {
  public:

    ShowerTrackFinder(const fhicl::ParameterSet& pset);
    
    ~ShowerTrackFinder(); 
    
    //Generic Direction Finder
    recob::Track findInitialTrack(const art::Ptr<recob::PFParticle>& pfparticle) override;
  
  private:  

    // Define standard art tool interface
    void configure(const fhicl::ParameterSet& pset) override;
    
  };
  
  
  ShowerTrackFinder::ShowerTrackFinder(const fhicl::ParameterSet& pset)
  {
    configure(pset);
  }
  
  ShowerTrackFinder::~ShowerTrackFinder()
  {
  }
  
  void ShowerTrackFinder::configure(const fhicl::ParameterSet& pset)
  {
    
  }
  
  recob::Track ShowerTrackFinder::findInitialTrack(const art::Ptr<recob::PFParticle>& pfparticle){
    std::cout << "hello world find track" << std::endl;
   
    recob::Track track;
    return track;
  }
  
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackFinder)
