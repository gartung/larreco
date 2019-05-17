//############################################################################
//### Name:        ShowerTrackFinder                                       ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the initial shower track using a rms   ###
//###              based method to define when the shower starts to        ###
//###              shower. This methd is derived from the EMShower_module  ###
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

//LArSoft Includes 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"

//C++ Includes 
#include <iostream>

namespace ShowerRecoTools{

  class ShowerTrackFinder:IShowerTool {
  public:

    ShowerTrackFinder(const fhicl::ParameterSet& pset);
    
    ~ShowerTrackFinder(); 
    
    //Generic Track Finder
    int findMetric(const art::Ptr<recob::PFParticle>& pfparticle,
		   art::Event& Event,
		   reco::shower::ShowerPropertyHolder& ShowerPropHolder
		   ) override;

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
  
  
  int ShowerTrackFinder::findMetric(const art::Ptr<recob::PFParticle>& pfparticle,
					    art::Event& Event,
					    reco::shower::ShowerPropertyHolder& ShowerPropHolder
					    ){
    std::cout << "hello world find track" << std::endl;
    
    recob::Track track;
    std::cout << "hello world find track end1" << std::endl;
    ShowerPropHolder.SetInitialTrack(track); 
    std::cout << "hello world find track end" << std::endl;
    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackFinder)
