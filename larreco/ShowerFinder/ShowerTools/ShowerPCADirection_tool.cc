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

//LArSoft Includes 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

//C++ Includes 
#include <iostream>

//Root Includes 
#include "TVector3.h"

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
    
  }

  int ShowerPCADirection::findMetric(const art::Ptr<recob::PFParticle>& pfparticle,
				      art::Event& Event,
				      reco::shower::ShowerPropertyHolder& ShowerPropHolder){
    
    std::cout << "hello world shower direction" << std::endl;
    TVector3 ShowerDirection;
    ShowerPropHolder.SetShowerDirection(ShowerDirection);
    return 0;
  }
}
  
DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerPCADirection)
  
