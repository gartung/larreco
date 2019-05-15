//############################################################################
//### Name:        ShowerPCADirection                                      ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the shower direction using PCA         ###
//###              methods. Derviced from PandoraShowers Method            ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerDirectionFinder.h"

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

  
  class ShowerPCADirection:IShowerDirectionFinder {
    
  public:
    
    ShowerPCADirection(const fhicl::ParameterSet& pset);
    
    ~ShowerPCADirection(); 
    
    //Generic Direction Finder
    TVector3 findDirection(const art::Ptr<recob::PFParticle>& pfparticle) override;
    
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
  
  TVector3 ShowerPCADirection::findDirection(const art::Ptr<recob::PFParticle>& pfparticle){
    std::cout << "hello world shower direction" << std::endl;
    TVector3 ShowerDirection;
    return ShowerDirection;
  }
  
}
  
DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerPCADirection)
  
