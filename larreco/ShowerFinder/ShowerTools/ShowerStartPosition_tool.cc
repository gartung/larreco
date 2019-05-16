//############################################################################
//### Name:        ShowerStartPosition                                     ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the start poistion                     ###
//###              methods.                                                ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerStartPositionFinder.h"

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

namespace ShowerRecoTools{


  class ShowerStartPosition:IShowerStartPositionFinder {

  public:
    
    ShowerStartPosition(const fhicl::ParameterSet& pset);
    
    ~ShowerStartPosition(); 
    
    // Define standard art tool interface
    void configure(const fhicl::ParameterSet& pset) override;
    
    //Generic Direction Finder
    TVector3 findShowerStartPosition(const art::Ptr<recob::PFParticle>& pfparticle) override;
    
    
  };
  
  
  ShowerStartPosition::ShowerStartPosition(const fhicl::ParameterSet& pset)
  {
    configure(pset);
  }
  
  ShowerStartPosition::~ShowerStartPosition()
  {
  }
  
  void ShowerStartPosition::configure(const fhicl::ParameterSet& pset)
  {
    
  }
  
  TVector3 ShowerStartPosition::findShowerStartPosition(const art::Ptr<recob::PFParticle>& pfparticle){
    std::cout << "hello world start position" << std::endl;
    TVector3 ShowerStartPosition;
    return ShowerStartPosition;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerStartPosition)

