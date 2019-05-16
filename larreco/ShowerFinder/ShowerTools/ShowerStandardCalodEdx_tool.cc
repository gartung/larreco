//############################################################################
//### Name:        ShowerStandardCalodEdx                                  ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the dEdx of the start track of the     ###
//###              shower using the standard calomitry module. Derived     ###
//###              from the EMShower_module.cc                             ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerdEdxFinder.h"

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
#include <vector> 

//Root Includes 
#include "TVector3.h"

namespace ShowerRecoTools{

  class ShowerStandardCalodEdx:IShowerdEdxFinder {
    
  public:

    ShowerStandardCalodEdx(const fhicl::ParameterSet& pset);

    ~ShowerStandardCalodEdx(); 
    
    //Generic Direction Finder
    std::vector<double> finddEdx(const art::Ptr<recob::PFParticle>& pfparticle) override;

  private:
    
    // Define standard art tool interface
    void configure(const fhicl::ParameterSet& pset) override;
    
    
  };
  
  
  ShowerStandardCalodEdx::ShowerStandardCalodEdx(const fhicl::ParameterSet& pset)
  {
    configure(pset);
  }
  
  ShowerStandardCalodEdx::~ShowerStandardCalodEdx()
  {
  }
  
  void ShowerStandardCalodEdx::configure(const fhicl::ParameterSet& pset)
  {
    
  }
  
  std::vector<double> ShowerStandardCalodEdx::finddEdx(const art::Ptr<recob::PFParticle>& pfparticle){
    std::cout << "hello world dEdx" << std::endl;
    std::vector<double>  ShowerStandardCalodEdx = {0,0,0};
    return ShowerStandardCalodEdx;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerStandardCalodEdx)

