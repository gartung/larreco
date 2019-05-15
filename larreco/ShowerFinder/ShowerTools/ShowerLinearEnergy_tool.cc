//############################################################################
//### Name:        ShowerLinearEnergy                                      ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the Energy of the shower. Derived      ###
//###              from the linear energy algorithm, written for           ###
//###              the EMShower_module.cc                                  ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerEnergyFinder.h"

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

namespace ShowerRecoTools {

  class ShowerLinearEnergy:IShowerEnergyFinder {
    
  public:

    ShowerLinearEnergy(const fhicl::ParameterSet& pset);
    
    ~ShowerLinearEnergy(); 
    
    //Generic Direction Finder
    std::vector<double> findEnergy(const art::Ptr<recob::PFParticle>& pfparticle) override;
    
  private:
    
    // Define standard art tool interface
    void configure(const fhicl::ParameterSet& pset) override;
    
  };
  
  
  ShowerLinearEnergy::ShowerLinearEnergy(const fhicl::ParameterSet& pset)
  {
    configure(pset);
  }
  
  ShowerLinearEnergy::~ShowerLinearEnergy()
  {
  }
  
  void ShowerLinearEnergy::configure(const fhicl::ParameterSet& pset)
  {
    
  }
  
  std::vector<double> ShowerLinearEnergy::findEnergy(const art::Ptr<recob::PFParticle>& pfparticle){
    std::cout << "hello world linear energy" << std::endl;
    std::vector<double>  ShowerLinearEnergy = {0,0,0};
    return ShowerLinearEnergy;
  }
  
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerLinearEnergy)
  

