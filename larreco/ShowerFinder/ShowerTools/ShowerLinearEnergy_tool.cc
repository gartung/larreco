//############################################################################
//### Name:        ShowerLinearEnergy                                      ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the Energy of the shower. Derived      ###
//###              from the linear energy algorithm, written for           ###
//###              the EMShower_module.cc                                  ###
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
#include <vector> 

//Root Includes 
#include "TVector3.h"

namespace ShowerRecoTools {

  class ShowerLinearEnergy:IShowerTool {
    
  public:

    ShowerLinearEnergy(const fhicl::ParameterSet& pset);
    
    ~ShowerLinearEnergy(); 
    
    //Generic Direction Finder
    int findMetric(const art::Ptr<recob::PFParticle>& pfparticle,
		    art::Event& Event,
		    reco::shower::ShowerPropertyHolder& ShowerPropHolder
		    ) override;
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
  
  int ShowerLinearEnergy::findMetric(const art::Ptr<recob::PFParticle>& pfparticle,
				     art::Event& Event,
				     reco::shower::ShowerPropertyHolder& ShowerPropHolder
				     ){
    std::cout << "hello world linear energy" << std::endl;
    std::vector<double>  ShowerLinearEnergy = {0,0,0};
    ShowerPropHolder.SetShowerEnergy(ShowerLinearEnergy);
    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerLinearEnergy)
  

