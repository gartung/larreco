//############################################################################
//### Name:        ShowerGenericTool                                       ###
//### Author:      You                                                     ###
//### Date:        13.05.19                                                ###
//### Description: Generic form of the shower tools                        ###
//###                                                                      ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/PFParticle.h"

namespace ShowerRecoTools {

  
  class ShowerGenericTool:IShowerTool {
    
  public:
    
    ShowerGenericTool(const fhicl::ParameterSet& pset);
    
    ~ShowerGenericTool(); 
    
    //Generic Direction Finder
    int CalculateProperty(const art::Ptr<recob::PFParticle>& pfparticle,
			  art::Event& Event,
			  reco::shower::ShowerPropertyHolder& ShowerPropHolder
			  ) override;

  private:
    
    // Define standard art tool interface
    void configure(const fhicl::ParameterSet& pset) override;


  };
  
  
  ShowerGenericTool::ShowerGenericTool(const fhicl::ParameterSet& pset)
  {
    configure(pset);
  }
  
  ShowerGenericTool::~ShowerGenericTool()
  {
  }
  
  void ShowerGenericTool::configure(const fhicl::ParameterSet& pset)
  {
  }

  int ShowerGenericTool::CalculateProperty(const art::Ptr<recob::PFParticle>& pfparticle,
					   art::Event& Event,
					   reco::shower::ShowerPropertyHolder& ShowerPropHolder){

    return 0;
  }
}
  
DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerGenericTool)
  
