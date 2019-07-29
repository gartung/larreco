//############################################################################
//### Name:        ShowerStartPosition                                     ###
//### Author:      Dominic Barker (dominic.barker@sheffield.ac.uk          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the start poistion                     ###
//###              methods.                                                ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "fhiclcpp/ParameterSet.h"  
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
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
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/RecoAlg/SBNShowerAlg.h"


//C++ Includes 
#include <iostream>
#include <vector>

//Root Includes 
#include "TVector3.h"
#include "TMath.h"

namespace ShowerRecoTools{


  class ShowerStartPosition: public IShowerTool {

  public:
    
    ShowerStartPosition(const fhicl::ParameterSet& pset);
    
    ~ShowerStartPosition(); 
    
    // Define standard art tool interface
    void configure(const fhicl::ParameterSet& pset) override;
    
    //Generic Direction Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
			  art::Event& Event,
			  reco::shower::ShowerElementHolder& ShowerEleHolder
			  ) override;
    
    art::InputTag fPFParticleModuleLabel;
    
  private:
    shower::SBNShowerAlg fSBNShowerAlg;

  };
  
  
  ShowerStartPosition::ShowerStartPosition(const fhicl::ParameterSet& pset)
    : fSBNShowerAlg(pset.get<fhicl::ParameterSet>("SBNShowerAlg"))
  {
    configure(pset);
  }
  
  ShowerStartPosition::~ShowerStartPosition()
  {
  }
  
  void ShowerStartPosition::configure(const fhicl::ParameterSet& pset)
  {
    fPFParticleModuleLabel      = pset.get<art::InputTag>("PFParticleModuleLabel","");
  }


  int ShowerStartPosition::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
					    art::Event& Event,
					    reco::shower::ShowerElementHolder& ShowerEleHolder
					    ){

    std::cout << "hello world start position" << std::endl;

    //Get the vertices.
    art::Handle<std::vector<recob::Vertex> > vtxHandle;
    std::vector<art::Ptr<recob::Vertex> > vertices;
    if (Event.getByLabel(fPFParticleModuleLabel, vtxHandle))
      art::fill_ptr_vector(vertices, vtxHandle);
    else {
      throw cet::exception("ShowerStartPosition") << "Could not get the pandora vertices. Something is not configured correctly. Please give the correct pandora module label. Stopping"; 
      return 1;
    } 

    //Get the spacepoints handle and the hit assoication
    art::Handle<std::vector<recob::SpacePoint> > spHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, spHandle)){
      throw cet::exception("ShowerStartPosition") << "Coquld not configure the spacepoint handle. Something is configured incorrectly. Stopping"; 
      return 1;
    } 
    art::FindManyP<recob::Hit> fmh(spHandle, Event, fPFParticleModuleLabel);
    if(!fmh.isValid()){
      throw cet::exception("ShowerStartPosition") << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }

    // Get the assocated pfParicle vertex PFParticles
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("ShowerStartPosition") << "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }
    art::FindManyP<recob::Vertex> fmv(pfpHandle, Event, fPFParticleModuleLabel);
    if(!fmv.isValid()){
      throw cet::exception("ShowerStartPosition") << "Vertex and PF particle association is somehow not valid. Stopping";
      return 1;
    }
    std::vector<art::Ptr<recob::Vertex> > vtx_cand = fmv.at(pfparticle.key());

    //If there is more than one then fail becuase I don't think that this can be the case 
    if(vtx_cand.size() > 1){
	throw cet::exception("ShowerStartPosition") << "There is more than one vertex associated with the pfparticle. Something the creator of this tool did not consider. Stopping";
      return 1;
    }
    
    //If there is only one vertex good news we just say that is the start of the shower.
    if(vtx_cand.size() == 1){
      art::Ptr<recob::Vertex> StartPositionVertex = vtx_cand[0];
      double xyz[3] = {-999,-999,-999};
      StartPositionVertex->XYZ(xyz);
      TVector3 ShowerStartPosition = {xyz[0], xyz[1], xyz[2]};
      TVector3 ShowerStartPositionErr = {-999, -999, -999};
      ShowerEleHolder.SetElement(ShowerStartPosition,ShowerStartPositionErr,"ShowerStartPosition");
      return 0;
    }

    //If we there have none then use the direction to find the neutrino vertex 
    if(ShowerEleHolder.CheckElement("ShowerDirection")){

      TVector3 ShowerDirection = {-999, -999, -999};
      ShowerEleHolder.GetElement("ShowerDirection",ShowerDirection);
      
      art::FindManyP<recob::SpacePoint> fmspp(pfpHandle, Event, fPFParticleModuleLabel);

      if (!fmspp.isValid()){
	throw cet::exception("ShowerStartPosition") << "Trying to get the spacepoints and failed. Something is not configured correctly. Stopping ";
	return 1;
      }

      //Get the spacepoints
      std::vector<art::Ptr<recob::SpacePoint> > spacePoints_pfp = fmspp.at(pfparticle.key());

      //Cannot continue if we have no spacepoints
      if(spacePoints_pfp.size() == 0){return 0;}

      //Get the Shower Center 
      TVector3 ShowerCentre = fSBNShowerAlg.ShowerCentre(spacePoints_pfp,fmh);

      //Order the Hits from the shower centre. The most negative will be the start position.
      fSBNShowerAlg.OrderShowerSpacePoints(spacePoints_pfp,ShowerCentre,ShowerDirection);

      //Set the start position.
      TVector3 ShowerStartPosition = fSBNShowerAlg.SpacePointPosition(spacePoints_pfp[0]);
      TVector3 ShowerStartPositionErr = {-999,-999,-999};
      ShowerEleHolder.SetElement(ShowerStartPosition,ShowerStartPositionErr,"ShowerStartPosition");
      return 0; 
    }
    
    mf::LogWarning("ShowerStartPosition") << "Start Position has not been set yet. If you are not calculating the start position again then maybe you should stop";
    return 0; 
  }


}
DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerStartPosition)

