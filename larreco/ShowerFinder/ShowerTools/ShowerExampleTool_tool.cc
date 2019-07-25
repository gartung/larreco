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
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"

namespace ShowerRecoTools {

  
  class ShowerExampleTool:IShowerTool {
    
  public:
    
    ShowerExampleTool(const fhicl::ParameterSet& pset);
    
    ~ShowerExampleTool(); 
    
    //Example Direction Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
			  art::Event& Event,
			  reco::shower::ShowerPropertyHolder& ShowerPropHolder
			  ) override;

  private:
    
    // Define standard art tool interface
    void configure(const fhicl::ParameterSet& pset) override;

    //Function to initialise the producer i.e produces<std::vector<recob::Vertex> >(); commands go here.
    void InitialiseProducers() override;

    //Function to add the assoctions
    int AddAssociations(art::Event& Event,
			reco::shower::ShowerElementHolder& ShowerEleHolder) override;

    //prehaps you want a fcl parameter. 
    art::InputTag fPFParticleModuleLabel;
    
    //Maybe an alg
    shower::SBNShowerAlg fSBNShowerAlg;


  };
  
  
  ShowerExampleTool::ShowerExampleTool(const fhicl::ParameterSet& pset)
    //Setup the algs and others here 
  : fSBNShowerAlg(pset.get<fhicl::ParameterSet>("SBNShowerAlg"))

  {
    configure(pset);
  }
  
  ShowerExampleTool::~ShowerExampleTool()
  {
  }
  
  void ShowerExampleTool::configure(const fhicl::ParameterSet& pset)
  {
    //Set up fcl params here 
    fPFParticleModuleLabel      = pset.get<art::InputTag>("PFParticleModuleLabel","");

    
  }

  void ShowerStartPosition::InitialiseProducers(){
    if(producerPtr == NULL){
      mf::LogWarning("ShowerStartPosition") << "The producer ptr has not been set" << std::endl;
      return;
    }

    //Do you create something and you want to save it the event. Initialsie here. For every event with have a vector of showers so each one has a vertex. This is what we are saving. Make sure to use the name "myvertex" later down the line.
    InitialiseProduct<std::vector<recob::Vertex> >("myvertex"); 
  
    //We can also do associations
    InitialiseProduct<art::Assns<recob::Shower, recob::Vertex> >("myvertexassan");


  }


  int ShowerExampleTool::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
					   art::Event& Event,
					   reco::shower::ShowerElementHolder& ShowerEleHolder){

    //In here calculate a shower or element (or multiple). It can be something used to create the recob::shower i.e. the direction. These have specific names so be careful to make these correctly. Alternative you can create something completely new e.g. recob::Vertex and add it the shower element holder

    //Now we are calculating the property of the shower like pfparticle. You have access to everything in the event. Maybe you want the vertex.
    art::Handle<std::vector<recob::Vertex> > vtxHandle;
    std::vector<art::Ptr<recob::Vertex> > vertices;
    if (Event.getByLabel(fPFParticleModuleLabel, vtxHandle))
      art::fill_ptr_vector(vertices, vtxHandle);
    else {
      throw cet::exception("ShowerStartPosition") << "Could not get the pandora vertices. Something is not configured correctly. Please give the correct pandora module label. Stopping"; 
      return 1;
    } 

    

    //Remember the module goes through the tools and if you want to (fcl param) it will loop over them twice. You can check to see if a element has been set with a specific name:
    bool shower_direction_set ShowerElementHolder.CheckElement("ShowerDirection");
    
    TVector3 ShowerDirection = {-999, -999, -999};
    
    //Then you can go and get that element if you want to use it and fill it in for you. 
    if(shower_direction_set){
      ShowerEleHolder.GetElement("ShowerDirection",ShowerDirection);
    }

    //Do some crazy physics - Some legacy code in here for ease. 
    recob::Vertex proposed_vertex = vertices[0]; 
    double xyz[3] = {-999,-999,-999};
    proposed_vertex->XYZ(xyz);

    if(ShowerDirection.X() < 0){
      xyz[0] = - xyz[0];
      xyz[1] = - xyz[1];
      xyz[2] = - xyz[2];
    }
    recob::Vertex vew_vertex = recob::Vertex(xyz);
    TVector3 recobshower_vertex = {xyz[0], xyz[1], xyz[2]};
    TVector3 recobshower_err = {xyz[0]*0.1, xyz[1]*0.1, xyz[2]*0.1}; 
    //You can set elements of the recob::shower just choose the right name (you can acess them later);
    ShowerEleHolder.SetElement(recobshower_verex,recobshower_err,"ShowerStartPosition");
    
    //Or you can set one of the save elements 
    ShowerEleHolder.SetElement(vew_vertex,"myvertex");

    //Or a new unsave one 
    ShowerEleHolder.SetElement(xyz,"xyz");

    //You can also read out what ptr are set and what elements are set: Coming soon.

    return 0;
  }

  void ShowerStartPosition::AddAssociations(art::Event& Event,
					   reco::shower::ShowerElementHolder& ShowerEleHolder
					   ){
    //Here you add elements to associations defined. You can get the art::Ptrs by  GetProducedElementPtr<T>. Then you can add single like a usally association using AddSingle<assn<T>. Ass below.
    const art::Ptr<recob::Vertex> vertexptr =  GetProducedElementPtr<recob::Vertex>("myvertex", ShowerEleHolder);
    const art::Ptr<recob::Shower> showerptr = GetProducedElementPtr<recob::Shower>("shower", ShowerEleHolder);
    AddSingle<art::Assns<recob::Shower, recob::Vertex> >(showerptr,vertexptr,"myvertexassan");
    
    return 0;
  }
  
DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerExampleTool)
  
