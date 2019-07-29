//################################################################################
//### Name: SBNShower                                                          ###
//### Author: Dominic Barker                                                   ###
//### Date: 15.05.19                                                           ###
//### Description: Generic Shower Charaterisation module which allows the      ###
//###              the user choose which tool to calculate shower metrics.     ###
//###              For a complete shower the tools must define (use the exact  ###
//###              name) the following  metrics in the shower property holder: ###
//###              ShowerStartPosition                                         ###
//###              ShowerDirection                                             ###
//###              ShowerEnergy                                                ###
//###              ShowerdEdx                                                  ###
//################################################################################

//Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Utilities/make_tool.h" 
#include "art_root_io/TFileDirectory.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Vertex.h"


//LArSoft includes
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"
#include "larreco/ShowerFinder/ShowerElementHolder.hh"
#include "larreco/ShowerFinder/ShowerProduedPtrsHolder.hh"

//Root Includes
#include "TVector3.h"

//C++ Includes 
#include <vector>

namespace reco {
  namespace shower {
    class SBNShower;
  }
}

//Class 

class reco::shower::SBNShower: public art::EDProducer {
public:

  SBNShower(fhicl::ParameterSet const& pset);

  void produce(art::Event& evt);
  void reconfigure(fhicl::ParameterSet const& p);

  template <class T >
  art::Ptr<T> GetProducedElementPtr(std::string InstanceName, reco::shower::ShowerElementHolder& ShowerEleHolder, int iter=-1);


private:

  //fcl object names 
  art::InputTag fPFParticleModuleLabel;
  bool          fSecondInteration;
  bool          fAllowPartialShowers;
  bool          fVerbose; 

  //fcl tools
  std::vector<std::unique_ptr<ShowerRecoTools::IShowerTool> > fShowerTools;
  std::vector<std::string>                                    fShowerToolNames;

  //map to the unique ptrs to
  reco::shower::ShowerProduedPtrsHolder uniqueproducerPtrs;

};

template <class T >
art::Ptr<T> reco::shower::SBNShower::GetProducedElementPtr(std::string InstanceName, reco::shower::ShowerElementHolder& ShowerEleHolder, int iter){
  
  bool check_element = ShowerEleHolder.CheckElement(InstanceName);
  if(!check_element){
    throw cet::exception("SBNShower") << "To get a element that does not exist" << std::endl;
    return art::Ptr<T>();
  }
  
  bool check_ptr = uniqueproducerPtrs.CheckUniqueProduerPtr(InstanceName);
  if(!check_ptr){
    throw cet::exception("SBNShower") << "Tried to get a ptr that does not exist" << std::endl;
    return art::Ptr<T>();
  }
  
  
  //Get the number of the shower we are on.
  int index;
  if(iter != -1){
    index = iter;
  }
  else{
    index = ShowerEleHolder.GetShowerNumber();
  }

  //Make the ptr
  art::Ptr<T> artptr = uniqueproducerPtrs.GetArtPtr<T>(InstanceName,index);
  return artptr;
}



reco::shower::SBNShower::SBNShower(fhicl::ParameterSet const& pset) :
  EDProducer{pset}
{
  this->reconfigure(pset);
  produces<std::vector<recob::Shower> >();
  produces<art::Assns<recob::Shower, recob::Hit> >();
  produces<art::Assns<recob::Shower, recob::Cluster> >();
  produces<art::Assns<recob::Shower, recob::SpacePoint> >();
  produces<art::Assns<recob::Shower, recob::PFParticle> >();

  // Output -- showers and associations with hits and clusters
  uniqueproducerPtrs.SetShowerUniqueProduerPtr(type<std::vector<recob::Shower> >(),"shower");
  uniqueproducerPtrs.SetShowerUniqueProduerPtr(type<art::Assns<recob::Shower, recob::Cluster> >(),"clusterAssociationsbase");
  uniqueproducerPtrs.SetShowerUniqueProduerPtr(type<art::Assns<recob::Shower, recob::Hit> >(),"hitAssociationsbase");
  uniqueproducerPtrs.SetShowerUniqueProduerPtr(type<art::Assns<recob::Shower, recob::SpacePoint > >(),"spShowerAssociationsbase");
  uniqueproducerPtrs.SetShowerUniqueProduerPtr(type<art::Assns<recob::Shower, recob::PFParticle> >(),"pfShowerAssociationsbase");

  uniqueproducerPtrs.PrintPtrs();

}

void reco::shower::SBNShower::reconfigure(fhicl::ParameterSet const& pset) {

  //Intialise the tools 
  auto const tool_psets = pset.get<std::vector<fhicl::ParameterSet>>("ShowerFinderTools");
  for (auto const& tool_pset : tool_psets) {
    fShowerTools.push_back(art::make_tool<ShowerRecoTools::IShowerTool>(tool_pset));
    std::string paramset = tool_pset.to_compact_string();
    std::size_t pos = paramset.find("tool_type:");
    fShowerToolNames.push_back(paramset.substr(pos+10));
    std::cout << "Tools List: " << paramset.substr(pos) << std::endl;
  }

  //  Initialise the EDProducer ptr in the tools 
  for(auto const& ShowerTool: fShowerTools){
    ShowerTool->SetPtr(this);
    ShowerTool->InitaliseProducerPtr(uniqueproducerPtrs);
    ShowerTool->InitialiseProducers();
  }

  //Initialise the other paramters.
  fPFParticleModuleLabel = pset.get<art::InputTag>("PFParticleModuleLabel","pandora");
  fSecondInteration      = pset.get<bool         >("SecondInteration",false);
  fAllowPartialShowers   = pset.get<bool         >("AllowPartialShowers",false);
  fVerbose               = pset.get<bool         >("Verbose",false);
}
    
void reco::shower::SBNShower::produce(art::Event& evt) {

  //Ptr makers for the products 
  uniqueproducerPtrs.SetPtrMakers(evt);

  //Get the PFParticles
  art::Handle<std::vector<recob::PFParticle> > pfpHandle;
  std::vector<art::Ptr<recob::PFParticle> > pfps;
  if (evt.getByLabel(fPFParticleModuleLabel, pfpHandle)){
    art::fill_ptr_vector(pfps, pfpHandle);
  }
  else {
    throw cet::exception("SBNShower") << "pfps not loaded." << std::endl;
  }
 
  //Handle to access the pandora hits assans
  art::Handle<std::vector<recob::Cluster> > clusterHandle;
  if (!evt.getByLabel(fPFParticleModuleLabel,clusterHandle)){
    throw cet::exception("SBNShower") << "pfp clusters is not loaded." << std::endl;
  }

  //Get the assoications to hits, clusters and spacespoints 
  art::FindManyP<recob::Hit> fmh(clusterHandle, evt, fPFParticleModuleLabel);
  art::FindManyP<recob::Cluster> fmcp(pfpHandle, evt, fPFParticleModuleLabel);
  art::FindManyP<recob::SpacePoint> fmspp(pfpHandle, evt, fPFParticleModuleLabel);

  if(!fmcp.isValid()){
    throw cet::exception("SBNShower") << "Find many clusters is not valid." << std::endl;
  }
  if(!fmh.isValid()){
    throw cet::exception("SBNShower") << "Find many hits is not valid." << std::endl;
  }
  if(!fmspp.isValid()){
    throw cet::exception("SBNShower") << "Find many spacepoints is not valid." << std::endl;
  }

  //Holder to pass to the functions, contains the 6 properties of the shower 
  // - Start Poistion
  // - Direction
  // - Initial Track
  // - Initial Track Hits
  // - Energy 
  // - dEdx 
  reco::shower::ShowerElementHolder selement_holder;

  int i=0;
  int shower_iter = 0;
  //Loop of the pf particles
  for(auto const& pfp: pfps){

    std::cout << "new particle" << std::endl;
    //Update the shower iterator
    selement_holder.SetShowerNumber(shower_iter);

    //loop only over showers.
    if(pfp->PdgCode() != 11){continue;}

    //Calculate the shower properties 
    //Loop over the shower tools
    int err = 0;
    for(auto const& fShowerTool: fShowerTools){

      std::cout << "on next tool:" <<  fShowerToolNames[i]  << std::endl;
      //Calculate the metric
      err = fShowerTool->CalculateElement(pfp,evt,selement_holder);
      if(err){
	mf::LogError("SBNShower") << "Error in shower tool: " << fShowerToolNames[i]  << " with code: " << err << std::endl;
	if(!fAllowPartialShowers && !fSecondInteration) break;
      }
      ++i;
    }

    //Should we do a second interaction now we have done a first pass of the calculation
    i=0;
    if(fSecondInteration){
      std::cout << "on second iter next tool:" <<  fShowerToolNames[i]  << std::endl;

      for(auto const& fShowerTool: fShowerTools){
	//Calculate the metric
	err = fShowerTool->CalculateElement(pfp,evt,selement_holder);
      
	if(err){
	  mf::LogError("SBNShower") << "Error in shower tool: " << fShowerToolNames[i]  << " with code: " << err << std::endl;
	  if(!fAllowPartialShowers) break;
	}
	++i;
      }
    }

    //If we want a full shower and we recieved an error call from a tool return;
    if(err && !fAllowPartialShowers){
      mf::LogError("SBNShower") << "Error on tool. Assuming all the shower products and properties were not setting and bailing." << std::endl;
      continue;
    }

    //If we are are not allowing partial shower check all the products to make the shower are correctly set
    if(!fAllowPartialShowers){
      bool accept = true; 
      accept = accept & selement_holder.CheckElement("ShowerStartPosition");
      accept = accept & selement_holder.CheckElement("ShowerDirection");
      accept = accept & selement_holder.CheckElement("ShowerEnergy");
      accept = accept & selement_holder.CheckElement("ShowerdEdx");
      if(!accept){continue;}
          
      //Check All of the products that have been asked to be checked.
      bool elements_are_set = selement_holder.CheckAllElementSaveTag();
      if(!elements_are_set){
	mf::LogError("SBNShower") << "Not all the elements in the property holder which should be set are not. Bailing. " << std::endl; 
	continue;
      }
      
      ///Check all the producers 
      bool producers_are_set = uniqueproducerPtrs.CheckAllProducedElements(selement_holder);
      if(!producers_are_set){
	mf::LogError("SBNShower") << "Not all the elements in the property holder which are produced are not set. Bailing. " << std::endl; 
	continue;
      }
    }

    //Get the properties 
    TVector3                           ShowerStartPosition  = {-999,-999,-999};
    TVector3                           ShowerDirection      = {-999,-999,-999};
    std::vector<double>                ShowerEnergy         = {-999,-999,-999};
    std::vector<double>                ShowerdEdx           = {-999,-999,-999};

    int                                BestPlane               = -999;
    TVector3                           ShowerStartPositionErr  = {-999,-999,-999};
    TVector3                           ShowerDirectionErr      = {-999,-999,-999};
    std::vector<double>                ShowerEnergyErr         = {-999,-999,-999};
    std::vector<double>                ShowerdEdxErr           = {-999,-999,-999};
    
    err = 0;
    if(selement_holder.CheckElement("ShowerStartPosition"))    err += selement_holder.GetElementAndError("ShowerStartPosition",ShowerStartPosition,ShowerStartPositionErr);
    if(selement_holder.CheckElement("ShowerDirection"))        err += selement_holder.GetElementAndError("ShowerDirection",ShowerDirection,ShowerDirectionErr);
    if(selement_holder.CheckElement("ShowerEnergy"))           err += selement_holder.GetElementAndError("ShowerEnergy",ShowerEnergy,ShowerEnergyErr);
    if(selement_holder.CheckElement("ShowerdEdx"))             err += selement_holder.GetElementAndError("ShowerdEdx",ShowerdEdx,ShowerdEdxErr  );
    if(selement_holder.CheckElement("ShowerBestPlane"))        err += selement_holder.GetElement("ShowerBestPlane",BestPlane);

    if(err){
      throw cet::exception("SBNShower")  << "Error in SBNShower Module. A Check on a shower property failed " << std::endl;
    }

    if(fVerbose){
      //Check the shower
      std::cout<<"Shower Vertex: X:"<<ShowerStartPosition.X()<<" Y: "<<ShowerStartPosition.Y()<<" Z: "<<ShowerStartPosition.Z()<<std::endl;
      std::cout<<"Shower Direction: X:"<<ShowerDirection.X()<<" Y: "<<ShowerDirection.Y()<<" Z: "<<ShowerDirection.Z()<<std::endl;
      std::cout<<"Shower dEdx: size: "<<ShowerdEdx.size()<<" Plane 0: "<<ShowerdEdx.at(0)<<" Plane 1: "<<ShowerdEdx.at(1)<<" Plane 2: "<<ShowerdEdx.at(2)<<std::endl;
      std::cout<<"Shower Energy: size: "<<ShowerEnergy.size()<<" Plane 0: "<<ShowerEnergy.at(0)<<" Plane 1: "<<ShowerEnergy.at(1)<<" Plane 2: "<<ShowerEnergy.at(2)<<std::endl;
      std::cout<<"Shower Best Plane: "<<BestPlane<<std::endl;

      //Print what has been created in the shower
      selement_holder.PrintElements();
    }

    //Make the shower 
    recob::Shower shower = recob::Shower(ShowerDirection, ShowerDirectionErr,ShowerStartPosition, ShowerDirectionErr,ShowerEnergy,ShowerEnergyErr,ShowerdEdx, ShowerdEdxErr, BestPlane, -999);
    selement_holder.SetElement(shower,"shower");
    ++shower_iter;
    art::Ptr<recob::Shower> ShowerPtr = this->GetProducedElementPtr<recob::Shower>("shower",selement_holder);

    //Associate the pfparticle 
    uniqueproducerPtrs.AddSingle<art::Assns<recob::Shower, recob::PFParticle>>(ShowerPtr,pfp,"pfShowerAssociationsbase");
        
    //Get the associated hits,clusters and spacepoints
    std::vector<art::Ptr<recob::Cluster> >    showerClusters    = fmcp.at(pfp.key());
    std::vector<art::Ptr<recob::SpacePoint> > showerSpacePoints = fmspp.at(pfp.key());

    //Add the hits for each "cluster"
    for(auto const& cluster: showerClusters){

      //Associate the clusters 
      std::vector<art::Ptr<recob::Hit> > ClusterHits = fmh.at(cluster.key());
      uniqueproducerPtrs.AddSingle<art::Assns<recob::Shower, recob::Cluster>>(ShowerPtr,cluster,"clusterAssociationsbase");
          
      //Associate the hits
      for(auto const& hit: ClusterHits){
	uniqueproducerPtrs.AddSingle<art::Assns<recob::Shower, recob::Hit>>(ShowerPtr, hit,"hitAssociationsbase");
      }
    }

    //Associate the spacepoints
    for(auto const& sp: showerSpacePoints){
      uniqueproducerPtrs.AddSingle<art::Assns<recob::Shower, recob::SpacePoint>>(ShowerPtr,sp,"spShowerAssociationsbase");
    }

    //Loop over the tool data products and add them.
    uniqueproducerPtrs.AddDataProducts(selement_holder);


    //AddAssociations
    int assn_err = 0;
    for(auto const& fShowerTool: fShowerTools){
      assn_err += fShowerTool->AddAssociations(evt,selement_holder);
    }
    if(!fAllowPartialShowers && assn_err > 0){
      mf::LogError("SBNShower") << "A association failed and you are not allowing partial showers. The event will not be added to the event " << std::endl; 
      continue;
    }
   

    //Reset the showerproperty holder.
    selement_holder.ClearAll();
  }
  
  //Put everything in the event.
  uniqueproducerPtrs.MoveAllToEvent(evt);

  //Reset the ptrs to the data products
  uniqueproducerPtrs.reset();

}

DEFINE_ART_MODULE(reco::shower::SBNShower)

