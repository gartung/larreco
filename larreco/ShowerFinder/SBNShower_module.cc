//################################################################################
//### Name: SBNShower                                                          ###
//### Author: Dominic Barker                                                   ###
//### Date: 15.05.19                                                           ###
//### Description: Generic Shower Charaterisation module which allows the      ###
//###              the user choose which tool to calculate shower metrics      ###
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
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft includes
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"
#include "larreco/ShowerFinder/ShowerPropertyHolder.h"

//Root Includes
#include "TVector3.h"

//C++ Includes 
#include <vector>

namespace reco {
  namespace shower {
    class SBNShower;
  }
}

class reco::shower::SBNShower: public art::EDProducer {
public:

  SBNShower(fhicl::ParameterSet const& pset);

  void produce(art::Event& evt);
  void reconfigure(fhicl::ParameterSet const& p);

private:

  //fcl object names 
  art::InputTag fPFParticleModuleLabel;
  bool          fSecondInteration;
  bool          fAllowPartialShowers;

  //fcl tools
  std::vector<std::unique_ptr<ShowerRecoTools::IShowerTool> > fShowerTools;
  std::vector<std::string>                                    fShowerToolNames;
};


reco::shower::SBNShower::SBNShower(fhicl::ParameterSet const& pset) :
  EDProducer{pset}
{
  this->reconfigure(pset);
  produces<std::vector<recob::Shower> >();
  produces<std::vector<recob::Track> >();
  produces<art::Assns<recob::Shower, recob::Hit> >();
  produces<art::Assns<recob::Shower, recob::Cluster> >();
  produces<art::Assns<recob::Shower, recob::Track> >();
  produces<art::Assns<recob::Shower, recob::SpacePoint> >();
  produces<art::Assns<recob::Shower, recob::PFParticle> >();
  produces<art::Assns<recob::Track, recob::Hit> >();
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


  //Initialise the other paramters.
  fPFParticleModuleLabel = pset.get<art::InputTag>("PFParticleModuleLabel","pandora");
  fSecondInteration      = pset.get<bool         >("SecondInteration",false);
  fAllowPartialShowers   = pset.get<bool         >("AllowPartialShowers",false);
}
    
void reco::shower::SBNShower::produce(art::Event& evt) {

  // Output -- showers and associations with hits and clusters
  auto showers                = std::make_unique<std::vector<recob::Shower> >();
  auto initialtracks          = std::make_unique<std::vector<recob::Track> >();
  auto clusterAssociations    = std::make_unique<art::Assns<recob::Shower, recob::Cluster> >();
  auto hitShowerAssociations  = std::make_unique<art::Assns<recob::Shower, recob::Hit> >();
  auto trackAssociations      = std::make_unique<art::Assns<recob::Shower, recob::Track> >();
  auto spShowerAssociations   = std::make_unique<art::Assns<recob::Shower, recob::SpacePoint> >();
  auto pfShowerAssociations   = std::make_unique<art::Assns<recob::Shower, recob::PFParticle> >();
  auto hittrackAssociations   = std::make_unique<art::Assns<recob::Track,  recob::Hit> >();

  //Ptr makers for the products 
  const art::PtrMaker<recob::Track> makeTrackPtr(evt);
  const art::PtrMaker<recob::Shower> makeShowerPtr(evt);

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
  reco::shower::ShowerPropertyHolder sprop_holder;

  int i=0;
  //Loop of the pf particles
  for(auto const& pfp: pfps){

    std::cout << "new particle" << std::endl;

    //loop only over showers.
    if(pfp->PdgCode() != 11){continue;}

    //Calculate the shower properties 
    //Loop over the shower tools
    for(auto const& fShowerTool: fShowerTools){

      std::cout << "on next tool:" <<  fShowerToolNames[i]  << std::endl;
      //Calculate the metric
      int err = fShowerTool->CalculateProperty(pfp,evt,sprop_holder);
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
	int err = fShowerTool->CalculateProperty(pfp,evt,sprop_holder);
      
	if(err){
	  mf::LogError("SBNShower") << "Error in shower tool: " << fShowerToolNames[i]  << " with code: " << err << std::endl;
	  if(!fAllowPartialShowers) break;
	}
	++i;
      }
    }

    if(!fAllowPartialShowers){
      bool accept = true; 
      accept = accept & sprop_holder.CheckShowerStartPosition();
      accept = accept & sprop_holder.CheckShowerDirection();
      accept = accept & sprop_holder.CheckShowerEnergy();
      accept = accept & sprop_holder.CheckShowerdEdx();
      accept = accept & sprop_holder.CheckInitialTrack();
      if(!accept){continue;}
    }


    //Get the properties 
    const TVector3                           ShowerStartPosition  = sprop_holder.GetShowerStartPosition();
    const TVector3                           ShowerDirection      = sprop_holder.GetShowerDirection();
    const std::vector<double>                ShowerEnergy         = sprop_holder.GetShowerEnergy();
    const std::vector<double>                ShowerdEdx           = sprop_holder.GetShowerdEdx();
    const recob::Track                       InitialTrack         = sprop_holder.GetInitialTrack();
    const std::vector<art::Ptr<recob::Hit> > InitialTrackHits     = sprop_holder.GetInitialTrackHits(); 
    const int                                BestPlane            = sprop_holder.GetBestPlane();
    
    //To Do
    const TVector3            ShowerDirectionErr              = sprop_holder.GetShowerDirectionErr();
    const TVector3            ShowerStartPositionErrvertexErr = sprop_holder.GetShowerStartPositionErr();
    const std::vector<double> ShowerEnergyErr                 = sprop_holder.GetShowerEnergyErr();
    const std::vector<double> ShowerdEdxErr                   = sprop_holder.GetShowerdEdxErr(); ;
    
    //Check the shower
    std::cout<<"Shower Vertex: X:"<<ShowerStartPosition.X()<<" Y: "<<ShowerStartPosition.Y()<<" Z: "<<ShowerStartPosition.Z()<<std::endl;
    std::cout<<"Shower Direction: X:"<<ShowerDirection.X()<<" Y: "<<ShowerDirection.Y()<<" Z: "<<ShowerDirection.Z()<<std::endl;
    std::cout<<"Shower Track: NumHits: "<<InitialTrack.NumberTrajectoryPoints()<<std::endl;
    std::cout<<"Shower dEdx: size: "<<ShowerdEdx.size()<<" Plane 0: "<<ShowerdEdx.at(0)<<" Plane 1: "<<ShowerdEdx.at(1)<<" Plane 2: "<<ShowerdEdx.at(2)<<std::endl;
    std::cout<<"Shower Energy: size: "<<ShowerEnergy.size()<<" Plane 0: "<<ShowerEnergy.at(0)<<" Plane 1: "<<ShowerEnergy.at(1)<<" Plane 2: "<<ShowerEnergy.at(2)<<std::endl;
    std::cout<<"Shower Best Plane: "<<BestPlane<<std::endl;

    //Make the shower 
    recob::Shower shower = recob::Shower(ShowerDirection, ShowerDirectionErr,ShowerStartPosition, ShowerDirectionErr,ShowerEnergy,ShowerEnergyErr,ShowerdEdx, ShowerdEdxErr, BestPlane, -999);
    showers->push_back(shower);
    art::Ptr<recob::Shower> ShowerPtr(makeShowerPtr(showers->size() - 1));

    //Make the track and the pointer 
    initialtracks->push_back(InitialTrack);
    art::Ptr<recob::Track> InitialTrackPtr(makeTrackPtr(initialtracks->size() - 1));
  
    //Associate the trakc and the shower 
    trackAssociations->addSingle(ShowerPtr,InitialTrackPtr);

    //Associate the pfparticle 
    pfShowerAssociations->addSingle(ShowerPtr,pfp);

    //Associate the track and the initial hits.
    for(auto const& hit: InitialTrackHits){
      hittrackAssociations->addSingle(InitialTrackPtr,hit);
    }

    //Get the associated hits,clusters and spacepoints
    std::vector<art::Ptr<recob::Cluster> >    showerClusters    = fmcp.at(pfp.key());
    std::vector<art::Ptr<recob::SpacePoint> > showerSpacePoints = fmspp.at(pfp.key());

    //Add the hits for each "cluster"
    for(auto const& cluster: showerClusters){

      //Associate the clusters 
      std::vector<art::Ptr<recob::Hit> > ClusterHits = fmh.at(cluster.key());
      clusterAssociations->addSingle(ShowerPtr,cluster);

      //Associate the hits
      for(auto const& hit: ClusterHits){
	 hitShowerAssociations->addSingle(ShowerPtr, hit);
      }
    }

    //Associate the spacepoints
    for(auto const& sp: showerSpacePoints){
      spShowerAssociations->addSingle(ShowerPtr,sp);
    }

    //Reset the showerproperty holder.
    sprop_holder.clear();
  }
  
  // Put in event
  evt.put(std::move(showers));
  evt.put(std::move(initialtracks));
  evt.put(std::move(hitShowerAssociations));
  evt.put(std::move(clusterAssociations));
  evt.put(std::move(trackAssociations));
  evt.put(std::move(spShowerAssociations));
  evt.put(std::move(hittrackAssociations));
  evt.put(std::move( pfShowerAssociations));
}

DEFINE_ART_MODULE(reco::shower::SBNShower)

