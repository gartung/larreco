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
  //  produces<std::vector<recob::SpacePoint> >();
  produces<art::Assns<recob::Shower, recob::Hit> >();
  produces<art::Assns<recob::Shower, recob::Cluster> >();
  produces<art::Assns<recob::Track, recob::Shower> >();
  produces<art::Assns<recob::Shower, recob::SpacePoint> >();
  produces<art::Assns<recob::SpacePoint, recob::Hit> >();
}

void reco::shower::SBNShower::reconfigure(fhicl::ParameterSet const& pset) {

  //Intialise the tools 
  const fhicl::ParameterSet& ShowerTools = pset.get<fhicl::ParameterSet>("ShowerFinderTools");
  for(const std::string& ShowerTool : ShowerTools.get_pset_names()){
    const fhicl::ParameterSet& ShowerToolParamSet = ShowerTools.get<fhicl::ParameterSet>(ShowerTool);
    fShowerTools.push_back(art::make_tool<ShowerRecoTools::IShowerTool>(ShowerToolParamSet));
    fShowerToolNames.push_back(ShowerToolParamSet.id().to_string());
  }

  //Initialise the other paramters.
  fPFParticleModuleLabel = pset.get<art::InputTag>("PFParticleModuleLabel","pandora");
  fSecondInteration      = pset.get<bool         >("SecondInteration",false);
}
    
void reco::shower::SBNShower::produce(art::Event& evt) {

  // Output -- showers and associations with hits and clusters
  std::unique_ptr<std::vector<recob::Shower> > showers(new std::vector<recob::Shower>);
  std::unique_ptr<std::vector<recob::Track> > initialtracks(new std::vector<recob::Track>);
  std::unique_ptr<art::Assns<recob::Shower, recob::Cluster> > clusterAssociations(new art::Assns<recob::Shower, recob::Cluster>);
  std::unique_ptr<art::Assns<recob::Shower, recob::Hit> > hitShowerAssociations(new art::Assns<recob::Shower, recob::Hit>);
  std::unique_ptr<art::Assns<recob::Track, recob::Shower> > trackAssociations(new art::Assns<recob::Track, recob::Shower>);
  std::unique_ptr<art::Assns<recob::Shower, recob::SpacePoint> > spShowerAssociations(new art::Assns<recob::Shower, recob::SpacePoint>);
  std::unique_ptr<art::Assns<recob::SpacePoint, recob::Hit> > hitSpAssociations(new art::Assns<recob::SpacePoint, recob::Hit>);

  //Ptr for initial track assans
  const art::PtrMaker<recob::Track> makeTrackPtr(evt);

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

  //Holder to pass to the functions, contains the 5 properties of the shower 
  // - Start Poistion
  // - Direction
  // - Initial Track
  // - Energy 
  // - dEdx 
  reco::shower::ShowerPropertyHolder sprop_holder;

  int i=0;
  //Loop of the pf particles
  for(auto const& pfp: pfps){
    
    //Calculate the shower properties 
    //Loop over the shower tools
    for(auto const& fShowerTool: fShowerTools){

      std::cout << "on next tool" << std::endl;
      //Calculate the metric
      int err = fShowerTool->findMetric(pfp,evt,sprop_holder);
      if(err){
	throw cet::exception("SBNShower") << "Error in shower tool: " << fShowerToolNames[i] << " with code: " << err << std::endl;
	break;
      }
      ++i
    }

    //Should we do a second interaction now we have done a first pass of the calculation
    if(fSecondInteration){
      for(auto const& fShowerTool: fShowerTools){
	//Calculate the metric
	int err = fShowerTool->findMetric(pfp,evt,sprop_holder);
      }
      if(err){
	throw cet::exception("SBNShower") << "Error in shower tool: " << fShowerTool->GetToolName() << " with code: " << err << std::endl;
	break;
      }
    }

    //Get the properties 
    const TVector3            ShowerStartPosition  = sprop_holder.GetShowerStartPosition();
    const TVector3            ShowerDirection      = sprop_holder.GetShowerDirection();
    const std::vector<double> ShowerEnergy         = sprop_holder.GetShowerEnergy();
    const std::vector<double> ShowerdEdx           = sprop_holder.GetShowerdEdx();
    const recob::Track        InitialTrack         = sprop_holder.GetInitialTrack();
    
    //To Do
    const TVector3            ShowerDirectionErr              = sprop_holder.GetShowerDirectionErr();
    const TVector3            ShowerStartPositionErrvertexErr = sprop_holder.GetShowerStartPositionErr();
    const std::vector<double> ShowerEnergyErr                 = sprop_holder.GetShowerEnergyErr();
    const std::vector<double> ShowerdEdxErr                   = sprop_holder.GetShowerdEdxErr(); ;

    //Make the shower 
    recob::Shower shower = recob::Shower(ShowerDirection, ShowerDirectionErr,ShowerStartPosition, ShowerDirectionErr,ShowerEnergy,ShowerEnergyErr,ShowerdEdx, ShowerdEdxErr, -999, -999);
    showers->push_back(shower);
    showers->back().set_id(showers->size()-1);
    
    //Make the track and the pointer 
    initialtracks->push_back(sprop_holder.GetInitialTrack());
    art::Ptr<recob::Track> initialtrack(makeTrackPtr(initialtracks->size() - 1));
    
    //Get the associated hits,clusters and spacepoints
    std::vector<art::Ptr<recob::Hit> >        showerHits;
    std::vector<art::Ptr<recob::Cluster> >    showerClusters    = fmcp.at(pfp.key());
    std::vector<art::Ptr<recob::SpacePoint> > showerSpacePoints = fmspp.at(pfp.key());

    //Add the hits for each "cluster"
    for (size_t iclu = 0; iclu<showerClusters.size(); ++iclu){
      std::vector<art::Ptr<recob::Hit> > ClusterHits = fmh.at(showerClusters[iclu].key());
      showerHits.insert(showerHits.end(), ClusterHits.begin(), ClusterHits.end());
    }
    
    //Create the associates
    util::CreateAssn(*this, evt, *(showers.get()), showerHits,        *(hitShowerAssociations.get()));
    util::CreateAssn(*this, evt, *(showers.get()), showerClusters,    *(clusterAssociations.get()));
    util::CreateAssn(*this, evt, *(showers.get()), initialtrack,      *(trackAssociations.get()));
    util::CreateAssn(*this, evt, *(showers.get()), showerSpacePoints, *(spShowerAssociations.get()));

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
  evt.put(std::move(hitSpAssociations));

}

DEFINE_ART_MODULE(reco::shower::SBNShower)

