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
#include "larreco/ShowerFinder/ShowerTools/IShowerDirectionFinder.h"
#include "larreco/ShowerFinder/ShowerTools/IShowerEnergyFinder.h"
#include "larreco/ShowerFinder/ShowerTools/IShowerInitialTrackFinder.h"
#include "larreco/ShowerFinder/ShowerTools/IShowerStartPositionFinder.h"
#include "larreco/ShowerFinder/ShowerTools/IShowerdEdxFinder.h"
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


  //fcl tool names
  std::unique_ptr<ShowerRecoTools::IShowerDirectionFinder>     fShowerDirectionFinder;
  std::unique_ptr<ShowerRecoTools::IShowerEnergyFinder>        fShowerEnergyFinder;
  std::unique_ptr<ShowerRecoTools::IShowerInitialTrackFinder>  fShowerInitialTrackFinder;
  std::unique_ptr<ShowerRecoTools::IShowerStartPositionFinder> fShowerStartPositionFinder;
  std::unique_ptr<ShowerRecoTools::IShowerdEdxFinder>          fShowerdEdxFinder;

  //fcl object names 
  art::InputTag fPFParticleModuleLabel;
  bool          fSecondInteration;

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
  fShowerDirectionFinder     = art::make_tool<ShowerRecoTools::IShowerDirectionFinder>    (pset.get<fhicl::ParameterSet>("ShowerDirectionFinder"));
  fShowerEnergyFinder        = art::make_tool<ShowerRecoTools::IShowerEnergyFinder>       (pset.get<fhicl::ParameterSet>("ShowerEnergyFinder"));
  fShowerInitialTrackFinder  = art::make_tool<ShowerRecoTools::IShowerInitialTrackFinder> (pset.get<fhicl::ParameterSet>("ShowerInitialTrackFinder"));
  fShowerStartPositionFinder = art::make_tool<ShowerRecoTools::IShowerStartPositionFinder>(pset.get<fhicl::ParameterSet>("ShowerStartPositionFinder"));
  fShowerdEdxFinder          = art::make_tool<ShowerRecoTools::IShowerdEdxFinder>         (pset.get<fhicl::ParameterSet>("ShowerdEdxFinder"));

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

  //Give the event and shower property information to the tools.
  fShowerStartPositionFinder->InitialiseEvent(evt,sprop_holder);
  fShowerInitialTrackFinder ->InitialiseEvent(evt,sprop_holder);
  fShowerDirectionFinder    ->InitialiseEvent(evt,sprop_holder);
  fShowerEnergyFinder       ->InitialiseEvent(evt,sprop_holder);
  fShowerdEdxFinder         ->InitialiseEvent(evt,sprop_holder);

  //Loop of the pf particles 
  for(auto const& pfp: pfps){
    
    //Calculate the shower properties 
    TVector3 ShowerStartPosition     = fShowerStartPositionFinder->findShowerStartPosition(pfp);
    sprop_holder.SetShowerStartPosition(ShowerStartPosition);

    recob::Track InitialTrack        = fShowerInitialTrackFinder ->findInitialTrack(pfp);
    sprop_holder.SetInitialTrack       (InitialTrack);

    TVector3 ShowerDirection         = fShowerDirectionFinder    ->findDirection(pfp);
    sprop_holder.SetShowerDirection    (ShowerDirection);
	
    std::vector<double> ShowerEnergy = fShowerEnergyFinder       ->findEnergy(pfp);
    sprop_holder.SetShowerEnergy       (ShowerEnergy);

    std::vector<double> ShowerdEdx   = fShowerdEdxFinder         ->finddEdx(pfp);
    sprop_holder.SetShowerdEdx         (ShowerdEdx);


    //Should we do a second interaction now we have done a first pass of the calculation
    if(fSecondInteration){

      ShowerStartPosition = fShowerStartPositionFinder->findShowerStartPosition(pfp);
      sprop_holder.SetShowerStartPosition(ShowerStartPosition);

      InitialTrack = fShowerInitialTrackFinder->findInitialTrack(pfp);
      sprop_holder.SetInitialTrack(InitialTrack);
      
      ShowerDirection = fShowerDirectionFinder->findDirection(pfp);
      sprop_holder.SetShowerDirection(ShowerDirection);
      
      ShowerEnergy = fShowerEnergyFinder->findEnergy(pfp);
      sprop_holder.SetShowerEnergy(ShowerEnergy);
      
      ShowerdEdx  = fShowerdEdxFinder->finddEdx(pfp);
      sprop_holder.SetShowerdEdx(ShowerdEdx);
    }

    //Get the properties 
    TVector3 ShowerStartPosition_f           = sprop_holder.GetShowerStartPosition();
    const recob::Track InitialTrack_f        = sprop_holder.GetInitialTrack();
    const TVector3 ShowerDirection_f         = sprop_holder.GetShowerDirection();
    const std::vector<double> ShowerEnergy_f = sprop_holder.GetShowerEnergy();
    const std::vector<double> ShowerdEdx_f   = sprop_holder.GetShowerdEdx();
    
    //To Do
    const TVector3 directionErr_f;
    const TVector3 vertexErr_f;
    const std::vector<double> totalEnergyErr_f;
    const std::vector<double> dEdxErr_f;

    //Make the shower 
    recob::Shower shower = recob::Shower(ShowerDirection_f, directionErr_f,ShowerStartPosition_f, vertexErr_f,ShowerEnergy_f,totalEnergyErr_f,ShowerdEdx_f, dEdxErr_f, -999, -999);
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

