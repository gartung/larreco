//############################################################################
//### Name:        ShowerDirectionCheater                                  ###
//### Author:      Ed Tyley                                                ###
//### Date:        16.07.19                                                ###
//### Description: Cheating tool using truth for shower direction          ###
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
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/RecoAlg/SBNShowerAlg.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larreco/RecoAlg/MCRecoUtils/ShowerUtils.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

//C++ Includes
#include <iostream>
#include <cmath>

//Root Includes
#include "TVector3.h"
#include "TMath.h"
#include "TPrincipal.h"
#include "TVector.h"
#include "TTree.h"

namespace ShowerRecoTools {

class ShowerDirectionCheater:IShowerTool {

  public:

    ShowerDirectionCheater(const fhicl::ParameterSet& pset);

    ~ShowerDirectionCheater();

    //Generic Direction Finder
    int CalculateProperty(const art::Ptr<recob::PFParticle>& pfparticle,
        art::Event& Event,
        reco::shower::ShowerPropertyHolder& ShowerPropHolder
        ) override;

  private:

    // Define standard art tool interface
    void configure(const fhicl::ParameterSet& pset) override;

    TVector3 ShowerPCAVector(std::vector<art::Ptr<recob::SpacePoint> >& spacePoints_pfp, art::FindManyP<recob::Hit>& fmh, TVector3& ShowerCentre);
    double RMSShowerGradient(std::vector<art::Ptr<recob::SpacePoint> >& sps, TVector3& ShowerCentre, TVector3& Direction);
    double CalculateRMS(std::vector<float> perps);

    //Algorithm functions
    shower::SBNShowerAlg fSBNShowerAlg;

    //Services
    detinfo::DetectorProperties const* fDetProp;
    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

    //fcl
    art::InputTag fPFParticleModuleLabel;
    float fNSegments;
    bool fRMSFlip;
    bool fVertexFlip;

    //TTree Branch variables
    TTree* Tree;
    float vertexDotProduct;
    float rmsGradient; 
};


ShowerDirectionCheater::ShowerDirectionCheater(const fhicl::ParameterSet& pset)
  : fSBNShowerAlg(pset.get<fhicl::ParameterSet>("SBNShowerAlg")),
  fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())
{
  configure(pset);
  if (vertexDotProduct||rmsGradient){
    Tree = tfs->make<TTree>("DebugTreeDirCheater", "DebugTree from shower direction cheater");
    if (fVertexFlip) Tree->Branch("vertexDotProduct",&vertexDotProduct);
    if (fRMSFlip)    Tree->Branch("rmsGradient",&rmsGradient);
  }
}

ShowerDirectionCheater::~ShowerDirectionCheater()
{
}

void ShowerDirectionCheater::configure(const fhicl::ParameterSet& pset)
{
  fPFParticleModuleLabel  = pset.get<art::InputTag>("PFParticleModuleLabel","");
  fNSegments              = pset.get<float>        ("NSegments");
  fRMSFlip                = pset.get<bool>         ("RMSFlip");
  fVertexFlip             = pset.get<bool>         ("VertexFlip");
}

int ShowerDirectionCheater::CalculateProperty(const art::Ptr<recob::PFParticle>& pfparticle, art::Event& Event, reco::shower::ShowerPropertyHolder& ShowerPropHolder){

  std::cout <<"#########################################\n"<<
    "hello world direction cheater\n" <<"#########################################\n"<< std::endl;

  const sim::ParticleList& particles = particleInventory->ParticleList();

  // Roll up showers if not already done:
  std::map<int,const simb::MCParticle*> trueParticles;
  std::map<int,std::vector<int> > showersMothers;
  // Make a map of track id to particle
  for (sim::ParticleList::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt){
    const simb::MCParticle *particle = particleIt->second;
    trueParticles[particle->TrackId()] = particle;
    //std::cout<<"Particle ID: "<<particle->TrackId()<<" and PDG: "<<particle->PdgCode()<<std::endl;
  }

  // Loop over daughters and find th`e mothers
  for (sim::ParticleList::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt){
    const simb::MCParticle *particle = particleIt->second;
    const simb::MCParticle *mother   = particle;
    // While the grand mother exists and is an electron or photon
    // Note the true mother will skip this loop and fill itself into the map
    while (mother->Mother()!=0 && trueParticles.find(mother->Mother())!= trueParticles.end() && (TMath::Abs(trueParticles[mother->Mother()]->PdgCode()==11) || TMath::Abs(trueParticles[mother->Mother()]->PdgCode()==22))){
      mother = trueParticles[mother->Mother()];
    }
    //std::cout<<"Mother Done"<<std::endl;
    // This if statement double checks mothers and discards non-shower primary particles
    showersMothers[mother->TrackId()].push_back(particle->TrackId());
  }                   

  //Get the hits from the shower:
  art::Handle<std::vector<recob::PFParticle> > pfpHandle;
  if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
    throw cet::exception("ShowerLinearEnergy") << "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
    return 1;
  }

  //Get the clusters
  art::Handle<std::vector<recob::Cluster> > clusHandle;
  if (!Event.getByLabel(fPFParticleModuleLabel, clusHandle)){
    throw cet::exception("ShowerLinearEnergy") << "Could not get the pandora clusters. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
    return 1;
  }
  art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleModuleLabel);
  std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(pfparticle.key());

  //Get the hit association
  art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleModuleLabel);

  std::vector<art::Ptr<recob::Hit> > showerHits;
  for(auto const& cluster: clusters){

    //Get the hits
    std::vector<art::Ptr<recob::Hit> > hits = fmhc.at(cluster.key());
    showerHits.insert(showerHits.end(),hits.begin(),hits.end());
  }

  //Get the true particle from the shower
  std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain(showersMothers,showerHits,2);

  if(ShowerTrackInfo.first==-99999) {
    std::cout<<"True Shower Not Found"<<std::endl;
    return 1;
  }
  const simb::MCParticle* trueParticle = trueParticles[ShowerTrackInfo.first];
  TVector3 trueDir = {trueParticle->Px(),trueParticle->Py(),trueParticle->Pz()};
  trueDir = trueDir.Unit();

  ShowerPropHolder.SetShowerDirection(trueDir);

  if (fRMSFlip || fVertexFlip){
    //Get the SpacePoints and hits
    art::FindManyP<recob::SpacePoint> fmspp(pfpHandle, Event, fPFParticleModuleLabel);

    if (!fmspp.isValid()){
      throw cet::exception("ShowerStartPosition") << "Trying to get the spacepoint and failed. Something is not configured correctly. Stopping ";
      return 1;
    }

    art::Handle<std::vector<recob::SpacePoint> > spHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, spHandle)){
      throw cet::exception("ShowerStartPosition") << "Could not configure the spacepoint handle. Something is configured incorrectly. Stopping";
      return 1;
    }
    art::FindManyP<recob::Hit> fmh(spHandle, Event, fPFParticleModuleLabel);
    if(!fmh.isValid()){
      throw cet::exception("ShowerStartPosition") << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }
    std::vector<art::Ptr<recob::SpacePoint> > spacePoints = fmspp.at(pfparticle.key());

    if (spacePoints.size()==0) return 1;

    //Get Shower Centre
    float TotalCharge;
    TVector3 ShowerCentre = fSBNShowerAlg.ShowerCentre(spacePoints, fmh, TotalCharge);


    //Check if we are pointing the correct direction or not, First try the start position
    if(ShowerPropHolder.CheckShowerStartPosition() && fVertexFlip){

      //Get the General direction as the vector between the start position and the centre
      TVector3 StartPositionVec = ShowerPropHolder.GetShowerStartPosition();
      TVector3 GeneralDir       = (ShowerCentre - StartPositionVec).Unit();

      //Dot product
      vertexDotProduct = trueDir.Dot(GeneralDir);

      //If the dotproduct is negative the Direction needs Flipping
      if(vertexDotProduct < 0){
        trueDir[0] = - trueDir[0];
        trueDir[1] = - trueDir[1];
        trueDir[2] = - trueDir[2];
      }

      //ShowerPropHolder.SetShowerDirection(trueDir);
    }

    if (fRMSFlip){
      //Otherwise Check against the RMS of the shower. Method adapated from EMShower Thanks Mike.
      rmsGradient = RMSShowerGradient(spacePoints,ShowerCentre,trueDir);
      if(rmsGradient < 0){

        trueDir[0] = - trueDir[0];
        trueDir[1] = - trueDir[1];
        trueDir[2] = - trueDir[2];
      }

      //ShowerPropHolder.SetShowerDirection(trueDir);
    }
    Tree->Fill();
  }    

  std::cout <<"#########################################\n"<<
    "direction cheater Done\n" <<"#########################################\n"<< std::endl;
  return 0;
}


//Function to calculate the RMS at segements of the shower and calculate the gradient of this. If negative then the direction is pointing the opposite way to the correct one
double ShowerDirectionCheater::RMSShowerGradient(std::vector<art::Ptr<recob::SpacePoint> >& sps, TVector3& ShowerCentre, TVector3& Direction){

  //Order the spacepoints
  fSBNShowerAlg.OrderShowerSpacePoints(sps,ShowerCentre,Direction);

  //Get the length of the shower.
  double minProj =fSBNShowerAlg.SpacePointProjection(sps[0],ShowerCentre,Direction);
  double maxProj =fSBNShowerAlg.SpacePointProjection(sps[sps.size()-1],ShowerCentre,Direction);

  double length = (maxProj-minProj);
  double segmentsize = length/fNSegments;

  std::map<int, std::vector<float> > len_segment_map;

  //Split the the spacepoints into segments.
  for(auto const& sp: sps){

    //Get the the projected length
    double len = fSBNShowerAlg.SpacePointProjection(sp,ShowerCentre,Direction);

    //Get the length to the projection
    double  len_perp = fSBNShowerAlg.SpacePointPerpendiular(sp,ShowerCentre,Direction,len);
    
    int sg_len = round(len/segmentsize);
    //TODO: look at this:
    //int sg_len = round(len/segmentsize+fNSegments/2); //Add to make positive
    len_segment_map[sg_len].push_back(len_perp);
  }

  int counter = 0;
  float sumx  = 0;
  float sumy  = 0;
  float sumx2 = 0;
  float sumxy = 0;

  //Get the rms of the segments and caclulate the gradient.
  for(auto const& segment: len_segment_map){
    if (segment.second.size()<2) continue;
      float RMS = CalculateRMS(segment.second);
      std::cout<<"RMS: "<<RMS<<std::endl;
      //Calculate the gradient using regression
      sumx  += segment.first;
      sumy  += RMS;
      sumx2 += segment.first * segment.first;
      sumxy += RMS * segment.first;
      ++counter;
  }
  
    double RMSgradient = (counter*sumxy - sumx*sumy)/(counter*sumx2 - sumx*sumx);
    //std::cout<<"sumX: "<<sumx<<" sumY: "<<sumy<<" sumx2: "<<sumx2<<" sumxy: "<<sumxy<<std::endl;
    std::cout<<"RMS Gradient: "<<RMSgradient<<std::endl;
    return RMSgradient;

  }

  double ShowerDirectionCheater::CalculateRMS(std::vector<float> perps){
    int counter = 0;
    double sum  = 0;
    for (const auto &perp : perps){
      sum= perp*perp;
      ++counter;
    }
    double rms = TMath::Sqrt(sum/(counter-1));

    return rms;
  }

}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerDirectionCheater)
