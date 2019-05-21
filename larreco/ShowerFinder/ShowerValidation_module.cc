////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       ShowerValidation                                                                          //
// Modlue Type: Analyser                                                                                  //
// File         ShowerValidation_module.cc                                                                //
// Author:      Dominic Barker dominic.barker@sheffield.ac.uk                                             //
//              Edward Tyley e.tyley@sheffield.ac.uk                                                      //
//                                                                                                        //
// Usage:       The Validation module takes an array of shower modules and perfoms truth matching on      //
//              reconstructed EM showers. It calculates the truth information from geant information      //
//              such as the diretion and starting position and energy deposited. Metrics are then defined //
//              to compare the efficiency of the shower reconstruction modules. These are presented       //
//              as histograms and in terms of energy. Energy plots look at the metrics for specific       //
//              energies between +- fEnergyWidth of the values described in the fcl table                 //
//                                                                                                        //
// Updates:     31.10.2018  Clustering Validation Added                                                   //
//              13.11.2018  Hit Validation, 2D Histograms and refactoring of the code                     //
//              22.01.2019  TTree output added                                                            //
//              18.03.2019  Added PFParticle Validation                                                   //
//                                                                                                        //
////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//Larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/MCRecoUtils/RecoUtils.h"
#include "larreco/RecoAlg/MCRecoUtils/ShowerUtils.h"

//Root Includes
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

//C++ Includes
#include <vector>
#include <iostream>


namespace ana {
  class ShowerValidation;
}

class ana::ShowerValidation : public art::EDAnalyzer {
public:

  ShowerValidation(const fhicl::ParameterSet& pset);

  void analyze(const art::Event& evt);
  void endJob();
  void beginJob();

  void initTree(TTree* Tree, std::string branchName, std::map<std::string,std::vector<float> >& Metric, std::vector<std::string> fShowerModuleLabels);
  void initClusterTree(TTree* Tree,  std::string branchName,  std::map<std::string,std::vector<std::vector<std::vector<float> > > >& Metric, std::vector<std::string> fShowerModuleLabels);

  void ClusterValidation(std::vector<art::Ptr<recob::Cluster> >& clusters,
			 const art::Event& evt,
			 art::Handle<std::vector<recob::Cluster> >& clusterHandle,
			 std::map<int,std::vector<int> >& ShowerMotherTrackIDs,
			 std::map<int,float>& MCTrack_Energy_map,
			 std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > > & MCTrack_hit_map,
			 int& TrueShowerID,
			 float& simenergy,
			 std::string & fShowerModuleLabel);

  void PFPValidation(std::vector<art::Ptr<recob::Cluster> >& clusters,
		     art::Ptr<recob::PFParticle>  pfp,
		     const art::Event& evt,
		     art::Handle<std::vector<recob::Cluster> >& clusterHandle,
		     std::map<int,std::vector<int> >& ShowerMotherTrackIDs,
		     std::map<int,float>& MCTrack_Energy_map,
		     std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > > & MCTrack_hit_map,
		     std::string & fShowerModuleLabel);


private:

  //fcl parameters
  std::string fGenieGenModuleLabel;
  std::string fLArGeantModuleLabel;
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  std::string fPFParticleLabel;

  bool  fUseBiggestShower;
  bool  fDrawCanvases;
  bool  fFillOnlyClosestShower;
  bool  fRemoveNonContainedParticles;
  bool  fPFPValidation;
  int   fVerbose;
  int   fMinHitSize;
  float fEnergyWidth;
  float fSimEnergyCut;
  float fDensityCut;
  float fMaxSimEnergy;

  std::vector<float>       fEnergies;
  std::vector<std::string> fShowerModuleLabels;
  std::vector<std::string> fHitModuleLabels;

  std::map<std::string,std::vector<float> > ShowerDirection_X_TreeVal;
  std::map<std::string,std::vector<float> > ShowerDirection_Y_TreeVal;
  std::map<std::string,std::vector<float> > ShowerDirection_Z_TreeVal;
  std::map<std::string,std::vector<float> > ShowerStart_X_TreeVal;
  std::map<std::string,std::vector<float> > ShowerStart_Y_TreeVal;
  std::map<std::string,std::vector<float> > ShowerStart_Z_TreeVal;
  std::map<std::string,std::vector<float> > ShowerLength_TreeVal;
  std::map<std::string,std::vector<float> > ShowerEnergyDiff_TreeVal;
  std::map<std::string,std::vector<float> > ShowerdEdx_TreeVal;
  std::map<std::string,std::vector<float> > EventSeggy_TreeVal;
  std::map<std::string,std::vector<float> > ShowerEnergyCompleteness_TreeVal;
  std::map<std::string,std::vector<float> > ShowerHitsPurity_TreeVal;
  std::map<std::string,std::vector<float> > ShowerHitsCompleteness_TreeVal;
  std::map<std::string,std::vector<float> > ShowerEnergyPurity_TreeVal;
  std::map<std::string,std::vector<float> > ShowerEnergy_TreeVal;
  std::map<std::string,std::vector<float> > ShowerHitNum_TreeVal;
  std::map<std::string,std::vector<float> > ShowerTotalEnergyDiff_TreeVal;
  std::map<std::string,std::vector<float> > ShowerMag_TreeVal;
  std::map<std::string,std::vector<float> > ShowerDirectionDiff_TreeVal;
  std::map<std::string,std::vector<float> > ShowerRecoEnergyVsTrueEnergyinRecoShower_TreeVal;
  std::map<std::string,std::vector<float> > ShowerTrueEnergy_TreeVal;
  std::map<std::string,std::vector<float> > TrueEnergy_TreeVal;
  std::map<std::string,std::vector<float> > TrueHitNum_TreeVal;
  std::map<std::string,std::vector<float> > ShowerBestPlane_TreeVal;
  std::map<std::string,std::vector<float> > GeoProjectionMatched_TreeVal;
  std::map<std::string,std::vector<float> > PandoraTrackNum_TreeVal;

  std::map<std::string,std::vector<float> > PFP_Num_Neutrinos_TreeVal;
  std::map<std::string,std::vector<float> > PFP_Vertex_Distance_X_TreeVal;
  std::map<std::string,std::vector<float> > PFP_Vertex_Distance_Y_TreeVal;
  std::map<std::string,std::vector<float> > PFP_Vertex_Distance_Z_TreeVal;
  std::map<std::string,std::vector<float> > PFP_Vertex_Distance_Mag_TreeVal;
  std::map<std::string,std::vector<float> > PFP_Tracks_TreeVal;
  std::map<std::string,std::vector<float> > PFP_Showers_TreeVal;
  std::map<std::string,std::vector<float> > PFP_Showers_Vertices_TreeVal;
  std::map<std::string,std::vector<float> > PFP_ProjectionMatched_TreeVal;
  std::map<std::string,std::vector<float> > PFP_HitsCompleteness_TreeVal;
  std::map<std::string,std::vector<float> > PFP_EnergyCompleteness_TreeVal;
  std::map<std::string,std::vector<float> > PFP_HitsPurity_TreeVal;
  std::map<std::string,std::vector<float> > PFP_EnergyPurity_TreeVal;
  std::map<std::string,std::vector<float> > PFP_Track_ProjectionMatched_TreeVal;
  std::map<std::string,std::vector<float> > PFP_Track_HitsCompleteness_TreeVal;
  std::map<std::string,std::vector<float> > PFP_Track_EnergyCompleteness_TreeVal;
  std::map<std::string,std::vector<float> > PFP_Track_HitsPurity_TreeVal;
  std::map<std::string,std::vector<float> > PFP_Track_EnergyPurity_TreeVal;
  std::map<std::string,std::vector<float> > PFP_Shower_ProjectionMatched_TreeVal;
  std::map<std::string,std::vector<float> > PFP_Shower_HitsCompleteness_TreeVal;
  std::map<std::string,std::vector<float> > PFP_Shower_EnergyCompleteness_TreeVal;
  std::map<std::string,std::vector<float> > PFP_Shower_HitsPurity_TreeVal;
  std::map<std::string,std::vector<float> > PFP_Shower_EnergyPurity_TreeVal;

  std::map<std::string,std::vector<float> > Event_NumTrueShowers_TreeVal;
  std::map<std::string,std::vector<float> > Event_NumTrueShowers_EnergyCut_TreeVal;
  std::map<std::string,std::vector<float> > Event_NumTrueShowers_DensityCut_TreeVal;
  std::map<std::string,std::vector<float> > Event_TrueShowerE_TreeVal;
  std::map<std::string,std::vector<float> > Event_TrueShowerE_EnergyCut_TreeVal;
  std::map<std::string,std::vector<float> > Event_TrueShowerE_DensityCut_TreeVal;
  std::map<std::string,std::vector<float> > PosDir_TreeVal;

  std::map<std::string,std::vector<std::vector<std::vector<float> > > >  ClusterProjectionMatchedEnergy_TreeVal;
  std::map<std::string,std::vector<std::vector<std::vector<float> > > >  ClusterCompletenessEnergy_TreeVal;
  std::map<std::string,std::vector<std::vector<std::vector<float> > > >  ClusterPurityEnergy_TreeVal;
  std::map<std::string,std::vector<std::vector<std::vector<float> > > >  ClusterCompletenessHits_TreeVal;
  std::map<std::string,std::vector<std::vector<std::vector<float> > > >  ClusterPurityHits_TreeVal;
  std::map<std::string,std::vector<std::vector<std::vector<float> > > >  ClusterCompPurityEnergy_TreeVal;
  std::map<std::string,std::vector<std::vector<std::vector<float> > > >  ClusterCompPurityHits_TreeVal;

  std::map<std::string,std::vector<std::vector<float> > > HitCompletenessEnergy_TreeVal;
  std::map<std::string,std::vector<std::string> > ShowerStartEndProcess_TreeVal;


  float EventRun_TreeVal;
  float EventSubrun_TreeVal;
  float EventNumber_TreeVal;

  //TTree
  TTree* Tree;

  //Service handles
  art::ServiceHandle<cheat::BackTrackerService> backtracker;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<art::TFileService> tfs;

  int numevents;
  int containmentCutNum;
  int numshowers;
  int numshowerspassTPC;
  int numshowerspassdensity;
  int numshoowerspassenergy;
  int numrecoshowers;
  int numrecoshowersana;

};

ana::ShowerValidation::ShowerValidation(const fhicl::ParameterSet& pset) : EDAnalyzer(pset){

  fGenieGenModuleLabel         = pset.get<std::string>("GenieGenModuleLabel");
  fLArGeantModuleLabel         = pset.get<std::string>("LArGeantModuleLabel");
  fHitsModuleLabel             = pset.get<std::string>("HitsModuleLabel");
  fTrackModuleLabel            = pset.get<std::string>("TrackModuleLabel");
  fPFParticleLabel             = pset.get<std::string>("PFParticleLabel");
  fShowerModuleLabels          = pset.get<std::vector<std::string> >("ShowerModuleLabels");
  fHitModuleLabels             = pset.get<std::vector<std::string> >("HitModuleLabels");
  fEnergies                    = pset.get<std::vector<float> >("Energies");
  fUseBiggestShower            = pset.get<bool>("UseBiggestShower");
  fDrawCanvases                = pset.get<bool>("DrawCanvases");
  fFillOnlyClosestShower       = pset.get<bool>("FillOnlyClosestShower");
  fRemoveNonContainedParticles = pset.get<bool>("RemoveNonContainedParticles");
  fPFPValidation               = pset.get<bool>("PFPValidation");
  fVerbose                     = pset.get<int>("Verbose");
  fMinHitSize                  = pset.get<int>("MinHitSize");
  fEnergyWidth                 = pset.get<float>("EnergyWidth");
  fSimEnergyCut                = pset.get<float>("SimEnergyCut");
  fDensityCut                  = pset.get<float>("DensityCut");
  fMaxSimEnergy                = pset.get<float>("MaxSimEnergy");
}

void ana::ShowerValidation::initTree(TTree* Tree, std::string branchName, std::map<std::string,std::vector<float> >& Metric,   std::vector<std::string> fShowerModuleLabels){
  for(unsigned int j=0; j<fShowerModuleLabels.size(); ++j){
    std::string branchString = branchName + "_" + fShowerModuleLabels[j];
    const char* branchChar   = branchString.c_str();
    Tree->Branch(branchChar,"std::vector<float>", &Metric[fShowerModuleLabels[j]], 32000, 0);
  }
}

void ana::ShowerValidation::initClusterTree(TTree* Tree, std::string branchName,  std::map<std::string,std::vector<std::vector<std::vector<float> > > >& Metric, std::vector<std::string> fShowerModuleLabels){
  for(unsigned int j=0; j<fShowerModuleLabels.size(); ++j){
    std::string branchString = branchName + "_" + fShowerModuleLabels[j];
    const char* branchChar   = branchString.c_str();
    Tree->Branch(branchChar,"std::map<std::string,std::vector<std::vector<std::vector<float> > > >"
		 , &Metric[fShowerModuleLabels[j]], 32000, 0);
  }
}



void ana::ShowerValidation::beginJob() {


  numevents = 0;
  containmentCutNum = 0;
  numshowers = 0;
  numshowerspassTPC = 0;
  numshowerspassdensity = 0;
  numshoowerspassenergy = 0;
  numrecoshowers = 0;
  numrecoshowersana = 0;

  Tree = tfs->make<TTree>("MetricTree", "Tree Holding all metric information");
  gInterpreter->GenerateDictionary("vector<vector<vector<float> > >","vector");
  Tree->Branch("EventRun", &EventRun_TreeVal, 32000, 0);
  Tree->Branch("EventSubrun", &EventSubrun_TreeVal, 32000, 0);
  Tree->Branch("EventNumber", &EventNumber_TreeVal, 32000, 0);



  initTree(Tree,"ShowerDirection_X",ShowerDirection_X_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerDirection_Y",ShowerDirection_Y_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerDirection_Z",ShowerDirection_Z_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerStart_X",ShowerStart_X_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerStart_Y",ShowerStart_Y_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerStart_Z",ShowerStart_Z_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerLength",ShowerLength_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerEnergyDiff",ShowerEnergyDiff_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerdEdx",ShowerdEdx_TreeVal,fShowerModuleLabels);
  initTree(Tree,"EventSeggy",EventSeggy_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerEnergyCompleteness",ShowerEnergyCompleteness_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerHitsPurity",ShowerHitsPurity_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerHitsCompleteness",ShowerHitsCompleteness_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerEnergyPurity",ShowerEnergyPurity_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerEnergy",ShowerEnergy_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerHitNum",ShowerHitNum_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerTotalEnergyDiff",ShowerTotalEnergyDiff_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerMag",ShowerMag_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerDirectionDiff",ShowerDirectionDiff_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerRecoEnergyVsTrueEnergyinRecoShower",ShowerRecoEnergyVsTrueEnergyinRecoShower_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerTrueEnergy",ShowerTrueEnergy_TreeVal,fShowerModuleLabels);
  initTree(Tree,"TrueEnergy",TrueEnergy_TreeVal,fShowerModuleLabels);
  initTree(Tree,"TrueHitNum",TrueHitNum_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerBestPlane",ShowerBestPlane_TreeVal,fShowerModuleLabels);
  initTree(Tree,"GeoProjectionMatched",GeoProjectionMatched_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PandoraTrackNum",PandoraTrackNum_TreeVal,fShowerModuleLabels);

  initTree(Tree,"PFP_Num_Neutrinos",PFP_Num_Neutrinos_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Vertex_Distance_X",PFP_Vertex_Distance_X_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Vertex_Distance_Y",PFP_Vertex_Distance_Y_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Vertex_Distance_Z",PFP_Vertex_Distance_Z_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Vertex_Distance_Mag",PFP_Vertex_Distance_Mag_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Tracks",PFP_Tracks_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Showers",PFP_Showers_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Showers_Vertces",PFP_Showers_Vertices_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_ProjectionMatched",PFP_ProjectionMatched_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_HitsCompleteness",PFP_HitsCompleteness_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_EnergyCompleteness",PFP_EnergyCompleteness_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_HitsPurity",PFP_HitsPurity_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_EnergyPurity",PFP_EnergyPurity_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Track_ProjectionMatched",PFP_Track_ProjectionMatched_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Track_HitsCompleteness",PFP_Track_HitsCompleteness_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Track_EnergyCompleteness",PFP_Track_EnergyCompleteness_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Track_HitsPurity",PFP_Track_HitsPurity_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Track_EnergyPurity",PFP_Track_EnergyPurity_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Shower_ProjectionMatched",PFP_Shower_ProjectionMatched_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Shower_HitsCompleteness",PFP_Shower_HitsCompleteness_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Shower_EnergyCompleteness",PFP_Shower_EnergyCompleteness_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Shower_HitsPurity",PFP_Shower_HitsPurity_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Shower_EnergyPurity",PFP_Shower_EnergyPurity_TreeVal,fShowerModuleLabels);
  initTree(Tree,"Event_NumTrueShowers",Event_NumTrueShowers_TreeVal,fShowerModuleLabels);
  initTree(Tree,"Event_NumTrueShowers_EnergyCut",Event_NumTrueShowers_EnergyCut_TreeVal,fShowerModuleLabels);
  initTree(Tree,"Event_NumTrueShowers_DensityCut",Event_NumTrueShowers_DensityCut_TreeVal,fShowerModuleLabels);
  initTree(Tree,"Event_TrueShowerE",Event_TrueShowerE_TreeVal,fShowerModuleLabels);
  initTree(Tree,"Event_TrueShowerE_EnergyCut",Event_TrueShowerE_EnergyCut_TreeVal,fShowerModuleLabels);
  initTree(Tree,"Event_TrueShowerE_DensityCut",Event_TrueShowerE_DensityCut_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PosDir",PosDir_TreeVal,fShowerModuleLabels);




  initTree(Tree,"ShowerdEdx",ShowerdEdx_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerMag",ShowerMag_TreeVal,fShowerModuleLabels);
  initTree(Tree,"ShowerDirectionDiff",ShowerDirectionDiff_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Tracks",PFP_Tracks_TreeVal,fShowerModuleLabels);
  initTree(Tree,"PFP_Showers",PFP_Showers_TreeVal,fShowerModuleLabels);

  for(unsigned int j=0; j<fHitModuleLabels.size(); ++j){
    std::cout << geom->Nplanes() << std::endl;
    HitCompletenessEnergy_TreeVal[fHitModuleLabels[j]].resize(geom->Nplanes());

    std::string hitString = "HitCompletenessEnergy_" + fShowerModuleLabels[j];
    const char* hitChar   = hitString.c_str();

    std::string processString = "ShowerStartEndProcess_" + fShowerModuleLabels[j];
    const char* processChar   = processString.c_str();

    Tree->Branch(hitChar,"std::vector<std::vector<float> > >", &HitCompletenessEnergy_TreeVal[fHitModuleLabels[j]], 32000, 0);

    Tree->Branch(processChar,"<std::string,std::vector<std::string>", &ShowerStartEndProcess_TreeVal[fShowerModuleLabels[j]], 32000, 0);

  }
}


void ana::ShowerValidation::analyze(const art::Event& evt) {

  EventRun_TreeVal = evt.run();
  EventSubrun_TreeVal = evt.subRun();
  EventNumber_TreeVal = evt.event();

  ++numevents;

  //Getting  MC truth information
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if(evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    {art::fill_ptr_vector(mclist, mctruthListHandle);}

  //Getting the SimWire Information
  //Get the SimChannels so that we can find the IDEs deposited on them.
  art::Handle<std::vector<sim::SimChannel> > simChannelHandle;
  std::vector<art::Ptr<sim::SimChannel> > simchannels;
  if(evt.getByLabel(fLArGeantModuleLabel,simChannelHandle))
    {art::fill_ptr_vector(simchannels, simChannelHandle);}

  //Getting the Hit Information
  art::Handle<std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hits;
  if(evt.getByLabel(fHitsModuleLabel,hitListHandle))
    {art::fill_ptr_vector(hits, hitListHandle);}

  //Get the track Information (hopfully you have pandora track)
  art::Handle<std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracks;
  if(evt.getByLabel(fTrackModuleLabel,trackListHandle))
    {art::fill_ptr_vector(tracks, trackListHandle);}

  //I think that doing getManyByType kind of initalises the handles giving every particle product id. Doing this allows us to find handles for the individal hits later.

  //Get all the hits
  std::vector<art::Handle<std::vector<recob::Hit> > > hitHandles;
  evt.getManyByType(hitHandles);

  //Get all the clusters
  std::vector<art::Handle<std::vector<recob::Cluster> > > clusterHandles;
  evt.getManyByType(clusterHandles);

  //Get all the pfparticles
  std::vector<art::Handle<std::vector<recob::PFParticle> > > pfpHandles;
  evt.getManyByType(pfpHandles);


  //###############################################
  //### Get the Truth information for the event ###
  //###############################################

  //List the particles in the event
  const sim::ParticleList& particles = particleInventory->ParticleList();

  //Loop over the particles
  std::map<int,const simb::MCParticle*> trueParticles;
  std::map<int,const simb::MCParticle*> trueInitialParticles;
  std::map<int,bool> mcparticlescontained;
  std::map<int,float> trueParticleEnergy;
  std::map<int,std::vector<int> > ShowersMothers; //Mothers are the key Daughters are in the vector.

  int num_of_showers =0;
  int num_of_showers_viaEcut = 0;
  int num_of_showers_viaDensitycut = 0;
  std::vector<int> E_of_showers;
  std::vector<int> E_of_showers_viaEcut;
  std::vector<int> E_of_showers_viaDensitycut;
  float simenergy=-99999999;

  std::map<int,std::map<geo::PlaneID,int> > MCWires_Track_Map = RecoUtils::NumberofMCWiresHitMap(simchannels);

  //Make a map of Track id and pdgcode
  for (sim::ParticleList::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt){
    const simb::MCParticle *particle = particleIt->second;
    trueParticleEnergy[particle->TrackId()] = 0;
    trueParticles[particle->TrackId()] = particle;
  }

  for (sim::ParticleList::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt){
    const simb::MCParticle *particle = particleIt->second;

    bool contained = true;

    //Particles with mother 0 are the initial particles (neutrino events this is the particles generated after the interaction. Keep note of these.
    simenergy = particle->E();
    trueInitialParticles[particle->TrackId()] = particle;

    //Check to see if the particle is contained and find the trajectory of the particle
    //Get the number of Traj points to loop over
    unsigned int TrajPoints = particle->NumberTrajectoryPoints();

    //Get the startpoistion so we can get the initial tpc.
    const TLorentzVector StartPositionTrajP = particle->Position(0);
    double start_vtx[3] = {StartPositionTrajP.X() ,StartPositionTrajP.Y(), StartPositionTrajP.Z()};
    geo::TPCID init_idtpc = geom->FindTPCAtPosition(start_vtx);

    //Loop over the trajectory points (they are in order). Loop to find the start point.
    for(unsigned int TrajPoints_it=0; TrajPoints_it<TrajPoints; ++TrajPoints_it){

      //Find the vertex of the vector
      const TLorentzVector PositionTrajP = particle->Position(TrajPoints_it);
      double vtx[3] = {PositionTrajP.X() ,PositionTrajP.Y(), PositionTrajP.Z()};

      //Find if the vertex is in the TPC. If so make it the start point.
      geo::TPCID idtpc = geom->FindTPCAtPosition(vtx);

      //Remove because this is crap and we need to think of something better.
      if(idtpc != init_idtpc){std::cout <<"Particle outside the TPC" << std::endl; contained=false; break;}
    }

    if(contained){mcparticlescontained[particle->TrackId()] = true;}
    else{mcparticlescontained[particle->TrackId()] = false;}

    if(fVerbose > 1){std::cout << "True Particle with track ID: " << particle->TrackId() << " Has code of: " << particle->PdgCode() << " and Energy of: " << particle->E() << " With Mother: " << particle->Mother() << " Proccess: " << particle->Process() << " End Process: "  << particle->EndProcess() << std::endl;}


    //Find The electron or photon mother of the shower.
    const simb::MCParticle * particle_temp = particle;
    int ShowerMotherID = particle_temp->TrackId();

    if(particle_temp->Mother() != 0){
      if(trueParticles.find(ShowerMotherID) == trueParticles.end() || trueParticles.find(particle_temp->Mother()) == trueParticles.end()){continue;}


      while(particle_temp->Mother() != 0 && (TMath::Abs(trueParticles[particle_temp->Mother()]->PdgCode()) == 11 || trueParticles[particle_temp->Mother()]->PdgCode() == 22)){

	ShowerMotherID = particle_temp->Mother();

	particle_temp =  trueParticles[particle_temp->Mother()];

	if(trueParticles.find(particle_temp->Mother()) == trueParticles.end()){break;}
      }

    }

    if(ShowersMothers.find(ShowerMotherID) == ShowersMothers.end() && (TMath::Abs(trueParticles[ShowerMotherID]->PdgCode()) == 11 || trueParticles[ShowerMotherID]->PdgCode() == 22)){ShowersMothers[ShowerMotherID].push_back(ShowerMotherID);}
  }

  std::map<int,int> Daughters_used;

  //Find the all the daughter particles associated with the mothers
  for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end(); ++showermother){

    ++numshowers;

    int Gen_Num = 0;
    int Daughter_id = showermother->first;
    while(Gen_Num!=0 || Daughters_used[showermother->first] != trueParticles[showermother->first]->NumberDaughters()){

      if(Daughters_used[Daughter_id] == trueParticles[Daughter_id]->NumberDaughters() || (trueParticles[Daughter_id]->PdgCode() != 22 && TMath::Abs(trueParticles[Daughter_id]->PdgCode()) != 11 )){
	if(trueParticles.find(Daughter_id) == trueParticles.end()){continue;}
	(showermother->second).push_back(Daughter_id);
	Daughter_id = trueParticles[Daughter_id]->Mother();
	--Gen_Num;
	++Daughters_used[Daughter_id];
	continue;
      }

      //Sometimes the duaghter is not in the particle list and so breaks the code. I presume its because its too small in energy to propergate by the simulation.
      if(trueParticles.find(trueParticles[Daughter_id]->Daughter(Daughters_used[Daughter_id])) != trueParticles.end()){
	Daughter_id = trueParticles[Daughter_id]->Daughter(Daughters_used[Daughter_id]);
	++Gen_Num;
      }
      else {
	++Daughters_used[Daughter_id];
      }
    }
  }


  if(fVerbose > 1){
    for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end(); ++showermother){
      const simb::MCParticle *motherparticle = trueParticles[showermother->first]; //test
      std::cout << " A Mother has track id: " << showermother->first << " with Energy: "<< motherparticle->E() << std::endl;
      for(std::vector<int>::iterator daughterID=(showermother->second).begin(); daughterID!=(showermother->second).end();++daughterID){
	std::cout << " Has a daughter of id: " << *daughterID << std::endl;
      }
    }
  }

  //Time to cut the true showers and make sure they are a shower.
  for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end();){

    //I've read that pair production starts to dominate at around ~100 MeV so to find how many showers we expect loop over the mother particle. Pi0=143.97 MeV min gammas = 71.985 MeV which is greater than that from electrons at ~100MeV so pi0 should always shower? So cut on anything below 100MeV in energy.

    //It ain't a shower I'm interested in if it didn't start with a pi0 or electron...probably.
    const simb::MCParticle *motherparticle = trueParticles[showermother->first];
    int pdgcode = motherparticle->PdgCode();
    if(TMath::Abs(pdgcode) == 11 || pdgcode == 22){
      ++num_of_showers;
      E_of_showers.push_back(motherparticle->E());
      if(motherparticle->E() > fSimEnergyCut){
	++num_of_showers_viaEcut;
	++numshoowerspassenergy;

	E_of_showers_viaEcut.push_back(motherparticle->E());

	//using the RecoUtil function calculate the number of hits that see a charge deposition from the track.
	std::map<geo::PlaneID,int> Hit_num_map = RecoUtils::NumberofHitsThatContainEnergyDepositedByTracks(showermother->second, hits);

	//Calculaute the number of wires hit.
	std::map<geo::PlaneID,int> Wire_num_map = ShowerUtils::NumberofWiresHitByShower(showermother->second, hits);

	int low_density=0;

	//Compare hit density on the collection plane;
	for(std::map<geo::PlaneID,int>::iterator Hitnum_iter=Hit_num_map.begin(); Hitnum_iter!=Hit_num_map.end(); ++Hitnum_iter){
	  if(Wire_num_map[Hitnum_iter->first] == 0){continue;}
	  double Hit_num = (Hitnum_iter->second);
	  double Wire_num = Wire_num_map[Hitnum_iter->first];
	  double Hit_Density = Hit_num/Wire_num;
	  //check with dom
	  if(Hit_Density > fDensityCut){
	    ++num_of_showers_viaDensitycut;
	    ++numshowerspassdensity;
	    E_of_showers_viaDensitycut.push_back(motherparticle->E());

	    ++low_density;
	    break;
	  }
	}

      	//If we don't have a density bigger than one in at least one plane then it aint a shower. This could be due to hit reco.
	if(low_density == 0){
	  if(fVerbose > 0){std::cout << "Mother removed with id: " << showermother->first << " becuase the density is too low in the hit reconstruction" << std::endl;}
      	  showermother = ShowersMothers.erase(showermother);
	  continue;
	}

	//Should we remove shower mothers if the  they are not contained in one TPC and one TPC only
	if(fRemoveNonContainedParticles){
	  bool showercontained = true;
	  for(std::vector<int>::iterator showerdaughter=(showermother->second).begin(); showerdaughter!=(showermother->second).end(); ++showerdaughter){
	    showercontained = showercontained && mcparticlescontained[*showerdaughter];
	  }
	  if(!showercontained){
	    if(fVerbose > 0){std::cout << "Mother removed with id: " << showermother->first << " becuase it was not contained" << std::endl;}
	    showermother = ShowersMothers.erase(showermother);
	    ++containmentCutNum;
	    std::cout<<motherparticle->E()<<std::endl;
	    std::cout<<"Event Killed by containment cut"<<std::endl;
	    return; // Comment this out if you only want to remove showers which are not contained, this removes whole events
	    continue;
	  }
	}

	++numshowerspassTPC;
	++showermother;
      }
      else {
	if(fVerbose > 0){std::cout << "Mother removed with id: " << showermother->first << " becuase the true energy is too low" << std::endl;}
	showermother = ShowersMothers.erase(showermother);
	continue;
      }
    }
    else {
      if(fVerbose > 0){std::cout << "Mother removed with id: " << showermother->first << " becuase it is not a electron or photon" << std::endl;}
      showermother = ShowersMothers.erase(showermother);
    }
  }


  //If there are no true showers we can't validate
  if(ShowersMothers.size() == 0){
    if(fVerbose > 0){std::cout << " No Mothers finishing" << std::endl;}
    return;
  }

  //If we are looking at a pion then fill the smallest and largest showers
  for (sim::ParticleList::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt){

    const simb::MCParticle *particle = particleIt->second;

    float SmallestShowerE_Before = 999;
    float BiggestShowerE_Before  = -999;
    float SmallestShowerE_After = 999;
    float BiggestShowerE_After  = -999;

    if(particle->PdgCode() == 111){

      int NumDaughters = particle->NumberDaughters();
      for(int daughter=0; daughter<NumDaughters; ++daughter){

	//Get the daughter track ID
	int DaughterID = particle->Daughter(daughter);

	//Get the largest and smallest energies before the cuts
	if(trueParticles[DaughterID]->E() < SmallestShowerE_Before){SmallestShowerE_Before= trueParticles[DaughterID]->E();}
	if(trueParticles[DaughterID]->E() > BiggestShowerE_Before) {BiggestShowerE_Before = trueParticles[DaughterID]->E();}

	//Get the largest and smallest energies after the cuts
	if(ShowersMothers.find(DaughterID) != ShowersMothers.end()){
	  if(trueParticles[DaughterID]->E() < SmallestShowerE_After){SmallestShowerE_After = trueParticles[DaughterID]->E();}
	  if(trueParticles[DaughterID]->E() > BiggestShowerE_After) {BiggestShowerE_After  = trueParticles[DaughterID]->E();}
	}
      }
    }

    if(particle->PdgCode() == 11 && particle->TrackId() == 0){
      SmallestShowerE_Before  = particle->E();
      BiggestShowerE_Before = particle->E();
      SmallestShowerE_After = particle->E();
      BiggestShowerE_After  = particle->E();
    }

  }

  //Get the MC Energy deposited for each MC track.
  std::map<int,float> MCTrack_Energy_map = RecoUtils::TrueEnergyDepositedFromMCTracks(simchannels);

  //Get the number of hits associated wit each Track. This is done for every hits handle in the event.
  std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > > MCTrack_hit_map;
  for(auto& handle : hitHandles){

    if(!handle.isValid()){
      mf::LogError("ShowerValidation") << "Bad hit handle from the all the hit handles" << std::endl;
      continue;
    }

    //RawHitFinder is a bit silly at the moment so ignore it for the time being - REMOVE
    if(handle.provenance()->moduleLabel() == "fasthit"){continue;}

    //Getting the Hit Information
    art::Handle<std::vector<recob::Hit> > hitHandle;
    std::vector<art::Ptr<recob::Hit> > hits_fromhandle;
    if(evt.getByLabel(handle.provenance()->moduleLabel(),hitHandle))
      {art::fill_ptr_vector(hits_fromhandle, hitHandle);}


    //Get a map of the number of hits per plane each track has deposited.
    MCTrack_hit_map[handle.id()] = RecoUtils::NumberofPlaneHitsPerTrack(hits_fromhandle);
  }

  //######################
  //### Hit Validation ###
  //######################

  //Loop over the hit handles
  for(unsigned int hitlab_it=0; hitlab_it<fHitModuleLabels.size(); ++hitlab_it){

    //Getting the Hit Information
    art::Handle<std::vector<recob::Shower> > hitValHandle;
    std::vector<art::Ptr<recob::Shower> > hits_fromValhandle;
    std::string fHitModuleLabel = fHitModuleLabels[hitlab_it];
    if(evt.getByLabel(fHitModuleLabel,hitValHandle))
      {art::fill_ptr_vector(hits_fromValhandle,hitValHandle);}

    //Calculate the total energy deposited
    float TotalEnergyDeposited = 0;
    for(auto track_iter: MCTrack_Energy_map){
      TotalEnergyDeposited += track_iter.second;
    }

    //Calculate how much energy was deposited in the hits
    for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){

      float TotalEnergyDepinHits = RecoUtils::TotalEnergyDepinHits(hits,plane_id.Plane);
      if(TotalEnergyDepinHits != 0){
	HitCompletenessEnergy_TreeVal[fHitModuleLabel][plane_id.Plane].push_back((TotalEnergyDepinHits)/TotalEnergyDeposited);
	if(fVerbose > 1){std::cout << "Hit Completeness:"  << (TotalEnergyDepinHits)/TotalEnergyDeposited << std::endl;}
      }
      else{
	HitCompletenessEnergy_TreeVal[fHitModuleLabel][plane_id.Plane].push_back(-99999);
      }
    }
  }


  if(fVerbose > 1){std::cout << "Hit Validation Complete" << std::endl;}

  //##############################################
  //Get the reconstructed info and Match with MC #
  //##############################################

  for(unsigned int shwrlab_it=0; shwrlab_it<fShowerModuleLabels.size(); ++shwrlab_it){

    if(ShowersMothers.size() == 0){
      if(fVerbose > 0){std::cout << "No Mothers to Match to the showers" << std::endl;}
      continue;
    }

    //Getting the Shower Information
    art::Handle<std::vector<recob::Shower> > showerListHandle;
    std::vector<art::Ptr<recob::Shower> > showers;
    std::string fShowerModuleLabel = fShowerModuleLabels[shwrlab_it];
    if(evt.getByLabel(fShowerModuleLabel,showerListHandle))
      {art::fill_ptr_vector(showers,showerListHandle);}

    if(showers.size() == 0){
      if(fVerbose){std::cout << "No Shower in the Event" << std::endl;}
      continue;
    }

    art::Handle<std::vector<recob::PFParticle> > pfpListHandle;
    std::vector<art::Ptr<recob::PFParticle> > pfps;
    if(evt.getByLabel(fPFParticleLabel,pfpListHandle))
      {art::fill_ptr_vector(pfps,pfpListHandle);}


    //Getting the Shower Information
    //Association between Showers and 2d Hits
    art::FindManyP<recob::Hit> fmh(showerListHandle, evt, fShowerModuleLabel);

    //Association between Showers and clusters
    art::FindManyP<recob::Cluster> fmch(showerListHandle, evt, fShowerModuleLabel);

    //Association between pfParticle and clusters
    art::FindManyP<recob::Cluster> fmpfc(pfpListHandle, evt, fPFParticleLabel);

    //Association between pfParticle and vertex
    art::FindManyP<recob::Vertex> fmpfv(pfpListHandle, evt, fPFParticleLabel);

    //Association between Showers and pfParticle
    art::FindManyP<recob::PFParticle> fmpf(showerListHandle, evt, fShowerModuleLabel);

    std::vector< art::Ptr<recob::Hit> > showerhits; //hits in the shower
    std::vector< art::Ptr<recob::Vertex> > neutrinoVertices;
    //######################
    //### PFP Validation ###
    //######################

    if (fPFPValidation){

      //Create a map between PFParticles and their IDs
      std::map<int, art::Ptr<recob::PFParticle> > pfpsMap;
      for (unsigned int i=0; i<pfps.size();++i){
	art::Ptr<recob::PFParticle>& pfp = pfps.at(i);
	pfpsMap[pfp->Self()] = pfp;
      }

      std::vector<int> pfpPrimaries;
      int pfpTrackCounter        = 0;
      int pfpShowerCounter       = 0;
      int pfpNeutrinoCounter     = 0;
      int pfpShowerVertexCounter = 0;
      if(fVerbose > 0) { std::cout<<"Number of PFP: "<<pfps.size()<<std::endl; };
      for (unsigned int pfp_iter=0; pfp_iter<pfps.size();++pfp_iter){
	art::Ptr<recob::PFParticle>& pfp = pfps.at(pfp_iter);

	if(fVerbose > 1){
	  std::cout<<"PFParticle: "<<pfp->Self()<<" with: PDG: "<<pfp->PdgCode()
		   <<" Parent: "<<pfp->Parent()<<std::endl;
	}
	if ((pfp->PdgCode()==12) ||(pfp->PdgCode()==14)){ //Find Neutrino and primary daughters
	  // Gives a vector of Duaghter particle ID's, get Daughters using PFParticlesMap
	  const std::vector<size_t> Daughters = pfp->Daughters();
	  for (unsigned int daughter_iter=0; daughter_iter< Daughters.size(); daughter_iter++) {
	    art::Ptr<recob::PFParticle>& Daughter =pfpsMap[Daughters.at(daughter_iter)];
	    if(fmpfc.isValid()){
	      art::Handle<std::vector<recob::Cluster > > clusterHandle;
	      evt.get(fmpfc.at(0).front().id(),clusterHandle);
	      if(clusterHandle.isValid()){
		std::vector< art::Ptr<recob::Cluster> > pfpClusters = fmpfc.at(Daughter.key());
		ana::ShowerValidation::PFPValidation(pfpClusters,Daughter,evt,clusterHandle,ShowersMothers,MCTrack_Energy_map,MCTrack_hit_map,fShowerModuleLabel);
	      }
	    }

	    // Split into shower like, 11, and track like, 13
	    if (Daughter->PdgCode() == 11) {
	      pfpPrimaries.push_back(Daughter->Self());
	      ++pfpShowerCounter;

	      art::Handle<std::vector<recob::Vertex > > vertexHandle;
	      evt.get(fmpfv.at(0).front().id(),vertexHandle);

	      if(vertexHandle.isValid()) {
		std::vector< art::Ptr<recob::Vertex> > pfpVertexVector = fmpfv.at(Daughter.key());
		pfpShowerVertexCounter = pfpVertexVector.size();
		PFP_Showers_Vertices_TreeVal[fShowerModuleLabel].push_back(pfpShowerVertexCounter);
      
	      }
	    }else if (Daughter->PdgCode() == 13) {
	      pfpPrimaries.push_back(Daughter->Self());
	      ++pfpTrackCounter;
	    }else {
	      std::cout<<"Something has gone horribly wrong, PFP PDG != 11||13"<<std::endl;
	    }
	  } // end loop over neutrino daughters

	  if(fmpfv.isValid()) {
	    art::Handle<std::vector<recob::Vertex > > vertexHandle;
	    evt.get(fmpfv.at(0).front().id(),vertexHandle);

	    if(vertexHandle.isValid()) {
	      std::vector< art::Ptr<recob::Vertex> > pfpVertexVector = fmpfv.at(pfp.key());
	      art::Ptr<recob::Vertex> pfpVertex = pfpVertexVector.at(0);
	      neutrinoVertices.push_back(pfpVertex);
	    }
	  }

	  ++pfpNeutrinoCounter;
	}
      }

      if(fVerbose > 0) { std::cout<<"Primary Tracks: "<<pfpTrackCounter<<" and Primary Showers: "<<pfpShowerCounter<<std::endl; };

      PFP_Num_Neutrinos_TreeVal[fShowerModuleLabel].push_back(pfpNeutrinoCounter);
      PFP_Tracks_TreeVal[fShowerModuleLabel].push_back(pfpTrackCounter);
      PFP_Showers_TreeVal[fShowerModuleLabel].push_back(pfpShowerCounter);
      if(fVerbose > 0) { std::cout<<"Number of PFP Neurtinos:"<<pfpNeutrinoCounter
				  <<" and vertex: "<<neutrinoVertices.size()<<std::endl; };
      //TODO: Make this more general than the vertex case
      for (unsigned int vertex_iter=0; vertex_iter<neutrinoVertices.size();++vertex_iter){
	art::Ptr<recob::Vertex> pfpVertex = neutrinoVertices.at(vertex_iter);

	// get the pfp neutrino vertex
	double pfpvtx[3];
	pfpVertex->XYZ(pfpvtx);

	// TODO: here i just take the first particle becuase i'm lazy, needs fixing
	std::map<int,const simb::MCParticle*>::iterator trueInitialParticleIt = trueInitialParticles.begin();
	const simb::MCParticle* initialParticle = (*trueInitialParticleIt).second;


	const TLorentzVector PositionTrajP = initialParticle->Position();
	double truevtx[3] = {PositionTrajP.X() ,PositionTrajP.Y(), PositionTrajP.Z()};

	double PFP_Start_Diff_X = pfpvtx[0] - truevtx[0];
	double PFP_Start_Diff_Y = pfpvtx[1] - truevtx[1];
	double PFP_Start_Diff_Z = pfpvtx[2] - truevtx[2];
	double PFP_Start_Diff = TMath::Sqrt(TMath::Power((pfpvtx[0] - truevtx[0]),2) + TMath::Power((pfpvtx[1] - truevtx[1]),2) + TMath::Power((pfpvtx[2] - truevtx[2]),2));

	PFP_Vertex_Distance_X_TreeVal[fShowerModuleLabel].push_back(PFP_Start_Diff_X);
	PFP_Vertex_Distance_Y_TreeVal[fShowerModuleLabel].push_back(PFP_Start_Diff_Y);
	PFP_Vertex_Distance_Z_TreeVal[fShowerModuleLabel].push_back(PFP_Start_Diff_Z);
	PFP_Vertex_Distance_Mag_TreeVal[fShowerModuleLabel].push_back(PFP_Start_Diff);

      }
    }

    //Get the ID of the shower hit module
    art::ProductID showerhit_productid = fmh.at(0).front().id();

    unsigned int max_hitnum=0;
    unsigned int biggest_shower_iter = 9999;

    for(unsigned int shower_iter = 0; shower_iter < showers.size(); ++shower_iter){

      ++numrecoshowers;

      //Get the shower
      art::Ptr<recob::Shower>& shower = showers.at(shower_iter);

      //Get the hits vector from the shower
      showerhits = fmh.at(shower.key());

      if(showerhits.size() > max_hitnum){
	max_hitnum = showerhits.size();
	biggest_shower_iter = shower_iter;
      }
    }

    std::map<int,float> MinStartDiff;
    std::map<int,float> MinStartX;
    std::map<int,float> MinStartY;
    std::map<int,float> MinStartZ;
    std::map<int,float> MinShowerDirection_Xdiff;
    std::map<int,float> MinShowerDirection_Ydiff;
    std::map<int,float> MinShowerDirection_Zdiff;
    std::map<int,float> MinShowerDirection_diff;
    std::map<int,int>   MinEvaluateShowerDirection;
    std::map<int,int>   MinEvaluateShowerStart;

    //Initialise the minimum values for each shower mother
    for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end();++showermother){
      MinStartDiff[showermother->first] = 99999;
      MinStartX[showermother->first] = -99999;
      MinStartY[showermother->first] = -99999;
      MinStartZ[showermother->first] = -99999;
      MinShowerDirection_Xdiff[showermother->first] = -99999;
      MinShowerDirection_Ydiff[showermother->first] = -99999;
      MinShowerDirection_Zdiff[showermother->first] = -99999;
      MinShowerDirection_diff[showermother->first] = -99999;
      MinEvaluateShowerDirection[showermother->first] = 0;
      MinEvaluateShowerStart[showermother->first] = 0;
    }

    //Loop over the showers in the event
    for(unsigned int shower_iter = 0; shower_iter < showers.size(); ++shower_iter){

      if(fUseBiggestShower == true){
	if(shower_iter != biggest_shower_iter){continue;}
      }

      //Get the shower
      art::Ptr<recob::Shower>& shower = showers.at(shower_iter);


      //#########################
      //### Shower Validation ###
      //#########################

      //Get the hits vector from the shower
      showerhits = fmh.at(shower.key());


      if((int) showerhits.size() < fMinHitSize) {continue;}

      int ShowerBest_Plane = shower->best_plane();
      if(std::isnan(ShowerBest_Plane)){ShowerBest_Plane = 0;}

      //Function from RecoUtils, finds the most probable track ID associated with the set of hits from there true energy depositons. The pair returns the energy as well.
      std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain(ShowersMothers,showerhits,ShowerBest_Plane);

      int ShowerTrackID = ShowerTrackInfo.first;
      double TrueEnergyDepWithinShower_FromTrueShower = ShowerTrackInfo.second;

      //Check to see if the shower was correctly matched.
      if(TrueEnergyDepWithinShower_FromTrueShower == -99999){
	continue;
	std::cout << "Reco Shower not matched to true shower. Think of reducing the energy threshold" << std::endl;
      }

      //Get the number of hits associated to the true particle.
      int TrueHitDep_FromTrueShower = 0;
      for(std::vector<int>::iterator track_id=ShowersMothers[ShowerTrackID].begin(); track_id!=ShowersMothers[ShowerTrackID].end(); ++track_id){
	for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){
	  TrueHitDep_FromTrueShower +=  MCTrack_hit_map[showerhit_productid][*track_id][plane_id];
	}
      }

      int NumberofHitsinRecoShower = showerhits.size();

      //Get the number of hits in the reco shower from the true shower.
      int TrueHitsDep_WithinRecoShower = 0;
      std::map<int,std::map<geo::PlaneID,int> >  MCTrack_showerhit_map = RecoUtils::NumberofPlaneHitsPerTrack(showerhits);
      for(auto const& planehit_map : MCTrack_showerhit_map){
	if(std::find(ShowersMothers[ShowerTrackID].begin(), ShowersMothers[ShowerTrackID].end(),planehit_map.first) == ShowersMothers[ShowerTrackID].end()){continue;}
	for(auto const& hit_num : planehit_map.second){
	  TrueHitsDep_WithinRecoShower += hit_num.second;
	}
      }


      double TrueEnergyDep_FromShower = 0;
      //Calculate the true Energy deposited By Shower
      for(std::vector<int>::iterator daughterID=ShowersMothers[ShowerTrackID].begin(); daughterID!=ShowersMothers[ShowerTrackID].end(); ++daughterID){
	TrueEnergyDep_FromShower += MCTrack_Energy_map[*daughterID];
      }

      double hitcompleteness = 0;
      if(TrueEnergyDep_FromShower != 0){
	hitcompleteness = ((double) TrueHitsDep_WithinRecoShower)/(double) TrueHitDep_FromTrueShower;
	ShowerHitsCompleteness_TreeVal[fShowerModuleLabel].push_back(hitcompleteness);
      }

      double hitpurity = 0;
      if(TrueHitsDep_WithinRecoShower != 0){
	hitpurity   =   (double) TrueHitsDep_WithinRecoShower/(double) NumberofHitsinRecoShower;
	ShowerHitsPurity_TreeVal[fShowerModuleLabel].push_back(hitpurity);
      }

      //Energy deposited within the set of Hits associated to the shower.
      double TrueEnergyDep_WithinRecoShower = 0;

      //Loop over the hits and find the IDEs
      if(fVerbose > 1){std::cout << "shower hits size: " << showerhits.size() << std::endl;}
      for(std::vector< art::Ptr<recob::Hit> >::iterator hitIt=showerhits.begin(); hitIt!=showerhits.end(); ++hitIt){

	//Get the plane ID
	geo::WireID wireid = (*hitIt)->WireID();
	int PlaneID = wireid.Plane;

	if(PlaneID != ShowerBest_Plane){continue;}

	//Split the Hit into its IDE for each track it associates with.
	std::vector<sim::TrackIDE> trackIDEs = backtracker->HitToTrackIDEs((*hitIt));
	for (unsigned int idIt = 0; idIt < trackIDEs.size(); ++idIt){

	  //Find the true total energy deposited in a set of hits.
	  TrueEnergyDep_WithinRecoShower += trackIDEs.at(idIt).energy;
	}
      }//Hit Loop

      double energycompleteness = 0;
      if(TrueEnergyDep_FromShower != 0){
	energycompleteness =  (TrueEnergyDepWithinShower_FromTrueShower)/TrueEnergyDep_FromShower;
	ShowerEnergyCompleteness_TreeVal[fShowerModuleLabel].push_back(energycompleteness);
      }

      double energypurity = 0;
      if(TrueEnergyDep_WithinRecoShower != 0){
	energypurity       =  TrueEnergyDepWithinShower_FromTrueShower/TrueEnergyDep_WithinRecoShower;
	ShowerEnergyPurity_TreeVal[fShowerModuleLabel].push_back(energypurity);
      }

      //Find the MCParticle this shower associates to
      const simb::MCParticle* MCShowerParticle = trueParticles.at(ShowerTrackID);

      //Find the Energy of the particle:
      float Energy = MCShowerParticle->E();
      ShowerTrueEnergy_TreeVal[fShowerModuleLabel].push_back(Energy*1000);

      ShowerStartEndProcess_TreeVal[fShowerModuleLabel].push_back(MCShowerParticle->EndProcess());

      //Get the number of Traj points to loop over
      unsigned int TrajPoints = MCShowerParticle->NumberTrajectoryPoints();

      // Select first traj point where the photon loses energy, last be default
      bool PhotonEnd = false;
      unsigned int PhotonEndTrajPoint = TrajPoints-1;
      for (unsigned int trajPoint=0; trajPoint<TrajPoints;trajPoint++){
	if (!PhotonEnd && MCShowerParticle->E(trajPoint) < (0.9*Energy)){
	  PhotonEndTrajPoint = trajPoint;
	  PhotonEnd = true;
	}
      }


      //Find the start and end points of the initial particle.
      TLorentzVector PositionTrajStart =  MCShowerParticle->Position(0);
      TLorentzVector PositionTrajEnd   =  MCShowerParticle->Position(PhotonEndTrajPoint);

      //The Start of position for the electron shower is PositionTrajStart but the start position for photon showers is at the end of the shower (the start of the e+- track)
      if(MCShowerParticle->PdgCode() == 22){
	PositionTrajStart = PositionTrajEnd;
      }

      //The three vector for track length is the shower direction
      //TVector3  TrueShowerDirection = (PositionTrajEnd - PositionTrajStart).Vect();
      TVector3  TrueShowerDirection(MCShowerParticle->Px(), MCShowerParticle->Py(),MCShowerParticle->Pz());

      //Initial track lentgh of the shower.
      double TrueTrackLength = TrueShowerDirection.Mag();

      //Get the information for the shower
      TVector3  ShowerDirection                  = shower->Direction();
      TVector3  ShowerStart                      = shower->ShowerStart();//cm
      double ShowerTrackLength                   = shower->Length();//cm
      std::vector< double >  ShowerEnergyPlanes  = shower->Energy();//MeV
      std::vector< double >  ShowerdEdX_vec      = shower->dEdx();//MeV/cm

      //Remove this
      //ShowerEnergyPlanes[2] = (ShowerEnergyPlanes[2] - 0.00155171)*0.00155171/4.39964 + 4.39964;

      //Bools to fill metric histrograms wheen needed.
      bool EvaluateShowerDirection       = false;
      bool EvaluateShowerStart           = false;
      bool EvaluateShowerLength          = false;
      bool EvaluateShowerEnergy          = false;
      bool EvaluateShowerdEdx            = false;
      bool EvalulateGeoProjectionMatched = false;

      //Evaulate 3D Shower Reconstruction Dependent Metrics
      if(!std::isnan(ShowerDirection.X()) || ShowerDirection.Mag() == 0) {EvaluateShowerDirection = true; ++MinEvaluateShowerDirection[ShowerTrackID];}
      if(!std::isnan(ShowerStart.X()))                                   {EvaluateShowerStart     = true; ++MinEvaluateShowerStart[ShowerTrackID];}
      if(!std::isnan(ShowerTrackLength) || shower->has_length())         {EvaluateShowerLength    = true;}
      if(ShowerEnergyPlanes.size() != 0){
	if(!std::isnan(ShowerEnergyPlanes.at(0))){
	  EvaluateShowerEnergy = true;
	}
      }
      if(ShowerdEdX_vec.size() != 0){
	if(!std::isnan(ShowerdEdX_vec.at(0))){
	  EvaluateShowerdEdx = true;
	}
      }

      //Get the angles between the direction
      float ShowerDirection_Xdiff = -99999;
      float ShowerDirection_Ydiff = -99999;
      float ShowerDirection_Zdiff = -99999;
      float ShowerDirection_diff  = -99999;

      if(TMath::Sqrt((TrueShowerDirection.Y()*TrueShowerDirection.Y() + TrueShowerDirection.Z()*TrueShowerDirection.Z()))*TMath::Sqrt((ShowerDirection.Y()*ShowerDirection.Y() + ShowerDirection.Z()*ShowerDirection.Z())) !=0){
	ShowerDirection_Xdiff = (TrueShowerDirection.Y()*ShowerDirection.Y() + TrueShowerDirection.Z()*ShowerDirection.Z())/(TMath::Sqrt((TrueShowerDirection.Y()*TrueShowerDirection.Y() + TrueShowerDirection.Z()*TrueShowerDirection.Z()))*TMath::Sqrt((ShowerDirection.Y()*ShowerDirection.Y() + ShowerDirection.Z()*ShowerDirection.Z())));
      }

      if(TMath::Sqrt((TrueShowerDirection.X()*TrueShowerDirection.X() + TrueShowerDirection.Z()*TrueShowerDirection.Z()))*TMath::Sqrt((ShowerDirection.X()*ShowerDirection.X() + ShowerDirection.Z()*ShowerDirection.Z())) !=0){
	ShowerDirection_Ydiff = (TrueShowerDirection.X()*ShowerDirection.X() + TrueShowerDirection.Z()*ShowerDirection.Z())/(TMath::Sqrt((TrueShowerDirection.X()*TrueShowerDirection.X() + TrueShowerDirection.Z()*TrueShowerDirection.Z()))*TMath::Sqrt((ShowerDirection.X()*ShowerDirection.X() + ShowerDirection.Z()*ShowerDirection.Z())));
      }

      if(TMath::Sqrt((TrueShowerDirection.Y()*TrueShowerDirection.Y() + TrueShowerDirection.X()*TrueShowerDirection.X()))*TMath::Sqrt((ShowerDirection.Y()*ShowerDirection.Y() + ShowerDirection.X()*ShowerDirection.X())) !=0){
	ShowerDirection_Zdiff = (TrueShowerDirection.Y()*ShowerDirection.Y() + TrueShowerDirection.X()*ShowerDirection.X())/(TMath::Sqrt((TrueShowerDirection.Y()*TrueShowerDirection.Y() + TrueShowerDirection.X()*TrueShowerDirection.X()))*TMath::Sqrt((ShowerDirection.Y()*ShowerDirection.Y() + ShowerDirection.X()*ShowerDirection.X())));
      }

      if(TrueShowerDirection.Mag() != 0 || ShowerDirection.Mag() !=0){
	ShowerDirection_diff  = TrueShowerDirection.Dot(ShowerDirection)/(TrueShowerDirection.Mag()*ShowerDirection.Mag());
      }

      //Get the Error in the position. intialised as 0,0,0 this is a problem here.
      double Start_diff = TMath::Sqrt(TMath::Power(PositionTrajStart.X()-ShowerStart.X(),2) + TMath::Power(PositionTrajStart.Y()-ShowerStart.Y(),2) + TMath::Power(PositionTrajStart.Z()-ShowerStart.Z(),2));

      //Fill the histograms.
      if(EvaluateShowerDirection){
	ShowerDirection_X_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_Xdiff);
	ShowerDirection_Y_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_Ydiff);
	ShowerDirection_Z_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_Zdiff);
	ShowerDirectionDiff_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_diff);
      }
      else{
	ShowerDirection_X_TreeVal[fShowerModuleLabel].push_back(-99999);
	ShowerDirection_Y_TreeVal[fShowerModuleLabel].push_back(-99999);
	ShowerDirection_Z_TreeVal[fShowerModuleLabel].push_back(-99999);
	ShowerDirectionDiff_TreeVal[fShowerModuleLabel].push_back(-99999);
      }

      if(EvaluateShowerStart){
	ShowerStart_X_TreeVal[fShowerModuleLabel].push_back(TMath::Abs(PositionTrajStart.X()-ShowerStart.X()));
	ShowerStart_Y_TreeVal[fShowerModuleLabel].push_back(TMath::Abs(PositionTrajStart.Y()-ShowerStart.Y()));
	ShowerStart_Z_TreeVal[fShowerModuleLabel].push_back(TMath::Abs(PositionTrajStart.Z()-ShowerStart.Z()));
	ShowerMag_TreeVal[fShowerModuleLabel].push_back(Start_diff);
      }
      else{
	ShowerStart_X_TreeVal[fShowerModuleLabel].push_back(-99999);
	ShowerStart_Y_TreeVal[fShowerModuleLabel].push_back(-99999);
	ShowerStart_Z_TreeVal[fShowerModuleLabel].push_back(-99999);
	ShowerMag_TreeVal[fShowerModuleLabel].push_back(-99999);
      }

      if(EvaluateShowerLength){
	ShowerLength_TreeVal[fShowerModuleLabel].push_back(TMath::Abs(TrueTrackLength-ShowerTrackLength));
      }
      else{
	ShowerLength_TreeVal[fShowerModuleLabel].push_back(-99999);
      }

      if(TrueEnergyDep_FromShower != 0 && EvaluateShowerEnergy){
 	ShowerEnergyDiff_TreeVal[fShowerModuleLabel].push_back(ShowerEnergyPlanes[ShowerBest_Plane]/TrueEnergyDep_FromShower);
 	ShowerTotalEnergyDiff_TreeVal[fShowerModuleLabel].push_back((ShowerEnergyPlanes[ShowerBest_Plane] - TrueEnergyDep_FromShower)/TrueEnergyDep_FromShower);
      }
      else{
	ShowerEnergyDiff_TreeVal[fShowerModuleLabel].push_back(-99999);
 	ShowerTotalEnergyDiff_TreeVal[fShowerModuleLabel].push_back(-99999);
      }

      if(TrueEnergyDepWithinShower_FromTrueShower != 0 && EvaluateShowerEnergy){
	ShowerRecoEnergyVsTrueEnergyinRecoShower_TreeVal[fShowerModuleLabel].push_back(ShowerEnergyPlanes[ShowerBest_Plane]/TrueEnergyDepWithinShower_FromTrueShower);
      }
      else{
	ShowerRecoEnergyVsTrueEnergyinRecoShower_TreeVal[fShowerModuleLabel].push_back(-99999);
      }

      if(EvaluateShowerEnergy){
	ShowerEnergy_TreeVal[fShowerModuleLabel].push_back(ShowerEnergyPlanes[ShowerBest_Plane]);
      }
      else{
	ShowerEnergy_TreeVal[fShowerModuleLabel].push_back(-99999);
      }

      ShowerHitNum_TreeVal[fShowerModuleLabel].push_back(showerhits.size());

      if(EvaluateShowerdEdx){
	ShowerdEdx_TreeVal[fShowerModuleLabel].push_back((ShowerdEdX_vec[ShowerBest_Plane]));
      }
      else{
	ShowerdEdx_TreeVal[fShowerModuleLabel].push_back(-99999);
      }

      //Fill the 2D histograms

      if(fVerbose > 0){
	std::cout << "#################################################" << std::endl;
	std::cout << "Global Event Information" << std::endl;
	std::cout << "#################################################" << std::endl;
	std::cout << "Shower Label: " <<  fShowerModuleLabel << std::endl;
	std::cout << "Number of Reco showers: " << showers.size() << std::endl;
	std::cout << "Energy Simulated: " << simenergy << std::endl;
	std::cout << "Number of True Showers Before Cuts: " <<  num_of_showers << std::endl;
	std::cout << "Number of True Showers That pass the E cut: " <<  num_of_showers_viaEcut << std::endl;
	std::cout << "Number of True Showers That pass the Density cut: " << num_of_showers_viaDensitycut << std::endl;
	std::cout << "Number of True Showers in the Event: " << ShowersMothers.size() << std::endl;
	std::cout << "#################################################" << std::endl;
	std::cout << "Reconstructed/Associated Truth Information" << std::endl;
	std::cout << "#################################################" << std::endl;
	std::cout << "Hit Size: " << showerhits.size() << std::endl;
	std::cout << "ShowerBest_Plane: " << ShowerBest_Plane << std::endl;
	if(EvaluateShowerDirection){
	  std::cout << "True Start: " << PositionTrajStart.X() << " Shower Start: " << ShowerStart.X() << std::endl;
	  std::cout << "X Poisition: " <<  ShowerStart.X() << "Y Position " << ShowerStart.Y() << " Z Poistion: " << ShowerStart.Z() << std::endl;
	}
	if(EvaluateShowerLength){
	  std::cout << "TrueTrackLength: " << TrueTrackLength << " ShowerTrackLength: " << ShowerTrackLength << std::endl;
	}
	if(EvaluateShowerEnergy){
	  std::cout << "Best Plane Reco Shower Energy " << ShowerEnergyPlanes[ShowerBest_Plane] << std::endl;
	}
	std::cout << "True Energy Deposited from true shower within Shower: " << TrueEnergyDepWithinShower_FromTrueShower << std::endl;
	std::cout << "True Energy Deposited From true associated Shower: " << TrueEnergyDep_FromShower << std::endl;
	std::cout << "True Energy Deposited by hits in shower: " <<  TrueEnergyDep_WithinRecoShower << std::endl;
	std::cout << "Energy Purity: " << energypurity << " Energy completeness: " << energycompleteness << std::endl;
	std::cout << "Hit Purity: " << hitpurity << "Hit Completeness: " << hitcompleteness << std::endl;
	std::cout << "#################################################" <<std::endl;
      }

      if(fVerbose > 1){std::cout << "Shower Validation Complete" << std::endl;}


      //##########################
      //### Cluster Validation ###
      //##########################

      //Firstly look inf the hits in the shower match in all given planes
      std::map<std::vector<double>, int> HitCoord_map;
      float numclusters = -9999;
      float geoprojectionmatched_score = -99999;
      for(std::vector< art::Ptr<recob::Hit> >::iterator hitIt=showerhits.begin(); hitIt!=showerhits.end(); ++hitIt){
	try{
	  const std::vector<sim::IDE> hitIDEs = backtracker->HitToAvgSimIDEs(*hitIt);
	  for(unsigned int hitIDE=0; hitIDE<hitIDEs.size(); ++hitIDE){
	    std::vector<double> hitcoord = {hitIDEs[hitIDE].x, hitIDEs[hitIDE].y, hitIDEs[hitIDE].z};
	    //std::cout << "x: " << hitIDEs[hitIDE].x << ", y: " << hitIDEs[hitIDE].y << ", z: " << hitIDEs[hitIDE].z << " plane: " << (*hitIt)->WireID().Plane << std::endl;;
	    ++HitCoord_map[hitcoord];
	  }
	}
	catch(...){
	  if(fVerbose>1){std::cout << "Noise Hit" << std::endl;}
	}
      }


      //Get the clusters associated to the shower.
      if(fmch.isValid()){

	art::Handle<std::vector<recob::Cluster > > clusterHandle;
	evt.get(fmch.at(shower.key()).front().id(),clusterHandle);
	if(clusterHandle.isValid()){
	  std::vector<art::Ptr<recob::Cluster> > showerclusters = fmch.at(shower.key());
	  ana::ShowerValidation::ClusterValidation(showerclusters,evt,clusterHandle,ShowersMothers,MCTrack_Energy_map,MCTrack_hit_map,ShowerTrackID,simenergy,fShowerModuleLabel);

	  //Loop over the hit coordinate map where there there is a hit on every plane give a 1.
	  numclusters = showerclusters.size();
	  int geoprojectionmatched = 0;
	  for(std::map<std::vector<double>, int>::iterator coord=HitCoord_map.begin(); coord!=HitCoord_map.end(); ++coord){
	    if(coord->second == numclusters){++geoprojectionmatched;}
	  }
	  if(numclusters > 0){
	    EvalulateGeoProjectionMatched = true;
	    geoprojectionmatched_score = (float) geoprojectionmatched/(float) HitCoord_map.size();
	  }
	}
	else{
	  mf::LogError("ShowerValidation") << "Cluster handle is stale. No clustering validation done" << std::endl;
	}
      }
      else if(fmpf.isValid()){

	//Find the Clusters associated to PF particle.
	art::Handle<std::vector<recob::PFParticle> > pfpHandle;
	evt.get(fmpf.at(shower.key()).front().id(),pfpHandle);
	if(pfpHandle.isValid()){
	  art::FindManyP<recob::Cluster> fmcpf(pfpHandle, evt, pfpHandle.provenance()->moduleLabel());
	  if(fmcpf.isValid()){
	    art::Handle<std::vector<recob::Cluster > > clusterHandle;
	    evt.get(fmcpf.at(0).front().id(),clusterHandle);
	    if(clusterHandle.isValid()){
	      std::vector< art::Ptr<recob::Cluster> > showerclusters = fmcpf.at(shower.key());
	      ana::ShowerValidation::ClusterValidation(showerclusters,evt,clusterHandle,ShowersMothers,MCTrack_Energy_map,MCTrack_hit_map,ShowerTrackID,simenergy,fShowerModuleLabel);
	      //Loop over the hit coordinate map where there there is a hit on every plane give a 1.
	      numclusters = showerclusters.size();

	      int geoprojectionmatched = 0;
	      for(std::map<std::vector<double>, int>::iterator coord=HitCoord_map.begin(); coord!=HitCoord_map.end(); ++coord){
		if(coord->second == numclusters){++geoprojectionmatched;}
	      }

	      if(numclusters > 0){
		EvalulateGeoProjectionMatched = true;
		geoprojectionmatched_score = (float) geoprojectionmatched/(float) HitCoord_map.size();
	      }
	    }
	    else{
	      mf::LogError("ShowerValidation") << "Cluster handle is stale. No clustering validation done" << std::endl;
	    }
	  }
	  else{
	    mf::LogError("ShowerValidation") << "No Assosoication between pf particles and clusters was found for shower made from a pf particle. No clustering validation was done." << std::endl;
	  }
	}
	else{
	  mf::LogError("ShowerValidation") << "pf particle handle is stale" << std::endl;
	}
      }
      else{
	mf::LogError("ShowerValidation") << "No cluster or pandora association" << std::endl;
      }

      if(EvalulateGeoProjectionMatched){

	GeoProjectionMatched_TreeVal[fShowerModuleLabel].push_back(geoprojectionmatched_score);

      }
      else{
	GeoProjectionMatched_TreeVal[fShowerModuleLabel].push_back(-99999);
      }

      ++numrecoshowersana;

      if(fVerbose > 1){std::cout << "Cluster Validation Complete" << std::endl;}

    }//Shower Loop

    //#########################
    //###   Event Metric    ###
    //#########################

    //Calculate the True Hit number
    for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end(); ++showermother){
      int TrueHitDep_FromTrueShowers = 0;
      for(std::vector<int>::iterator track_id=(showermother->second).begin(); track_id!=(showermother->second).end(); ++track_id){
	for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){
	  TrueHitDep_FromTrueShowers +=  MCTrack_hit_map[showerhit_productid][*track_id][plane_id];
	}
      }
    }

    //This is rather pandora track specific for tuning pandora so lets add it if the shower module used pandora and the track module Pandora track is used. Very horrid and sorry for that.
    PandoraTrackNum_TreeVal[fShowerModuleLabel].push_back(tracks.size());

    Event_NumTrueShowers_TreeVal[fShowerModuleLabel].push_back(num_of_showers);
    Event_NumTrueShowers_EnergyCut_TreeVal[fShowerModuleLabel].push_back(num_of_showers_viaEcut);
    Event_NumTrueShowers_DensityCut_TreeVal[fShowerModuleLabel].push_back(num_of_showers_viaDensitycut);
    for (uint i=0; i<E_of_showers.size(); i++){
      Event_TrueShowerE_TreeVal[fShowerModuleLabel].push_back(E_of_showers.at(i));
    }
    for(uint i=0; i<E_of_showers_viaEcut.size(); i++){
      Event_TrueShowerE_EnergyCut_TreeVal[fShowerModuleLabel].push_back(E_of_showers_viaEcut.at(i));
    }
    for(uint i=0; i<E_of_showers_viaDensitycut.size(); i++){
      Event_TrueShowerE_DensityCut_TreeVal[fShowerModuleLabel].push_back(E_of_showers_viaDensitycut.at(i));
    }

    //Whats the segementyness of the event.
    EventSeggy_TreeVal[fShowerModuleLabel].push_back((float)showers.size()/(float)num_of_showers_viaDensitycut);
    TrueEnergy_TreeVal[fShowerModuleLabel].push_back(simenergy*1000);



  }//Shower Module labels

  //Fill the tree
  Tree->Fill();

  for(unsigned int shwrlab_it=0; shwrlab_it<fShowerModuleLabels.size(); ++shwrlab_it){
    ShowerDirection_X_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerDirection_Y_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerDirection_Z_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerStart_X_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerStart_Y_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerStart_Z_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerLength_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerEnergyDiff_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerdEdx_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    EventSeggy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerEnergyCompleteness_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerHitsPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerHitsCompleteness_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerEnergyPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerEnergy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerHitNum_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerTotalEnergyDiff_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerMag_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerDirectionDiff_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerRecoEnergyVsTrueEnergyinRecoShower_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerTrueEnergy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    TrueEnergy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    TrueHitNum_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PosDir_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PandoraTrackNum_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();

    PFP_Num_Neutrinos_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_Vertex_Distance_X_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_Vertex_Distance_Y_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_Vertex_Distance_Z_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_Vertex_Distance_Mag_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_Tracks_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_Showers_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_Showers_Vertices_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_ProjectionMatched_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_HitsCompleteness_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_EnergyCompleteness_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_HitsPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_EnergyPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_Track_ProjectionMatched_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_Track_HitsCompleteness_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_Track_EnergyCompleteness_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_Track_HitsPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_Track_EnergyPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_Shower_ProjectionMatched_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_Shower_HitsCompleteness_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_Shower_EnergyCompleteness_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_Shower_HitsPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PFP_Shower_EnergyPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    Event_NumTrueShowers_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    Event_NumTrueShowers_EnergyCut_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    Event_NumTrueShowers_DensityCut_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    Event_TrueShowerE_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    Event_TrueShowerE_EnergyCut_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    Event_TrueShowerE_DensityCut_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();

    ClusterProjectionMatchedEnergy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ClusterCompletenessEnergy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ClusterPurityEnergy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ClusterCompletenessHits_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ClusterPurityHits_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ClusterCompPurityEnergy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ClusterCompPurityHits_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();

  }

  for(unsigned int hitlab_it=0; hitlab_it<fHitModuleLabels.size(); ++hitlab_it){
    for(unsigned int plane_it=0; plane_it<geom->Nplanes(); ++plane_it){
      HitCompletenessEnergy_TreeVal[fHitModuleLabels[hitlab_it]][plane_it].clear();
    }
  }

  return;
}


void ana::ShowerValidation::ClusterValidation(std::vector< art::Ptr<recob::Cluster> >& clusters,
					      const art::Event& evt,
					      art::Handle<std::vector<recob::Cluster> >& clusterHandle,
					      std::map<int,std::vector<int> >& ShowerMotherTrackIDs,
					      std::map<int,float>& MCTrack_Energy_map,
					      std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > >& MCTrack_hit_map,
					      int& TrueShowerID,
					      float& simenergy,
					      std::string& fShowerModuleLabel){


  //Initialise Trees
  std::vector<std::vector<float> > ClusterProjectionMatchedEnergy_TreeVec((int)geom->Nplanes());
  std::vector<std::vector<float> > ClusterCompletenessEnergy_TreeVec((int)geom->Nplanes());
  std::vector<std::vector<float> > ClusterPurityEnergy_TreeVec((int)geom->Nplanes());
  std::vector<std::vector<float> > ClusterCompletenessHits_TreeVec((int)geom->Nplanes());
  std::vector<std::vector<float> > ClusterPurityHits_TreeVec((int)geom->Nplanes());
  std::vector<std::vector<float> > ClusterCompPurityEnergy_TreeVec((int)geom->Nplanes());
  std::vector<std::vector<float> > ClusterCompPurityHits_TreeVec((int)geom->Nplanes());

  //Get the associated hits
  art::FindManyP<recob::Hit> fmhc(clusterHandle, evt, clusterHandle.provenance()->moduleLabel());

  //Holder for cluster its
  std::vector< art::Ptr<recob::Hit> > clusterhits;
  art::Handle<std::vector<recob::Hit > > hitHandle;

  if(!hitHandle.isValid()){
    mf::LogError("ShowerValidation") << "Hits handle is stale. No clustering validation done" << std::endl;
    return;
  }

  //Get the Hits Handle used for this cluster type WARNING
  evt.get(fmhc.at(clusters.at(0).key()).front().id(),hitHandle);

  //Get the hits vector from the shower
  for(auto const& cluster : clusters){

    clusterhits = fmhc.at(cluster.key());
    //Function from RecoUtils, finds the most probable track ID associated with the set of hits from there true energy depositons. The pair returns the energy in the hits as well.
    std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain(ShowerMotherTrackIDs,clusterhits, cluster->Plane().Plane);

    //Make sure the cluster has been matched.
    if(ShowerTrackInfo.second == -99999){
      std::cout << "Reco cluster not matched to a True shower" << std::endl;
      continue;
    }

    float TotalTrueEnergy = 0;
    float signalhits      = 0;
    float totalhits       = 0;

    for(std::vector<int>::iterator daughterID=ShowerMotherTrackIDs[ShowerTrackInfo.first].begin(); daughterID!=ShowerMotherTrackIDs[ShowerTrackInfo.first].end(); ++daughterID){

      //Calculate the true Energy deposited By Shower
      TotalTrueEnergy += MCTrack_Energy_map[*daughterID];
      // Get the map of planeIDs to number of hits
      std::map<geo::PlaneID, int> hitPlaneMap = RecoUtils::NumberofPlaneHitsFromTrack(*daughterID, clusterhits);

      for(std::map<geo::PlaneID, int>::iterator hitPlaneMapit = hitPlaneMap.begin();  hitPlaneMapit != hitPlaneMap.end();  hitPlaneMapit++){
	std::cout<<"Plane ID: "<<(*hitPlaneMapit).first<<std::endl;
	std::cout<<"SignalHits: "<<(*hitPlaneMapit).second<<std::endl;
	std::cout<<"TotalHits: "<<MCTrack_hit_map[hitHandle.id()][*daughterID][(*hitPlaneMapit).first]<<std::endl;

	//Count how many hits are from the true shower.
	signalhits += (*hitPlaneMapit).second;
	//Count how many hits are missed in the plane from the true track id.
	totalhits += MCTrack_hit_map[hitHandle.id()][*daughterID][(*hitPlaneMapit).first];
      }
    }

    int projection_match = -99999;

    //Have we matched the 2D cluster to the correct shower correctly. In terms of Energy depositions:
    if(ShowerTrackInfo.first == TrueShowerID){projection_match = 1;}
    else{projection_match = 0;}

    //Calculate the purity and completeness metrics
    float completeness_hits   = -99999;
    float purity_hits         = -99999;
    float completeness_energy = -99999;
    float purity_energy       = -99999;

    float TotalEnergyDepinHits = RecoUtils::TotalEnergyDepinHits(clusterhits,cluster->Plane().Plane);

    ClusterProjectionMatchedEnergy_TreeVec[cluster->Plane().Plane].push_back(projection_match);

    if(totalhits != 0){
      completeness_hits = (signalhits)/totalhits;
      ClusterCompletenessHits_TreeVec[cluster->Plane().Plane].push_back(completeness_hits);
    }

    if(clusterhits.size() != 0){
      purity_hits = signalhits/clusterhits.size();
      ClusterPurityHits_TreeVec[cluster->Plane().Plane].push_back(purity_hits);

    }

    if(TotalTrueEnergy != 0){
      completeness_energy = (ShowerTrackInfo.second)/TotalTrueEnergy;
      ClusterCompletenessEnergy_TreeVec[cluster->Plane().Plane].push_back(completeness_energy);
    }

    if(TotalEnergyDepinHits != 0){
      purity_energy = ShowerTrackInfo.second/TotalEnergyDepinHits;
      ClusterPurityEnergy_TreeVec[cluster->Plane().Plane].push_back(purity_energy);
    }

    if(totalhits != 0 && clusterhits.size() != 0){
      ClusterCompPurityHits_TreeVec[cluster->Plane().Plane].push_back(completeness_hits*purity_hits);
    }

    if(TotalTrueEnergy != 0 && TotalEnergyDepinHits != 0){
      ClusterCompPurityEnergy_TreeVec[cluster->Plane().Plane].push_back(completeness_energy*purity_energy);
    }

    if(fVerbose>1){
      std::cout << "#################################################"    << std::endl;
      std::cout << "               Cluster Metrics                   "    << std::endl;
      std::cout << "#################################################"    << std::endl;
      std::cout << "Projection matched:          " << projection_match    << std::endl;
      std::cout << "Cluster hit completeness:    " << completeness_hits   << std::endl;
      std::cout << "Cluster hit purity:          " << purity_hits         << std::endl;
      std::cout << "Cluster energy completeness: " << completeness_energy << std::endl;
      std::cout << "Cluster energy purity:       " << purity_energy       << std::endl;
      std::cout << "#################################################"    << std::endl;
    }
  }//Cluster Loop

  ClusterProjectionMatchedEnergy_TreeVal[fShowerModuleLabel].push_back(ClusterProjectionMatchedEnergy_TreeVec);
  ClusterCompletenessEnergy_TreeVal[fShowerModuleLabel].push_back(ClusterCompletenessEnergy_TreeVec);
  ClusterPurityEnergy_TreeVal[fShowerModuleLabel].push_back(ClusterPurityEnergy_TreeVec);
  ClusterCompletenessHits_TreeVal[fShowerModuleLabel].push_back(ClusterCompletenessHits_TreeVec);
  ClusterPurityHits_TreeVal[fShowerModuleLabel].push_back(ClusterPurityHits_TreeVec);
  ClusterCompPurityEnergy_TreeVal[fShowerModuleLabel].push_back(ClusterCompPurityEnergy_TreeVec);
  ClusterCompPurityHits_TreeVal[fShowerModuleLabel].push_back(ClusterCompPurityHits_TreeVec);

  return;
}


void ana::ShowerValidation::PFPValidation(std::vector<art::Ptr<recob::Cluster> >& clusters,
					  art::Ptr<recob::PFParticle>  pfp,
					  const art::Event& evt,
					  art::Handle<std::vector<recob::Cluster> >& clusterHandle,
					  std::map<int,std::vector<int> >& ShowerMotherTrackIDs,
					  std::map<int,float>& MCTrack_Energy_map,
					  std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > > & MCTrack_hit_map,
					  std::string & fShowerModuleLabel){


  //Get the associated hits
  art::FindManyP<recob::Hit> fmhc(clusterHandle, evt, clusterHandle.provenance()->moduleLabel());

  //Holder for cluster hits
  std::vector< art::Ptr<recob::Hit> > clusterhits;

  //Get the Hits Handle used for this cluster type
  art::Handle<std::vector<recob::Hit > > hitHandle;
  evt.get(fmhc.at(clusters.at(0).key()).front().id(),hitHandle);

  if(!hitHandle.isValid()){
    mf::LogError("ShowerValidation") << "Hits handle is stale. No pfp validation done" << std::endl;
    return;
  }

  float TotalTrueEnergy      = 0;
  float TotalEnergyDepinHits = 0;
  float signalhits           = 0;
  float totalhits            = 0;
  float pfphits              = 0;
  float pfpenergy            = 0;

  //Calculate the purity and completeness metrics
  float completeness_hits   = -99999;
  float purity_hits         = -99999;
  float completeness_energy = -99999;
  float purity_energy       = -99999;

  int projectionMatched = 1;
  std::vector<int> projected_IDs;

  //Get the hits vector from the shower
  for(auto const& cluster : clusters){
    clusterhits = fmhc.at(cluster.key());
    // if shower-like use showerutils, else if track-like use recoutils
    if (pfp->PdgCode() == 11) {
      //Function from RecoUtils, finds the most probable track ID associated with the set of hits from there true energy depositons. The pair returns the energy in the hits as well.
      std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain(ShowerMotherTrackIDs,clusterhits, cluster->Plane().Plane);

      //Make sure the cluster has been matched.
      if(ShowerTrackInfo.second == -99999){
	std::cout << "Reco cluster not matched to a True shower" << std::endl;
	continue;
      }

      for(std::vector<int>::iterator daughterID=ShowerMotherTrackIDs[ShowerTrackInfo.first].begin(); daughterID!=ShowerMotherTrackIDs[ShowerTrackInfo.first].end(); ++daughterID){

	//Calculate the true Energy deposited By Shower
	TotalTrueEnergy += MCTrack_Energy_map[*daughterID];

	// Get the map of planeIDs to number of hits
	std::map<geo::PlaneID, int> hitPlaneMap = RecoUtils::NumberofPlaneHitsFromTrack(*daughterID, clusterhits);

	for(std::map<geo::PlaneID, int>::iterator hitPlaneMapit = hitPlaneMap.begin();  hitPlaneMapit != hitPlaneMap.end();  hitPlaneMapit++){
	  //std::cout<<"Plane ID: "<<(*hitPlaneMapit).first<<std::endl;
       	  //std::cout<<"SignalHits: "<<(*hitPlaneMapit).second<<std::endl;
	  //std::cout<<"TotalHits: "<<MCTrack_hit_map[hitHandle.id()][*daughterID][(*hitPlaneMapit).first]<<std::endl;

	  //Count how many hits are from the true shower.
	  signalhits += (*hitPlaneMapit).second;
	  //Count how many hits are missed in the plane from the true track id.
	  totalhits += MCTrack_hit_map[hitHandle.id()][*daughterID][(*hitPlaneMapit).first];

	}
      }
      projected_IDs.push_back(ShowerTrackInfo.first);

      TotalEnergyDepinHits += RecoUtils::TotalEnergyDepinHits(clusterhits,cluster->Plane().Plane);

      pfphits   += clusterhits.size();
      pfpenergy += ShowerTrackInfo.second;
    }else if (pfp->PdgCode() == 13) {
      int PFPTrackInfo = RecoUtils::TrueParticleIDFromTotalRecoHits(clusterhits); //check which recoutil to use with dom

      if(PFPTrackInfo == -99999){
	std::cout << "Reco cluster not matched to a True shower" << std::endl;
	continue;
      }

      //Calculate the true Energy deposited By Shower
      TotalTrueEnergy += MCTrack_Energy_map[PFPTrackInfo];

      // Get the map of planeIDs to number of hits
      std::map<geo::PlaneID, int> hitPlaneMap = RecoUtils::NumberofPlaneHitsFromTrack(PFPTrackInfo, clusterhits);

      for(std::map<geo::PlaneID, int>::iterator hitPlaneMapit = hitPlaneMap.begin();  hitPlaneMapit != hitPlaneMap.end();  hitPlaneMapit++){
	std::cout<<"Plane ID: "<<(*hitPlaneMapit).first<<std::endl;
	std::cout<<"SignalHits: "<<(*hitPlaneMapit).second<<std::endl;
	std::cout<<"TotalHits: "<<MCTrack_hit_map[hitHandle.id()][PFPTrackInfo][(*hitPlaneMapit).first]<<std::endl;

	//Count how many hits are from the true shower.
	signalhits += (*hitPlaneMapit).second;
	//Count how many hits are missed in the plane from the true track id.
	totalhits += MCTrack_hit_map[hitHandle.id()][PFPTrackInfo][(*hitPlaneMapit).first];

      }


      //Calculate the true Energy deposited By Shower
      TotalTrueEnergy += MCTrack_Energy_map[PFPTrackInfo];


      projected_IDs.push_back(PFPTrackInfo);

      TotalEnergyDepinHits += RecoUtils::TotalEnergyDepinHits(clusterhits,cluster->Plane().Plane);

      pfphits      += clusterhits.size();
      pfpenergy    += RecoUtils::TotalEnergyDepinHitsFromTrack(clusterhits, PFPTrackInfo, cluster->Plane().Plane);//check with dom

    }else {
      std::cout<<"Something has gone horribly wrong, PFP PDG != 11||13"<<std::endl;
    }
  }

  //Have we matched the 2D cluster to the correct shower correctly. In terms of Energy depositions:
  //if(ShowerTrackInfo.first == TrueShowerID){projection_match = 1;}
  //else{projection_match = 0;}

  for (auto ID : projected_IDs){
    if (ID != projected_IDs.at(0)) {projectionMatched=0;}
  }

  PFP_ProjectionMatched_TreeVal[fShowerModuleLabel].push_back(projectionMatched);
  if (pfp->PdgCode() == 11) { PFP_Shower_ProjectionMatched_TreeVal[fShowerModuleLabel].push_back(projectionMatched);
  }else if (pfp->PdgCode() == 13) { PFP_Track_ProjectionMatched_TreeVal[fShowerModuleLabel].push_back(projectionMatched);
  }

  //ClusterProjectionMatchedEnergy_TreeVec[cluster->Plane().Plane].push_back(projection_match);
  if(totalhits != 0){
    completeness_hits = (signalhits)/totalhits;
    PFP_HitsCompleteness_TreeVal[fShowerModuleLabel].push_back(completeness_hits);
    if (pfp->PdgCode() == 11) { PFP_Shower_HitsCompleteness_TreeVal[fShowerModuleLabel].push_back(completeness_hits);
    }else if (pfp->PdgCode() == 13) { PFP_Track_HitsCompleteness_TreeVal[fShowerModuleLabel].push_back(completeness_hits);
    }
  }

  if(pfphits != 0){
    purity_hits = signalhits/pfphits;
    PFP_HitsPurity_TreeVal[fShowerModuleLabel].push_back(purity_hits);
    if (pfp->PdgCode() == 11) { PFP_Shower_HitsPurity_TreeVal[fShowerModuleLabel].push_back(purity_hits);
    }else if (pfp->PdgCode() == 13) { PFP_Track_HitsPurity_TreeVal[fShowerModuleLabel].push_back(purity_hits);
    }
  }

  if(TotalTrueEnergy != 0){
    completeness_energy = pfpenergy/TotalTrueEnergy;
    PFP_EnergyCompleteness_TreeVal[fShowerModuleLabel].push_back(completeness_energy);
    if (pfp->PdgCode() == 11) { PFP_Shower_EnergyCompleteness_TreeVal[fShowerModuleLabel].push_back(completeness_energy);
    }else if (pfp->PdgCode() == 13) { PFP_Track_EnergyCompleteness_TreeVal[fShowerModuleLabel].push_back(completeness_energy);
    }
  }

  if(TotalEnergyDepinHits != 0){
    purity_energy = pfpenergy/TotalEnergyDepinHits;
    PFP_EnergyPurity_TreeVal[fShowerModuleLabel].push_back(purity_energy);
    if (pfp->PdgCode() == 11) { PFP_Shower_EnergyPurity_TreeVal[fShowerModuleLabel].push_back(purity_energy);
    }else if (pfp->PdgCode() == 13) { PFP_Track_EnergyPurity_TreeVal[fShowerModuleLabel].push_back(purity_energy);
    }
  }


  if(fVerbose>1){
    std::cout << "#################################################"    << std::endl;
    std::cout << "                 PFP Metrics                     "    << std::endl;
    std::cout << "#################################################"    << std::endl;
    std::cout << "Projection matched:          " << projectionMatched   << std::endl;
    std::cout << "PFP hit completeness:    " << completeness_hits       << std::endl;
    std::cout << "PFP hit purity:          " << purity_hits             << std::endl;
    std::cout << "PFP energy completeness: " << completeness_energy     << std::endl;
    std::cout << "PFP energy purity:       " << purity_energy           << std::endl;
    std::cout << "#################################################"    << std::endl;
  }

  // need to add position difference

  return;
}

void ana::ShowerValidation::endJob() {

  std::cout << "Number of events ran over: " <<  numevents  << std::endl;
  std::cout << "Number of events removed by the containment cut: " << containmentCutNum <<std::endl;
  std::cout << "Number of initial MC Showers: " <<  numshowers  << std::endl;
  std::cout << "Number of MC showers that pass the energy cut: " <<  numshoowerspassenergy  << std::endl;
  std::cout << "Number of MC showers that pass the density cut: " << numshowerspassdensity  << std::endl;
  std::cout << "Number of MC showers that pass the containment cut: " <<  numshowerspassTPC  << std::endl;
  std::cout << "Number of reco showers: " << numrecoshowers  << std::endl;
  std::cout << "Number of reco showers analysed: " << numrecoshowersana  << std::endl;

}



DEFINE_ART_MODULE(ana::ShowerValidation)
