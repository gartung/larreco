#ifndef CELLTREETRUTH_MODULE
#define CELLTREETRUTH_MODULE

// LArSoft includes
#include "lardataobj/MCBase/MCStep.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "larcore/Geometry/Geometry.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes                                                              
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TUnixSystem.h"
#include "TDatabasePDG.h"
#include "TObjArray.h"
#include "TTimeStamp.h"
#include "TVector3.h"
#include "TObjArray.h"

// C++ Includes                                                                
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdio>

using namespace std;

namespace wc {

  class CellTreeTruth : public art::EDAnalyzer{
  public:

    explicit CellTreeTruth(fhicl::ParameterSet const& pset);
    virtual ~CellTreeTruth();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);
    void reconfigure(fhicl::ParameterSet const& pset);
    void initOutput();
    void reset();

    void processMCNeutrino(const art::Event& evt);
    void processMCParticle(const art::Event& evt);
    void processMCTrack(const art::Event& evt);
    void processMCShower(const art::Event& evt);
    void processSimChannel(const art::Event& evt);
    void processBeamOpFlash(const art::Event& evt);
    void processCosmicOpFlash(const art::Event& evt);

  private:

    std::string fOutFileName;

    TFile *fOutFile;
    TTree *fEventTree;

    int fRun;
    int fSubRun;
    int fEvent;
    double fEventTime;

    // MC Neutrino;
    bool fSaveMCNeutrino;
    std::string fMCNeutrinoLabel;
    TTree *fNeutrinoTree;
    vector<int> mcneutrino_nuStatus;
    vector<int> mcneutrino_nuTrackId;
    vector<int> mcneutrino_nuPdg;
    vector<int> mcneutrino_nuMother;
    vector<string> mcneutrino_nuProcess;
    vector<string> mcneutrino_nuEndProcess;
    vector<double> mcneutrino_nuMass;
    vector<double> mcneutrino_nuWeight;
    vector<double> mcneutrino_nuVertexX;
    vector<double> mcneutrino_nuVertexY;
    vector<double> mcneutrino_nuVertexZ;
    vector<double> mcneutrino_nuEnergy;
    //vector<double> mcneutrino_nuGvtxX;
    //vector<double> mcneutrino_nuGvtxY;
    //vector<double> mcneutrino_nuGvtxZ;
    vector<int> mcneutrino_nuRescatter;
    vector<int> mcneutrino_leptonStatus;
    vector<int> mcneutrino_leptonTrackId;
    vector<int> mcneutrino_leptonPdg;
    vector<int> mcneutrino_leptonMother;
    vector<string> mcneutrino_leptonProcess;
    vector<string> mcneutrino_leptonEndProcess;
    vector<double> mcneutrino_leptonMass;
    vector<double> mcneutrino_leptonWeight;
    vector<double> mcneutrino_leptonVertexX;
    vector<double> mcneutrino_leptonVertexY;
    vector<double> mcneutrino_leptonVertexZ;
    vector<double> mcneutrino_leptonEnergy;
    //vector<double> mcneutrino_leptonGvtxX;
    //vector<double> mcneutrino_leptonGvtxY;
    //vector<double> mcneutrino_leptonGvtxZ;
    vector<int> mcneutrino_leptonRescatter;
    vector<int> mcneutrino_ccnc;
    vector<int> mcneutrino_mode;
    vector<int> mcneutrino_interactionType;
    vector<int> mcneutrino_target;
    vector<int> mcneutrino_hitNuc;
    vector<int> mcneutrino_hitQuark;
    vector<double> mcneutrino_w;
    vector<double> mcneutrino_x;
    vector<double> mcneutrino_y;
    vector<double> mcneutrino_qSqr;
    vector<double> mcneutrino_pt;
    vector<double> mcneutrino_theta;

    // MC PARTICLE
    bool fSaveMCParticle;
    std::string fMCParticleLabel;
    TTree *fParticleTree;
    int Nmcparticle;
    vector<int> mcparticle_id;
    vector<int> mcparticle_statusCode;
    vector<int> mcparticle_pdg;
    vector<int> mcparticle_mother;
    vector<std::string> mcparticle_process;
    vector<std::string> mcparticle_endProcess;
    vector<int> mcparticle_ndaughters;
    vector<int> mcparticle_daughterId;
    vector<int> mcparticle_ntrajpts;
    vector<double> mcparticle_x;                                      
    vector<double> mcparticle_y;                              
    vector<double> mcparticle_z;                              
    vector<double> mcparticle_t;
    vector<double> mcparticle_px;
    vector<double> mcparticle_py;
    vector<double> mcparticle_pz;
    vector<double> mcparticle_e;
    vector<double> mcparticle_p;
    vector<double> mcparticle_pt;                              
    vector<double> mcparticle_endX;
    vector<double> mcparticle_endY;
    vector<double> mcparticle_endZ;
    vector<double> mcparticle_endT;
    vector<double> mcparticle_mass;
    vector<double> mcparticle_endPX;
    vector<double> mcparticle_endPY;
    vector<double> mcparticle_endPZ;
    vector<double> mcparticle_endE;
    vector<double> mcparticle_GvtxX;
    vector<double> mcparticle_GvtxY;
    vector<double> mcparticle_GvtxZ;
    vector<double> mcparticle_GvtxPX;
    vector<double> mcparticle_GvtxPY;
    vector<double> mcparticle_GvtxPZ;
    vector<double> mcparticle_GvtxE;
    vector<int> mcparticle_rescatter;
    vector<double> mcparticle_weight;

    // MC TRACK
    bool fSaveMCTrack;
    std::string fMCTrackLabel;
    TTree *fTrackTree;
    int Nmctrack;
    vector<int> mctrack_pdg;
    vector<int> mctrack_motherPdg;
    vector<int> mctrack_ancestorPdg;
    vector<int> mctrack_id;
    vector<int> mctrack_motherId;
    vector<int> mctrack_ancestorId;
    vector<std::string> mctrack_process;
    vector<std::string> mctrack_motherProcess;
    vector<std::string> mctrack_ancestorProcess;
    vector<double> mctrack_startX;
    vector<double> mctrack_startY;
    vector<double> mctrack_startZ;
    vector<double> mctrack_startT;
    vector<double> mctrack_motherStartX;
    vector<double> mctrack_motherStartY;
    vector<double> mctrack_motherStartZ;
    vector<double> mctrack_motherStartT;
    vector<double> mctrack_ancestorStartX;
    vector<double> mctrack_ancestorStartY;
    vector<double> mctrack_ancestorStartZ;
    vector<double> mctrack_ancestorStartT;
    vector<double> mctrack_endX;
    vector<double> mctrack_endY;
    vector<double> mctrack_endZ;
    vector<double> mctrack_endT;
    vector<double> mctrack_motherEndX;
    vector<double> mctrack_motherEndY;
    vector<double> mctrack_motherEndZ;
    vector<double> mctrack_motherEndT;
    vector<double> mctrack_ancestorEndX;
    vector<double> mctrack_ancestorEndY;
    vector<double> mctrack_ancestorEndZ;
    vector<double> mctrack_ancestorEndT;
    vector<double> mctrack_startPX;
    vector<double> mctrack_startPY;
    vector<double> mctrack_startPZ;
    vector<double> mctrack_startE;
    vector<double> mctrack_motherStartPX;
    vector<double> mctrack_motherStartPY;
    vector<double> mctrack_motherStartPZ;
    vector<double> mctrack_motherStartE;
    vector<double> mctrack_ancestorStartPX;
    vector<double> mctrack_ancestorStartPY;
    vector<double> mctrack_ancestorStartPZ;
    vector<double> mctrack_ancestorStartE;
    vector<double> mctrack_endPX;
    vector<double> mctrack_endPY;
    vector<double> mctrack_endPZ;
    vector<double> mctrack_endE;
    vector<double> mctrack_motherEndPX;
    vector<double> mctrack_motherEndPY;
    vector<double> mctrack_motherEndPZ;
    vector<double> mctrack_motherEndE;
    vector<double> mctrack_ancestorEndPX;
    vector<double> mctrack_ancestorEndPY;
    vector<double> mctrack_ancestorEndPZ;
    vector<double> mctrack_ancestorEndE;

    // MC SHOWER
    bool fSaveMCShower;
    std::string fMCShowerLabel;
    TTree *fShowerTree;
    int Nmcshower;
    vector<int> mcshower_pdg;
    vector<int> mcshower_motherPdg;
    vector<int> mcshower_ancestorPdg;
    vector<int> mcshower_id;
    vector<int> mcshower_motherId;
    vector<int> mcshower_ancestorId;
    vector<std::string> mcshower_process;
    vector<std::string> mcshower_motherProcess;
    vector<std::string> mcshower_ancestorProcess;
    vector<double> mcshower_startX;
    vector<double> mcshower_startY;
    vector<double> mcshower_startZ;
    vector<double> mcshower_startT;
    vector<double> mcshower_motherStartX;
    vector<double> mcshower_motherStartY;
    vector<double> mcshower_motherStartZ;
    vector<double> mcshower_motherStartT;
    vector<double> mcshower_ancestorStartX;
    vector<double> mcshower_ancestorStartY;
    vector<double> mcshower_ancestorStartZ;
    vector<double> mcshower_ancestorStartT;
    vector<double> mcshower_endX;
    vector<double> mcshower_endY;
    vector<double> mcshower_endZ;
    vector<double> mcshower_endT;
    vector<double> mcshower_motherEndX;
    vector<double> mcshower_motherEndY;
    vector<double> mcshower_motherEndZ;
    vector<double> mcshower_motherEndT;
    vector<double> mcshower_ancestorEndX;
    vector<double> mcshower_ancestorEndY;
    vector<double> mcshower_ancestorEndZ;
    vector<double> mcshower_ancestorEndT;
    vector<double> mcshower_startPX;
    vector<double> mcshower_startPY;
    vector<double> mcshower_startPZ;
    vector<double> mcshower_startE;
    vector<double> mcshower_motherStartPX;
    vector<double> mcshower_motherStartPY;
    vector<double> mcshower_motherStartPZ;
    vector<double> mcshower_motherStartE;
    vector<double> mcshower_ancestorStartPX;
    vector<double> mcshower_ancestorStartPY;
    vector<double> mcshower_ancestorStartPZ;
    vector<double> mcshower_ancestorStartE;
    vector<double> mcshower_endPX;
    vector<double> mcshower_endPY;
    vector<double> mcshower_endPZ;
    vector<double> mcshower_endE;
    vector<double> mcshower_motherEndPX;
    vector<double> mcshower_motherEndPY;
    vector<double> mcshower_motherEndPZ;
    vector<double> mcshower_motherEndE;
    vector<double> mcshower_ancestorEndPX;
    vector<double> mcshower_ancestorEndPY;
    vector<double> mcshower_ancestorEndPZ;
    vector<double> mcshower_ancestorEndE;
    vector<double> mcshower_detProfileX;
    vector<double> mcshower_detProfileY;
    vector<double> mcshower_detProfileZ;
    vector<double> mcshower_detProfileT;
    vector<double> mcshower_detProfilePX;
    vector<double> mcshower_detProfilePY;
    vector<double> mcshower_detProfilePZ;
    vector<double> mcshower_detProfileE;
    vector<double> mcshower_dEdx;

    // SIM CHANNEL
    std::string fSimChannelLabel; 
    int Nsimchannel;
    vector<int> simchannel_channelID;
    vector<int> simchannel_tdc;
    vector<int> simchannel_nEnergyDeposits;
    vector<int> simchannel_id;
    vector<float> simchannel_nElectrons;
    vector<float> simchannel_energy;
    vector<float> simchannel_x;
    vector<float> simchannel_y;
    vector<float> simchannel_z; 

    // OPFLASH --- BEAM
    bool fBeamSaveOpFlash;
    std::string fBeamOpFlashLabel;
    TTree *fBeamOpFlashTree;
    float opBeamMultPEThresh; // fhicl param
    int ofBeam_nFlash;
    vector<float> ofBeam_t;
    vector<float> ofBeam_peTotal;
    vector<int> ofBeam_multiplicity;
    TClonesArray *fBeamPEperOpDet;

    // OPFLASH --- COSMIC
    bool fCosmicSaveOpFlash;
    std::string fCosmicOpFlashLabel;
    TTree *fCosmicOpFlashTree;
    float opCosmicMultPEThresh; // fhicl param
    int ofCosmic_nFlash;
    vector<float> ofCosmic_t;
    vector<float> ofCosmic_peTotal;
    vector<int> ofCosmic_multiplicity;
    TClonesArray *fCosmicPEperOpDet;

    art::ServiceHandle<geo::Geometry> fGeometry;

    // dummy variables not filled and not used in imaging
    int fCalib_nChannel;
    std::vector<int> fCalib_channelId;
    TClonesArray *fCalib_wf;
  }; // class

  //-------------------------------------------------------------------
  CellTreeTruth::CellTreeTruth(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
  {
    std::cout << "CellTreeTruth constructor" << std::endl;
    reconfigure(parameterSet);
    initOutput();
  }

  //-------------------------------------------------------------------
  CellTreeTruth::~CellTreeTruth()
  {
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::reconfigure(fhicl::ParameterSet const& p)
  {
    std::cout << "Reconfigure" << std::endl;
    fMCNeutrinoLabel = p.get<std::string>("MCNeutrinoLabel");
    fMCParticleLabel = p.get<std::string>("MCParticleLabel");
    fMCTrackLabel = p.get<std::string>("MCTrackLabel");
    fMCShowerLabel = p.get<std::string>("MCShowerLabel");
    fSimChannelLabel = p.get<std::string>("SimChannelLabel");
    fBeamOpFlashLabel = p.get<std::string>("BeamOpFlashLabel");
    fCosmicOpFlashLabel = p.get<std::string>("CosmicOpFlashLabel");

    fSaveMCNeutrino = p.get<bool>("saveMCNeutrino");
    fSaveMCParticle = p.get<bool>("saveMCParticle");
    fSaveMCTrack = p.get<bool>("saveMCTrack");
    fSaveMCShower = p.get<bool>("saveMCShower");
    //fSaveSimChannel = p.get<bool>("saveSimChannel");
    fBeamSaveOpFlash = p.get<bool>("saveBeamOpFlash");
    fCosmicSaveOpFlash = p.get<bool>("saveCosmicOpFlash");
   
    fOutFileName = p.get<std::string>("outFile");

    opBeamMultPEThresh   = p.get<float>("opBeamMultPEThresh");
    opCosmicMultPEThresh   = p.get<float>("opCosmicMultPEThresh");
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::initOutput()
  {
    std::cout << "Initialize output" << std::endl;
    TDirectory* tmpDir = gDirectory;
    fOutFile = new TFile(fOutFileName.c_str(), "recreate");

    TNamed version("version","4.0");
    version.Write();

    TDirectory* subDir = fOutFile->mkdir("Event");
    subDir->cd();
    fEventTree = new TTree("Sim","Event Tree from Simulation");
    fEventTree->Branch("runNo", &fRun);
    fEventTree->Branch("subRunNo", &fSubRun);
    fEventTree->Branch("eventNo", &fEvent);    
    fEventTree->Branch("eventTime", &fEventTime);

    // MC Neutrino
    fNeutrinoTree = new TTree("MCNeutrino","MC Neutrino Tree from Output of LArG4");
    fNeutrinoTree->Branch("mcneutrino_nuStatus",             &mcneutrino_nuStatus);
    fNeutrinoTree->Branch("mcneutrino_nuTrackId",            &mcneutrino_nuTrackId);
    fNeutrinoTree->Branch("mcneutrino_nuPdg",                &mcneutrino_nuPdg);
    fNeutrinoTree->Branch("mcneutrino_nuMother",             &mcneutrino_nuMother);
    fNeutrinoTree->Branch("mcneutrino_nuProcess",            &mcneutrino_nuProcess);
    fNeutrinoTree->Branch("mcneutrino_nuEndProcess",         &mcneutrino_nuEndProcess);
    fNeutrinoTree->Branch("mcneutrino_nuMass",               &mcneutrino_nuMass);
    fNeutrinoTree->Branch("mcneutrino_nuWeight",             &mcneutrino_nuWeight);
    fNeutrinoTree->Branch("mcneutrino_nuVertexX",            &mcneutrino_nuVertexX);
    fNeutrinoTree->Branch("mcneutrino_nuVertexY",            &mcneutrino_nuVertexY);
    fNeutrinoTree->Branch("mcneutrino_nuVertexZ",            &mcneutrino_nuVertexZ);
    fNeutrinoTree->Branch("mcneutrino_nuEnergy",             &mcneutrino_nuEnergy);     
    //fNeutrinoTree->Branch("mcneutrino_nuGvtxX",              &mcneutrino_nuGvtxX);
    //fNeutrinoTree->Branch("mcneutrino_nuGvtxY",              &mcneutrino_nuGvtxY);
    //fNeutrinoTree->Branch("mcneutrino_nuGvtxZ",              &mcneutrino_nuGvtxZ);
    fNeutrinoTree->Branch("mcneutrino_nuRescatter",          &mcneutrino_nuRescatter);
    fNeutrinoTree->Branch("mcneutrino_leptonStatus",         &mcneutrino_leptonStatus);
    fNeutrinoTree->Branch("mcneutrino_leptonTrackId",        &mcneutrino_leptonTrackId);
    fNeutrinoTree->Branch("mcneutrino_leptonPdg",            &mcneutrino_leptonPdg);
    fNeutrinoTree->Branch("mcneutrino_leptonMother",         &mcneutrino_leptonMother);
    fNeutrinoTree->Branch("mcneutrino_leptonProcess",        &mcneutrino_leptonProcess);
    fNeutrinoTree->Branch("mcneutrino_leptonEndProcess",     &mcneutrino_leptonEndProcess);
    fNeutrinoTree->Branch("mcneutrino_leptonMass",           &mcneutrino_leptonMass);
    fNeutrinoTree->Branch("mcneutrino_leptonWeight",         &mcneutrino_leptonWeight);
    fNeutrinoTree->Branch("mcneutrino_leptonVertexX",        &mcneutrino_leptonVertexX);
    fNeutrinoTree->Branch("mcneutrino_leptonVertexY",        &mcneutrino_leptonVertexY);
    fNeutrinoTree->Branch("mcneutrino_leptonVertexZ",        &mcneutrino_leptonVertexZ);
    fNeutrinoTree->Branch("mcneutrino_leptonEnergy",         &mcneutrino_leptonEnergy); 
    //fNeutrinoTree->Branch("mcneutrino_leptonGvtxX",          &mcneutrino_leptonGvtxX);
    //fNeutrinoTree->Branch("mcneutrino_leptonGvtxY",          &mcneutrino_leptonGvtxY);
    //fNeutrinoTree->Branch("mcneutrino_leptonGvtxZ",          &mcneutrino_leptonGvtxZ);
    fNeutrinoTree->Branch("mcneutrino_leptonRescatter",      &mcneutrino_leptonRescatter);
    fNeutrinoTree->Branch("mcneutrino_ccnc",                 &mcneutrino_ccnc);
    fNeutrinoTree->Branch("mcneutrino_mode",                 &mcneutrino_mode);
    fNeutrinoTree->Branch("mcneutrino_interactionType",      &mcneutrino_interactionType);
    fNeutrinoTree->Branch("mcneutrino_target",               &mcneutrino_target);
    fNeutrinoTree->Branch("mcneutrino_hitNuc",               &mcneutrino_hitNuc);
    fNeutrinoTree->Branch("mcneutrino_hitQuark",             &mcneutrino_hitQuark);
    fNeutrinoTree->Branch("mcneutrino_w",                    &mcneutrino_w);
    fNeutrinoTree->Branch("mcneutrino_x",                    &mcneutrino_x);
    fNeutrinoTree->Branch("mcneutrino_y",                    &mcneutrino_y);
    fNeutrinoTree->Branch("mcneutrino_qSqr",                 &mcneutrino_qSqr);
    fNeutrinoTree->Branch("mcneutrino_pt",                   &mcneutrino_pt);
    fNeutrinoTree->Branch("mcneutrino_theta",                &mcneutrino_theta);

    // MC PARTICLE
    fParticleTree = new TTree("MCParticle","MC Particle Tree from Output of LArG4");
    fParticleTree->Branch("Nmcparticle",               &Nmcparticle);
    fParticleTree->Branch("mcparticle_id",             &mcparticle_id);
    fParticleTree->Branch("mcparticle_statusCode",     &mcparticle_statusCode);
    fParticleTree->Branch("mcparticle_pdg",            &mcparticle_pdg);
    fParticleTree->Branch("mcparticle_mother",         &mcparticle_mother);
    fParticleTree->Branch("mcparticle_process",        &mcparticle_process);
    fParticleTree->Branch("mcparticle_endProcess",     &mcparticle_endProcess);           
    fParticleTree->Branch("mcparticle_ndaughters",     &mcparticle_ndaughters);
    fParticleTree->Branch("mcparticle_daughterId",     &mcparticle_daughterId);
    fParticleTree->Branch("mcparticle_ntrajpts",       &mcparticle_ntrajpts);
    fParticleTree->Branch("mcparticle_x",              &mcparticle_x);
    fParticleTree->Branch("mcparticle_y",              &mcparticle_y);
    fParticleTree->Branch("mcparticle_z",              &mcparticle_z);
    fParticleTree->Branch("mcparticle_t",              &mcparticle_t);
    fParticleTree->Branch("mcparticle_px",             &mcparticle_px);
    fParticleTree->Branch("mcparticle_py",             &mcparticle_py);
    fParticleTree->Branch("mcparticle_pz",             &mcparticle_pz);
    fParticleTree->Branch("mcparticle_e",              &mcparticle_e);
    fParticleTree->Branch("mcparticle_p",              &mcparticle_p);
    fParticleTree->Branch("mcparticle_pt",             &mcparticle_pt);
    fParticleTree->Branch("mcparticle_endX",           &mcparticle_endX);
    fParticleTree->Branch("mcparticle_endY",           &mcparticle_endY);
    fParticleTree->Branch("mcparticle_endZ",           &mcparticle_endZ);
    fParticleTree->Branch("mcparticle_endT",           &mcparticle_endT);
    fParticleTree->Branch("mcparticle_mass",           &mcparticle_mass);
    fParticleTree->Branch("mcparticle_endPX",          &mcparticle_endPX);
    fParticleTree->Branch("mcparticle_endPY",          &mcparticle_endPY);
    fParticleTree->Branch("mcparticle_endPZ",          &mcparticle_endPZ);
    fParticleTree->Branch("mcparticle_endE",           &mcparticle_endE);
    fParticleTree->Branch("mcparticle_GvtxX",          &mcparticle_GvtxX);
    fParticleTree->Branch("mcparticle_GvtxY",          &mcparticle_GvtxY);
    fParticleTree->Branch("mcparticle_GvtxZ",          &mcparticle_GvtxZ);
    fParticleTree->Branch("mcparticle_GvtxPX",         &mcparticle_GvtxPX);
    fParticleTree->Branch("mcparticle_GvtxPY",         &mcparticle_GvtxPY);
    fParticleTree->Branch("mcparticle_GvtxPZ",         &mcparticle_GvtxPZ);
    fParticleTree->Branch("mcparticle_GvtxE",          &mcparticle_GvtxE);
    fParticleTree->Branch("mcparticle_rescatter",      &mcparticle_rescatter);
    fParticleTree->Branch("mcparticle_weight",         &mcparticle_weight);

    // MC TRACK
    fTrackTree = new TTree("MCTrack","MC Track Tree from Output of LArG4");
    fTrackTree->Branch("Nmctrack",                  &Nmctrack);
    fTrackTree->Branch("mctrack_pdg",               &mctrack_pdg);
    fTrackTree->Branch("mctrack_motherPdg",         &mctrack_motherPdg);
    fTrackTree->Branch("mctrack_ancestorPdg",       &mctrack_ancestorPdg);
    fTrackTree->Branch("mctrack_id",                &mctrack_id);
    fTrackTree->Branch("mctrack_motherId",          &mctrack_motherId);
    fTrackTree->Branch("mctrack_ancestorId",        &mctrack_ancestorId);
    fTrackTree->Branch("mctrack_process",           &mctrack_process);
    fTrackTree->Branch("mctrack_motherProcess",     &mctrack_motherProcess);
    fTrackTree->Branch("mctrack_ancestorProcess",   &mctrack_ancestorProcess);
    fTrackTree->Branch("mctrack_startX",            &mctrack_startX);
    fTrackTree->Branch("mctrack_startY",            &mctrack_startY);
    fTrackTree->Branch("mctrack_startZ",            &mctrack_startZ);
    fTrackTree->Branch("mctrack_startT",            &mctrack_startT);
    fTrackTree->Branch("mctrack_motherStartX",      &mctrack_motherStartX);
    fTrackTree->Branch("mctrack_motherStartY",      &mctrack_motherStartY);
    fTrackTree->Branch("mctrack_motherStartZ",      &mctrack_motherStartZ);
    fTrackTree->Branch("mctrack_motherStartT",      &mctrack_motherStartT);
    fTrackTree->Branch("mctrack_ancestorStartX",    &mctrack_ancestorStartX);
    fTrackTree->Branch("mctrack_ancestorStartY",    &mctrack_ancestorStartY);
    fTrackTree->Branch("mctrack_ancestorStartZ",    &mctrack_ancestorStartZ);
    fTrackTree->Branch("mctrack_ancestorStartT",    &mctrack_ancestorStartT);
    fTrackTree->Branch("mctrack_endX",              &mctrack_endX);
    fTrackTree->Branch("mctrack_endY",              &mctrack_endY);
    fTrackTree->Branch("mctrack_endZ",              &mctrack_endZ);
    fTrackTree->Branch("mctrack_endT",              &mctrack_endT);
    fTrackTree->Branch("mctrack_motherEndX",        &mctrack_motherEndX);
    fTrackTree->Branch("mctrack_motherEndY",        &mctrack_motherEndY);
    fTrackTree->Branch("mctrack_motherEndZ",        &mctrack_motherEndZ);
    fTrackTree->Branch("mctrack_motherEndT",        &mctrack_motherEndT);
    fTrackTree->Branch("mctrack_ancestorEndX",      &mctrack_ancestorEndX);
    fTrackTree->Branch("mctrack_ancestorEndY",      &mctrack_ancestorEndY);
    fTrackTree->Branch("mctrack_ancestorEndZ",      &mctrack_ancestorEndZ);
    fTrackTree->Branch("mctrack_ancestorEndT",      &mctrack_ancestorEndT);
    fTrackTree->Branch("mctrack_startPX",           &mctrack_startPX);
    fTrackTree->Branch("mctrack_startPY",           &mctrack_startPY);
    fTrackTree->Branch("mctrack_startPZ",           &mctrack_startPZ);
    fTrackTree->Branch("mctrack_startE",            &mctrack_startE);
    fTrackTree->Branch("mctrack_motherStartPX",     &mctrack_motherStartPX);
    fTrackTree->Branch("mctrack_motherStartPY",     &mctrack_motherStartPY);
    fTrackTree->Branch("mctrack_motherStartPZ",     &mctrack_motherStartPZ);
    fTrackTree->Branch("mctrack_motherStartE",      &mctrack_motherStartE);
    fTrackTree->Branch("mctrack_ancestorStartPX",   &mctrack_ancestorStartPX);
    fTrackTree->Branch("mctrack_ancestorStartPY",   &mctrack_ancestorStartPY);
    fTrackTree->Branch("mctrack_ancestorStartPZ",   &mctrack_ancestorStartPZ);
    fTrackTree->Branch("mctrack_ancestorStartE",    &mctrack_ancestorStartE);
    fTrackTree->Branch("mctrack_endPX",             &mctrack_endPX);
    fTrackTree->Branch("mctrack_endPY",             &mctrack_endPY);
    fTrackTree->Branch("mctrack_endPZ",             &mctrack_endPZ);
    fTrackTree->Branch("mctrack_endE",              &mctrack_endE);
    fTrackTree->Branch("mctrack_motherEndPX",       &mctrack_motherEndPX);
    fTrackTree->Branch("mctrack_motherEndPY",       &mctrack_motherEndPY);
    fTrackTree->Branch("mctrack_motherEndPZ",       &mctrack_motherEndPZ);
    fTrackTree->Branch("mctrack_motherEndE",        &mctrack_motherEndE);
    fTrackTree->Branch("mctrack_ancestorEndPX",     &mctrack_ancestorEndPX);
    fTrackTree->Branch("mctrack_ancestorEndPY",     &mctrack_ancestorEndPY);
    fTrackTree->Branch("mctrack_ancestorEndPZ",     &mctrack_ancestorEndPZ);
    fTrackTree->Branch("mctrack_ancestorEndE",      &mctrack_ancestorEndE);

    // MC SHOWER
    fShowerTree = new TTree("MCShower","MC Shower Tree from Output of LArG4");
    fShowerTree->Branch("Nmcshower",                  &Nmcshower);
    fShowerTree->Branch("mcshower_pdg",               &mcshower_pdg);
    fShowerTree->Branch("mcshower_motherPdg",         &mcshower_motherPdg);
    fShowerTree->Branch("mcshower_ancestorPdg",       &mcshower_ancestorPdg);
    fShowerTree->Branch("mcshower_id",                &mcshower_id);
    fShowerTree->Branch("mcshower_motherId",          &mcshower_motherId);
    fShowerTree->Branch("mcshower_ancestorId",        &mcshower_ancestorId);
    fShowerTree->Branch("mcshower_process",           &mcshower_process);
    fShowerTree->Branch("mcshower_motherProcess",     &mcshower_motherProcess);
    fShowerTree->Branch("mcshower_ancestorProcess",   &mcshower_ancestorProcess);
    fShowerTree->Branch("mcshower_startX",            &mcshower_startX);
    fShowerTree->Branch("mcshower_startY",            &mcshower_startY);
    fShowerTree->Branch("mcshower_startZ",            &mcshower_startZ);
    fShowerTree->Branch("mcshower_startT",            &mcshower_startT);
    fShowerTree->Branch("mcshower_motherStartX",      &mcshower_motherStartX);
    fShowerTree->Branch("mcshower_motherStartY",      &mcshower_motherStartY);
    fShowerTree->Branch("mcshower_motherStartZ",      &mcshower_motherStartZ);
    fShowerTree->Branch("mcshower_motherStartT",      &mcshower_motherStartT);
    fShowerTree->Branch("mcshower_ancestorStartX",    &mcshower_ancestorStartX);
    fShowerTree->Branch("mcshower_ancestorStartY",    &mcshower_ancestorStartY);
    fShowerTree->Branch("mcshower_ancestorStartZ",    &mcshower_ancestorStartZ);
    fShowerTree->Branch("mcshower_ancestorStartT",    &mcshower_ancestorStartT);
    fShowerTree->Branch("mcshower_endX",              &mcshower_endX);
    fShowerTree->Branch("mcshower_endY",              &mcshower_endY);
    fShowerTree->Branch("mcshower_endZ",              &mcshower_endZ);
    fShowerTree->Branch("mcshower_endT",              &mcshower_endT);
    fShowerTree->Branch("mcshower_motherEndX",        &mcshower_motherEndX);
    fShowerTree->Branch("mcshower_motherEndY",        &mcshower_motherEndY);
    fShowerTree->Branch("mcshower_motherEndZ",        &mcshower_motherEndZ);
    fShowerTree->Branch("mcshower_motherEndT",        &mcshower_motherEndT);
    fShowerTree->Branch("mcshower_ancestorEndX",      &mcshower_ancestorEndX);
    fShowerTree->Branch("mcshower_ancestorEndY",      &mcshower_ancestorEndY);
    fShowerTree->Branch("mcshower_ancestorEndZ",      &mcshower_ancestorEndZ);
    fShowerTree->Branch("mcshower_ancestorEndT",      &mcshower_ancestorEndT);
    fShowerTree->Branch("mcshower_startPX",           &mcshower_startPX);
    fShowerTree->Branch("mcshower_startPY",           &mcshower_startPY);
    fShowerTree->Branch("mcshower_startPZ",           &mcshower_startPZ);
    fShowerTree->Branch("mcshower_startE",            &mcshower_startE);
    fShowerTree->Branch("mcshower_motherStartPX",     &mcshower_motherStartPX);
    fShowerTree->Branch("mcshower_motherStartPY",     &mcshower_motherStartPY);
    fShowerTree->Branch("mcshower_motherStartPZ",     &mcshower_motherStartPZ);
    fShowerTree->Branch("mcshower_motherStartE",      &mcshower_motherStartE);
    fShowerTree->Branch("mcshower_ancestorStartPX",   &mcshower_ancestorStartPX);
    fShowerTree->Branch("mcshower_ancestorStartPY",   &mcshower_ancestorStartPY);
    fShowerTree->Branch("mcshower_ancestorStartPZ",   &mcshower_ancestorStartPZ);
    fShowerTree->Branch("mcshower_ancestorStartE",    &mcshower_ancestorStartE);
    fShowerTree->Branch("mcshower_endPX",             &mcshower_endPX);
    fShowerTree->Branch("mcshower_endPY",             &mcshower_endPY);
    fShowerTree->Branch("mcshower_endPZ",             &mcshower_endPZ);
    fShowerTree->Branch("mcshower_endE",              &mcshower_endE);
    fShowerTree->Branch("mcshower_motherEndPX",       &mcshower_motherEndPX);
    fShowerTree->Branch("mcshower_motherEndPY",       &mcshower_motherEndPY);
    fShowerTree->Branch("mcshower_motherEndPZ",       &mcshower_motherEndPZ);
    fShowerTree->Branch("mcshower_motherEndE",        &mcshower_motherEndE);
    fShowerTree->Branch("mcshower_ancestorEndPX",     &mcshower_ancestorEndPX);
    fShowerTree->Branch("mcshower_ancestorEndPY",     &mcshower_ancestorEndPY);
    fShowerTree->Branch("mcshower_ancestorEndPZ",     &mcshower_ancestorEndPZ);
    fShowerTree->Branch("mcshower_ancestorEndE",      &mcshower_ancestorEndE);
    fShowerTree->Branch("mcshower_detProfileX",       &mcshower_detProfileX);
    fShowerTree->Branch("mcshower_detProfileY",       &mcshower_detProfileY);
    fShowerTree->Branch("mcshower_detProfileZ",       &mcshower_detProfileZ);
    fShowerTree->Branch("mcshower_detProfileT",       &mcshower_detProfileT);
    fShowerTree->Branch("mcshower_detProfilePX",      &mcshower_detProfilePX);
    fShowerTree->Branch("mcshower_detProfilePY",      &mcshower_detProfilePY);
    fShowerTree->Branch("mcshower_detProfilePZ",      &mcshower_detProfilePZ);
    fShowerTree->Branch("mcshower_detProfileE",       &mcshower_detProfileE);
    fShowerTree->Branch("mcshower_dEdx",              &mcshower_dEdx);

    // SIM CHANNEL
    //fSimChTree = new TTree("SimChannel","SimChannel Tree from Output of LArG4");
    fEventTree->Branch("simide_size",           &Nsimchannel);
    fEventTree->Branch("simide_channelIdY",     &simchannel_channelID);
    fEventTree->Branch("simide_tdc",            &simchannel_tdc);
    fEventTree->Branch("simide_nEnergyDeposits", &simchannel_nEnergyDeposits);
    fEventTree->Branch("simide_trackId",        &simchannel_id);
    fEventTree->Branch("simide_numElectrons",   &simchannel_nElectrons);
    fEventTree->Branch("simide_energy",         &simchannel_energy);
    fEventTree->Branch("simide_x",              &simchannel_x);
    fEventTree->Branch("simide_y",              &simchannel_y);
    fEventTree->Branch("simide_z",              &simchannel_z);

    // OPFLASH --- BEAM ---
    fBeamOpFlashTree = new TTree("BeamOpFlash","Reco Beam OpFlash Tree");
    fBeamOpFlashTree->Branch("ofBeam_nFlash", &ofBeam_nFlash);
    fBeamOpFlashTree->Branch("ofBeam_t", &ofBeam_t);                              
    fBeamOpFlashTree->Branch("ofBeam_peTotal", &ofBeam_peTotal);                     
    fBeamOpFlashTree->Branch("ofBeam_multiplicity", &ofBeam_multiplicity);            
    fBeamPEperOpDet = new TClonesArray("TH1F");
    fBeamOpFlashTree->Branch("Beam_pe_opdet", &fBeamPEperOpDet, 256000, 0);

    // OPFLASH --- COSMIC ---
    fCosmicOpFlashTree = new TTree("CosmicOpFlash","Reco Cosmic OpFlash Tree");
    fCosmicOpFlashTree->Branch("ofCosmic_nFlash", &ofCosmic_nFlash);
    fCosmicOpFlashTree->Branch("ofCosmic_t", &ofCosmic_t);                              
    fCosmicOpFlashTree->Branch("ofCosmic_peTotal", &ofCosmic_peTotal);                     
    fCosmicOpFlashTree->Branch("ofCosmic_multiplicity", &ofCosmic_multiplicity);            
    fCosmicPEperOpDet = new TClonesArray("TH1F");
    fCosmicOpFlashTree->Branch("Cosmic_pe_opdet", &fCosmicPEperOpDet, 256000, 0);

    // dummy variables not used
    fEventTree->Branch("calib_nChannel", &fCalib_nChannel);           
    fEventTree->Branch("calib_channelId" , &fCalib_channelId);        
    fCalib_wf = new TClonesArray("TH1F");
    fEventTree->Branch("calib_wf", &fCalib_wf, 256000, 0);

    gDirectory = tmpDir;
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::beginJob()
  {
    std::cout << "Begin job" << std::endl;
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::endJob()
  {
    std::cout << "End job" << std::endl;
    TDirectory* tmpDir = gDirectory;
    fOutFile->cd("/Event");
    fEventTree->Write();
    if(fSaveMCNeutrino == true) { fNeutrinoTree->Write(); }
    if(fSaveMCParticle == true) { fParticleTree->Write(); }
    if(fSaveMCTrack == true) { fTrackTree->Write(); }
    if(fSaveMCShower == true) { fShowerTree->Write(); }
    //if(fSaveSimChannel == true) { fSimChTree->Write(); }
    if(fBeamSaveOpFlash == true) { fBeamOpFlashTree->Write(); }
    if(fCosmicSaveOpFlash == true) { fCosmicOpFlashTree->Write(); }
    gDirectory = tmpDir;
    fOutFile->Close();
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::beginRun(const art::Run&)
  {
    std::cout << "Begin run" << std::endl;
    mf::LogInfo("CellTreeTruth") << "begin run";
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::analyze(const art::Event& event)
  {
    std::cout << "Analyze" << std::endl;
    reset();

    fRun = event.run();
    fSubRun = event.subRun();
    fEvent = event.id().event();
    art::Timestamp ts = event.time();
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    fEventTime = tts.AsDouble();
    processSimChannel(event); 
    fEventTree->Fill();

    if(fSaveMCNeutrino == true) { processMCNeutrino(event); fNeutrinoTree->Fill(); }
    if(fSaveMCParticle == true) { processMCParticle(event); fParticleTree->Fill(); }
    if(fSaveMCTrack == true) { processMCTrack(event); fTrackTree->Fill(); }
    if(fSaveMCShower == true) { processMCShower(event); fShowerTree->Fill(); }
    if(fBeamSaveOpFlash == true) { processBeamOpFlash(event); fBeamOpFlashTree->Fill(); }
    if(fCosmicSaveOpFlash == true) { processCosmicOpFlash(event); fCosmicOpFlashTree->Fill(); }
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::reset()
  {
    std::cout << "Reset" << std::endl;
    if(fSaveMCNeutrino == true){
      mcneutrino_nuStatus.clear();
      mcneutrino_nuTrackId.clear();
      mcneutrino_nuPdg.clear();
      mcneutrino_nuMother.clear();
      mcneutrino_nuProcess.clear();
      mcneutrino_nuEndProcess.clear();
      mcneutrino_nuMass.clear();
      mcneutrino_nuWeight.clear();
      mcneutrino_nuVertexX.clear();
      mcneutrino_nuVertexY.clear();
      mcneutrino_nuVertexZ.clear();
      mcneutrino_nuEnergy.clear(); 
      //mcneutrino_nuGvtxX.clear();
      //mcneutrino_nuGvtxY.clear();
      //mcneutrino_nuGvtxZ.clear();
      mcneutrino_nuRescatter.clear();
      mcneutrino_leptonStatus.clear();
      mcneutrino_leptonTrackId.clear();
      mcneutrino_leptonPdg.clear();
      mcneutrino_leptonMother.clear();
      mcneutrino_leptonProcess.clear();
      mcneutrino_leptonEndProcess.clear();
      mcneutrino_leptonMass.clear();
      mcneutrino_leptonWeight.clear();
      mcneutrino_leptonVertexX.clear();
      mcneutrino_leptonVertexY.clear();
      mcneutrino_leptonVertexZ.clear();
      mcneutrino_leptonEnergy.clear();       
      //mcneutrino_leptonGvtxX.clear();
      //mcneutrino_leptonGvtxY.clear();
      //mcneutrino_leptonGvtxZ.clear();
      mcneutrino_leptonRescatter.clear();
      mcneutrino_ccnc.clear();
      mcneutrino_mode.clear();
      mcneutrino_interactionType.clear();
      mcneutrino_target.clear();
      mcneutrino_hitNuc.clear();
      mcneutrino_hitQuark.clear();
      mcneutrino_w.clear();
      mcneutrino_x.clear();
      mcneutrino_y.clear();
      mcneutrino_qSqr.clear();
      mcneutrino_pt.clear();
      mcneutrino_theta.clear();
    }

    if(fSaveMCParticle == true){
      Nmcparticle = 0;    
      mcparticle_id.clear();
      mcparticle_statusCode.clear();
      mcparticle_pdg.clear();
      mcparticle_mother.clear();
      mcparticle_process.clear();
      mcparticle_endProcess.clear();
      mcparticle_ndaughters.clear();
      mcparticle_daughterId.clear();
      mcparticle_ntrajpts.clear();
      mcparticle_x.clear();
      mcparticle_y.clear();
      mcparticle_z.clear();
      mcparticle_t.clear();
      mcparticle_px.clear();
      mcparticle_py.clear();
      mcparticle_pz.clear();
      mcparticle_e.clear();
      mcparticle_p.clear();
      mcparticle_pt.clear();
      mcparticle_endX.clear();
      mcparticle_endY.clear();
      mcparticle_endZ.clear();
      mcparticle_endT.clear();
      mcparticle_mass.clear();
      mcparticle_endPX.clear();
      mcparticle_endPY.clear();
      mcparticle_endPZ.clear();
      mcparticle_endE.clear();
      mcparticle_GvtxX.clear();
      mcparticle_GvtxY.clear();
      mcparticle_GvtxZ.clear();
      mcparticle_GvtxPX.clear();
      mcparticle_GvtxPY.clear();
      mcparticle_GvtxPZ.clear();
      mcparticle_GvtxE.clear();
      mcparticle_rescatter.clear();
      mcparticle_weight.clear();
    }

    if(fSaveMCTrack == true){
    mctrack_pdg.clear();
    mctrack_motherPdg.clear();
    mctrack_ancestorPdg.clear();
    mctrack_id.clear();
    mctrack_motherId.clear();
    mctrack_ancestorId.clear();
    mctrack_process.clear();
    mctrack_motherProcess.clear();
    mctrack_ancestorProcess.clear();
    mctrack_startX.clear();
    mctrack_startY.clear();
    mctrack_startZ.clear();
    mctrack_startT.clear();
    mctrack_motherStartX.clear();
    mctrack_motherStartY.clear();
    mctrack_motherStartZ.clear();
    mctrack_motherStartT.clear();
    mctrack_ancestorStartX.clear();
    mctrack_ancestorStartY.clear();
    mctrack_ancestorStartZ.clear();
    mctrack_ancestorStartT.clear();
    mctrack_endX.clear();
    mctrack_endY.clear();
    mctrack_endZ.clear();
    mctrack_endT.clear();
    mctrack_motherEndX.clear();
    mctrack_motherEndY.clear();
    mctrack_motherEndZ.clear();
    mctrack_motherEndT.clear();
    mctrack_ancestorEndX.clear();
    mctrack_ancestorEndY.clear();
    mctrack_ancestorEndZ.clear();
    mctrack_ancestorEndT.clear();
    mctrack_startPX.clear();
    mctrack_startPY.clear();
    mctrack_startPZ.clear();
    mctrack_startE.clear();
    mctrack_motherStartPX.clear();
    mctrack_motherStartPY.clear();
    mctrack_motherStartPZ.clear();
    mctrack_motherStartE.clear();
    mctrack_ancestorStartPX.clear();
    mctrack_ancestorStartPY.clear();
    mctrack_ancestorStartPZ.clear();
    mctrack_ancestorStartE.clear();
    mctrack_endPX.clear();
    mctrack_endPY.clear();
    mctrack_endPZ.clear();
    mctrack_endE.clear();
    mctrack_motherEndPX.clear();
    mctrack_motherEndPY.clear();
    mctrack_motherEndPZ.clear();
    mctrack_motherEndE.clear();
    mctrack_ancestorEndPX.clear();
    mctrack_ancestorEndPY.clear();
    mctrack_ancestorEndPZ.clear();
    mctrack_ancestorEndE.clear();
    }

    if(fSaveMCShower == true){
    mcshower_pdg.clear(); 
    mcshower_motherPdg.clear(); 
    mcshower_ancestorPdg.clear(); 
    mcshower_id.clear();
    mcshower_motherId.clear();
    mcshower_ancestorId.clear();
    mcshower_process.clear();
    mcshower_motherProcess.clear();
    mcshower_ancestorProcess.clear();
    mcshower_startX.clear();
    mcshower_startY.clear();
    mcshower_startZ.clear();
    mcshower_startT.clear();
    mcshower_motherStartX.clear();
    mcshower_motherStartY.clear();
    mcshower_motherStartZ.clear();
    mcshower_motherStartT.clear();
    mcshower_ancestorStartX.clear();
    mcshower_ancestorStartY.clear();
    mcshower_ancestorStartZ.clear();
    mcshower_ancestorStartT.clear();
    mcshower_endX.clear();
    mcshower_endY.clear();
    mcshower_endZ.clear();
    mcshower_endT.clear();
    mcshower_motherEndX.clear();
    mcshower_motherEndY.clear();
    mcshower_motherEndZ.clear();
    mcshower_motherEndT.clear();
    mcshower_ancestorEndX.clear();
    mcshower_ancestorEndY.clear();
    mcshower_ancestorEndZ.clear();
    mcshower_ancestorEndT.clear();
    mcshower_startPX.clear();
    mcshower_startPY.clear();
    mcshower_startPZ.clear();
    mcshower_startE.clear();
    mcshower_motherStartPX.clear();
    mcshower_motherStartPY.clear();
    mcshower_motherStartPZ.clear();
    mcshower_motherStartE.clear();
    mcshower_ancestorStartPX.clear();
    mcshower_ancestorStartPY.clear();
    mcshower_ancestorStartPZ.clear();
    mcshower_ancestorStartE.clear();
    mcshower_endPX.clear();
    mcshower_endPY.clear();
    mcshower_endPZ.clear();
    mcshower_endE.clear();
    mcshower_motherEndPX.clear();
    mcshower_motherEndPY.clear();
    mcshower_motherEndPZ.clear();
    mcshower_motherEndE.clear();
    mcshower_ancestorEndPX.clear();
    mcshower_ancestorEndPY.clear();
    mcshower_ancestorEndPZ.clear();
    mcshower_ancestorEndE.clear();
    mcshower_detProfileX.clear();
    mcshower_detProfileY.clear();
    mcshower_detProfileZ.clear();
    mcshower_detProfileT.clear();
    mcshower_detProfilePX.clear();
    mcshower_detProfilePY.clear();
    mcshower_detProfilePZ.clear();
    mcshower_detProfileE.clear();
    mcshower_dEdx.clear();
    }

    //if(fSaveSimChannel == true){
    Nsimchannel = 0;
    simchannel_channelID.clear();
    simchannel_tdc.clear();
    simchannel_nEnergyDeposits.clear();
    simchannel_id.clear();
    simchannel_nElectrons.clear();
    simchannel_energy.clear();
    simchannel_x.clear();
    simchannel_y.clear();
    simchannel_z.clear();
    //}

    if(fBeamSaveOpFlash == true){
      ofBeam_t.clear();
      ofBeam_peTotal.clear();
      ofBeam_multiplicity.clear();
      fBeamPEperOpDet->Delete();
    }

    if(fCosmicSaveOpFlash == true){
      ofCosmic_t.clear();
      ofCosmic_peTotal.clear();
      ofCosmic_multiplicity.clear();
      fCosmicPEperOpDet->Delete();
    }

    // dummy variables not used
    fCalib_channelId.clear();
    fCalib_wf->Clear();
  }

  //-------------------------------------------------------------------

  //-------------------------------------------------------------------
  void CellTreeTruth::processMCNeutrino(const art::Event& event)
  {
    std::cout << "Process MCNeutrino" << std::endl;
    art::Handle< std::vector<simb::MCTruth> > mctruthHandle;
    event.getByLabel(fMCNeutrinoLabel, mctruthHandle);
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    art::fill_ptr_vector(mclist, mctruthHandle);
    art::Ptr<simb::MCTruth> mctruth;

    if(mclist.size() > 0){
      mctruth = mclist.at(0);
      simb::MCNeutrino nu = mctruth->GetNeutrino();
      //### neutrino ###                                               
      mcneutrino_nuStatus.push_back(nu.Nu().StatusCode());
      mcneutrino_nuTrackId.push_back(nu.Nu().TrackId());
      mcneutrino_nuPdg.push_back(nu.Nu().PdgCode());
      mcneutrino_nuMother.push_back(nu.Nu().Mother());
      mcneutrino_nuProcess.push_back(nu.Nu().Process());
      mcneutrino_nuEndProcess.push_back(nu.Nu().EndProcess());
      mcneutrino_nuMass.push_back(nu.Nu().Mass());
      mcneutrino_nuWeight.push_back(nu.Nu().Weight());
      mcneutrino_nuVertexX.push_back(nu.Nu().Vx());
      mcneutrino_nuVertexY.push_back(nu.Nu().Vy());
      mcneutrino_nuVertexZ.push_back(nu.Nu().Vz());
      mcneutrino_nuEnergy.push_back(nu.Nu().E()); 
      //mcneutrino_nuGvtxX.push_back(nu.Nu().GetGvtx().X());
      //mcneutrino_nuGvtxY.push_back(nu.Nu().GetGvtx().Y());
      //mcneutrino_nuGvtxZ.push_back(nu.Nu().GetGvtx().Z());
      mcneutrino_nuRescatter.push_back(nu.Nu().Rescatter());
      //### lepton ###                                                
      mcneutrino_leptonStatus.push_back(nu.Lepton().StatusCode());
      mcneutrino_leptonTrackId.push_back(nu.Lepton().TrackId());
      mcneutrino_leptonPdg.push_back(nu.Lepton().PdgCode());
      mcneutrino_leptonMother.push_back(nu.Lepton().Mother());
      mcneutrino_leptonProcess.push_back(nu.Lepton().Process());
      mcneutrino_leptonEndProcess.push_back(nu.Lepton().EndProcess());
      mcneutrino_leptonMass.push_back(nu.Lepton().Mass());
      mcneutrino_leptonWeight.push_back(nu.Lepton().Weight());
      mcneutrino_leptonVertexX.push_back(nu.Lepton().Vx());
      mcneutrino_leptonVertexY.push_back(nu.Lepton().Vy());
      mcneutrino_leptonVertexZ.push_back(nu.Lepton().Vz());
      mcneutrino_leptonEnergy.push_back(nu.Lepton().E()); 
      //mcneutrino_leptonGvtxX.push_back(nu.Lepton().GetGvtx().X());
      //mcneutrino_leptonGvtxY.push_back(nu.Lepton().GetGvtx().Y());
      //mcneutrino_leptonGvtxZ.push_back(nu.Lepton().GetGvtx().Z());
      mcneutrino_leptonRescatter.push_back(nu.Lepton().Rescatter());
      //### other ###                                   
      mcneutrino_mode.push_back(nu.Mode());
      mcneutrino_interactionType.push_back(nu.InteractionType());
      mcneutrino_target.push_back(nu.Target());
      mcneutrino_hitNuc.push_back(nu.HitNuc());
      mcneutrino_hitQuark.push_back(nu.HitQuark());
      mcneutrino_w.push_back(nu.W());
      mcneutrino_x.push_back(nu.X());
      mcneutrino_y.push_back(nu.Y());
      mcneutrino_qSqr.push_back(nu.QSqr());
      mcneutrino_pt.push_back(nu.Pt());
      mcneutrino_theta.push_back(nu.Theta());
    }
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::processMCParticle(const art::Event& event)
  {
    std::cout << "Process MCParticle" << std::endl;
    art::Handle<std::vector<simb::MCParticle> > particleHandle;
    if(! event.getByLabel(fMCParticleLabel, particleHandle)){
      cout << "WARNING: no label " << fMCParticleLabel << endl;
      return;
    }
   
    std::vector<art::Ptr<simb::MCParticle> > particles;
    art::fill_ptr_vector(particles, particleHandle);
   
    int i=0; // track index in saved MCParticles
    for(auto const& particle : particles){
      mcparticle_id.push_back(particle->TrackId());
      mcparticle_statusCode.push_back(particle->StatusCode());
      mcparticle_pdg.push_back(particle->PdgCode());
      mcparticle_mother.push_back(particle->Mother());
      mcparticle_process.push_back(particle->Process());
      mcparticle_endProcess.push_back(particle->EndProcess());
      mcparticle_ndaughters.push_back(particle->NumberDaughters());
      for(int a=0; a<particle->NumberDaughters(); a++){
	mcparticle_daughterId.push_back(particle->Daughter(a));
      }
      mcparticle_ntrajpts.push_back(particle->NumberTrajectoryPoints());
      for(int b=0; b<(int)particle->NumberTrajectoryPoints(); b++){
	mcparticle_x.push_back(particle->Vx(b));
	mcparticle_y.push_back(particle->Vy(b));
	mcparticle_z.push_back(particle->Vz(b));
	mcparticle_t.push_back(particle->T(b));
	mcparticle_px.push_back(particle->Px(b));
	mcparticle_py.push_back(particle->Py(b));
	mcparticle_pz.push_back(particle->Pz(b));
	mcparticle_e.push_back(particle->E(b));
	mcparticle_p.push_back(particle->P(b));
	mcparticle_pt.push_back(particle->Pt(b));
      }
      mcparticle_endX.push_back(particle->EndX());
      mcparticle_endY.push_back(particle->EndY());
      mcparticle_endZ.push_back(particle->EndZ());
      mcparticle_endT.push_back(particle->EndT());
      mcparticle_mass.push_back(particle->Mass());
      mcparticle_endPX.push_back(particle->EndPx());
      mcparticle_endPY.push_back(particle->EndPy());
      mcparticle_endPZ.push_back(particle->EndPz());
      mcparticle_endE.push_back(particle->EndE()); 
      mcparticle_GvtxX.push_back(particle->GetGvtx().X());
      mcparticle_GvtxY.push_back(particle->GetGvtx().Y());
      mcparticle_GvtxZ.push_back(particle->GetGvtx().Z());
      mcparticle_GvtxPX.push_back(particle->GetGvtx().Px());
      mcparticle_GvtxPY.push_back(particle->GetGvtx().Py());
      mcparticle_GvtxPZ.push_back(particle->GetGvtx().Pz());
      mcparticle_GvtxE.push_back(particle->GetGvtx().Energy());
      mcparticle_rescatter.push_back(particle->Rescatter());
      mcparticle_weight.push_back(particle->Weight());
    } 
    Nmcparticle = i;
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::processMCTrack(const art::Event& event)
  {
    std::cout << "Process MCTrack" << std::endl;
    art::Handle< std::vector<sim::MCTrack> > trackHandle;
    if(! event.getByLabel(fMCTrackLabel, trackHandle)){
      cout << "WARNING: no label " << fMCTrackLabel << endl;
      return;
    }
    std::vector<art::Ptr<sim::MCTrack> > tracks;
    art::fill_ptr_vector(tracks, trackHandle);
    
    Nmctrack = (int)tracks.size();
    for(auto const& track : tracks){
      // origin
      mctrack_pdg.push_back(track->PdgCode());
      mctrack_id.push_back(track->TrackID());
      mctrack_process.push_back(track->Process());
      mctrack_motherPdg.push_back(track->MotherPdgCode());
      mctrack_motherId.push_back(track->MotherTrackID());
      mctrack_motherProcess.push_back(track->MotherProcess());
      mctrack_ancestorPdg.push_back(track->AncestorPdgCode());
      mctrack_ancestorId.push_back(track->AncestorTrackID());
      mctrack_ancestorProcess.push_back(track->AncestorProcess());
      mctrack_startX.push_back(track->Start().Position().X());
      mctrack_startY.push_back(track->Start().Position().Y());
      mctrack_startZ.push_back(track->Start().Position().Z());
      mctrack_startT.push_back(track->Start().Position().T());
      mctrack_motherStartX.push_back(track->MotherStart().Position().X());
      mctrack_motherStartY.push_back(track->MotherStart().Position().Y());
      mctrack_motherStartZ.push_back(track->MotherStart().Position().Z());
      mctrack_motherStartT.push_back(track->MotherStart().Position().T());
      mctrack_ancestorStartX.push_back(track->AncestorStart().Position().X());
      mctrack_ancestorStartY.push_back(track->AncestorStart().Position().Y());
      mctrack_ancestorStartZ.push_back(track->AncestorStart().Position().Z());
      mctrack_ancestorStartT.push_back(track->AncestorStart().Position().T());
      mctrack_endX.push_back(track->End().Position().X());
      mctrack_endY.push_back(track->End().Position().Y());
      mctrack_endZ.push_back(track->End().Position().Z());
      mctrack_endT.push_back(track->End().Position().T());
      mctrack_motherEndX.push_back(track->MotherEnd().Position().X());
      mctrack_motherEndY.push_back(track->MotherEnd().Position().Y());
      mctrack_motherEndZ.push_back(track->MotherEnd().Position().Z());
      mctrack_motherEndT.push_back(track->MotherEnd().Position().T());
      mctrack_ancestorEndX.push_back(track->AncestorEnd().Position().X());
      mctrack_ancestorEndY.push_back(track->AncestorEnd().Position().Y());
      mctrack_ancestorEndZ.push_back(track->AncestorEnd().Position().Z());
      mctrack_ancestorEndT.push_back(track->AncestorEnd().Position().T());
      mctrack_startPX.push_back(track->Start().Momentum().Px());
      mctrack_startPY.push_back(track->Start().Momentum().Py());
      mctrack_startPZ.push_back(track->Start().Momentum().Pz());
      mctrack_startE.push_back(track->Start().Momentum().E());
      mctrack_motherStartPX.push_back(track->MotherStart().Momentum().Px());
      mctrack_motherStartPY.push_back(track->MotherStart().Momentum().Py());
      mctrack_motherStartPZ.push_back(track->MotherStart().Momentum().Pz());
      mctrack_motherStartE.push_back(track->MotherStart().Momentum().E());
      mctrack_ancestorStartPX.push_back(track->AncestorStart().Momentum().Px());
      mctrack_ancestorStartPY.push_back(track->AncestorStart().Momentum().Py());
      mctrack_ancestorStartPZ.push_back(track->AncestorStart().Momentum().Pz());
      mctrack_ancestorStartE.push_back(track->AncestorStart().Momentum().E());
      mctrack_endPX.push_back(track->End().Momentum().Px());
      mctrack_endPY.push_back(track->End().Momentum().Py());
      mctrack_endPZ.push_back(track->End().Momentum().Pz());
      mctrack_endE.push_back(track->End().Momentum().E());
      mctrack_motherEndPX.push_back(track->MotherEnd().Momentum().Px());
      mctrack_motherEndPY.push_back(track->MotherEnd().Momentum().Py());
      mctrack_motherEndPZ.push_back(track->MotherEnd().Momentum().Pz());
      mctrack_motherEndE.push_back(track->MotherEnd().Momentum().E());
      mctrack_ancestorEndPX.push_back(track->AncestorEnd().Momentum().Px());
      mctrack_ancestorEndPY.push_back(track->AncestorEnd().Momentum().Py());
      mctrack_ancestorEndPZ.push_back(track->AncestorEnd().Momentum().Pz());
      mctrack_ancestorEndE.push_back(track->AncestorEnd().Momentum().E());      
    }
  }
  
  //-------------------------------------------------------------------
  void CellTreeTruth::processMCShower(const art::Event& event)
  {
    std::cout << "Process MCShower" << std::endl;
    art::Handle< std::vector<sim::MCShower> > showerHandle;
    if(! event.getByLabel(fMCShowerLabel, showerHandle)){
      cout << "WARNING: no label " << fMCShowerLabel << endl;
      return;
    }
    std::vector<art::Ptr<sim::MCShower> > showers;
    art::fill_ptr_vector(showers, showerHandle);
    
    Nmcshower = (int)showers.size();
    for(auto const& shower : showers){
      mcshower_pdg.push_back(shower->PdgCode());
      mcshower_id.push_back(shower->TrackID());
      mcshower_process.push_back(shower->Process());
      mcshower_motherPdg.push_back(shower->MotherPdgCode());
      mcshower_motherId.push_back(shower->MotherTrackID());
      mcshower_motherProcess.push_back(shower->MotherProcess());
      mcshower_ancestorPdg.push_back(shower->AncestorPdgCode());
      mcshower_ancestorId.push_back(shower->AncestorTrackID());
      mcshower_ancestorProcess.push_back(shower->AncestorProcess());
      mcshower_startX.push_back(shower->Start().Position().X());
      mcshower_startY.push_back(shower->Start().Position().Y());
      mcshower_startZ.push_back(shower->Start().Position().Z());
      mcshower_startT.push_back(shower->Start().Position().T());
      mcshower_motherStartX.push_back(shower->MotherStart().Position().X());
      mcshower_motherStartY.push_back(shower->MotherStart().Position().Y());
      mcshower_motherStartZ.push_back(shower->MotherStart().Position().Z());
      mcshower_motherStartT.push_back(shower->MotherStart().Position().T());
      mcshower_ancestorStartX.push_back(shower->AncestorStart().Position().X());
      mcshower_ancestorStartY.push_back(shower->AncestorStart().Position().Y());
      mcshower_ancestorStartZ.push_back(shower->AncestorStart().Position().Z());
      mcshower_ancestorStartT.push_back(shower->AncestorStart().Position().T());
      mcshower_endX.push_back(shower->End().Position().X());
      mcshower_endY.push_back(shower->End().Position().Y());
      mcshower_endZ.push_back(shower->End().Position().Z());
      mcshower_endT.push_back(shower->End().Position().T());
      mcshower_motherEndX.push_back(shower->MotherEnd().Position().X());
      mcshower_motherEndY.push_back(shower->MotherEnd().Position().Y());
      mcshower_motherEndZ.push_back(shower->MotherEnd().Position().Z());
      mcshower_motherEndT.push_back(shower->MotherEnd().Position().T());
      mcshower_ancestorEndX.push_back(shower->AncestorEnd().Position().X());
      mcshower_ancestorEndY.push_back(shower->AncestorEnd().Position().Y());
      mcshower_ancestorEndZ.push_back(shower->AncestorEnd().Position().Z());
      mcshower_ancestorEndT.push_back(shower->AncestorEnd().Position().T());
      mcshower_startPX.push_back(shower->Start().Momentum().Px());
      mcshower_startPY.push_back(shower->Start().Momentum().Py());
      mcshower_startPZ.push_back(shower->Start().Momentum().Pz());
      mcshower_startE.push_back(shower->Start().Momentum().E());
      mcshower_motherStartPX.push_back(shower->MotherStart().Momentum().Px());
      mcshower_motherStartPY.push_back(shower->MotherStart().Momentum().Py());
      mcshower_motherStartPZ.push_back(shower->MotherStart().Momentum().Pz());
      mcshower_motherStartE.push_back(shower->MotherStart().Momentum().E());
      mcshower_ancestorStartPX.push_back(shower->AncestorStart().Momentum().Px());
      mcshower_ancestorStartPY.push_back(shower->AncestorStart().Momentum().Py());
      mcshower_ancestorStartPZ.push_back(shower->AncestorStart().Momentum().Pz());
      mcshower_ancestorStartE.push_back(shower->AncestorStart().Momentum().E());
      mcshower_endPX.push_back(shower->End().Momentum().Px());
      mcshower_endPY.push_back(shower->End().Momentum().Py());
      mcshower_endPZ.push_back(shower->End().Momentum().Pz());
      mcshower_endE.push_back(shower->End().Momentum().E());
      mcshower_motherEndPX.push_back(shower->MotherEnd().Momentum().Px());
      mcshower_motherEndPY.push_back(shower->MotherEnd().Momentum().Py());
      mcshower_motherEndPZ.push_back(shower->MotherEnd().Momentum().Pz());
      mcshower_motherEndE.push_back(shower->MotherEnd().Momentum().E());
      mcshower_ancestorEndPX.push_back(shower->AncestorEnd().Momentum().Px());
      mcshower_ancestorEndPY.push_back(shower->AncestorEnd().Momentum().Py());
      mcshower_ancestorEndPZ.push_back(shower->AncestorEnd().Momentum().Pz());
      mcshower_ancestorEndE.push_back(shower->AncestorEnd().Momentum().E());
      mcshower_detProfileX.push_back(shower->DetProfile().Position().X());
      mcshower_detProfileY.push_back(shower->DetProfile().Position().Y());
      mcshower_detProfileZ.push_back(shower->DetProfile().Position().Z());
      mcshower_detProfileT.push_back(shower->DetProfile().Position().T());
      mcshower_detProfilePX.push_back(shower->DetProfile().Momentum().Px());
      mcshower_detProfilePY.push_back(shower->DetProfile().Momentum().Py());
      mcshower_detProfilePZ.push_back(shower->DetProfile().Momentum().Pz());
      mcshower_detProfileE.push_back(shower->DetProfile().Momentum().E());
      mcshower_dEdx.push_back(shower->dEdx());
    }
  }  

  //---------------------------------------------------------------------
  void CellTreeTruth::processSimChannel(const art::Event& event)
  {
    std::cout << "Process SimChannel" << std::endl;
    art::Handle< std::vector<sim::SimChannel> > simchannelHandle;
    if(! event.getByLabel(fSimChannelLabel, simchannelHandle)){
      cout << "WARNING: no label " << fSimChannelLabel << endl;
      return;
    }
    std::vector<art::Ptr<sim::SimChannel> > simch;
    art::fill_ptr_vector(simch, simchannelHandle);

    Nsimchannel = 0;
    //for(auto const& sc : simch){
    for(auto const& sc : (*simchannelHandle)){
      auto channelNumber = sc.Channel();
      auto const& timeSlices = sc.TDCIDEMap();
      for(auto const& timeSlice : timeSlices){
	auto const& energyDeposits = timeSlice.second;
	simchannel_nEnergyDeposits.push_back((int)energyDeposits.size());
	for(auto const& energyDeposit : energyDeposits){
	  Nsimchannel++;
	  simchannel_channelID.push_back(channelNumber);
	  simchannel_tdc.push_back(timeSlice.first);
	  simchannel_id.push_back((int)energyDeposit.trackID);
	  simchannel_nElectrons.push_back(energyDeposit.numElectrons);
	  simchannel_energy.push_back(energyDeposit.energy);
	  simchannel_x.push_back(energyDeposit.x);
	  simchannel_y.push_back(energyDeposit.y);
	  simchannel_z.push_back(energyDeposit.z); 
	}
      }
    }
  }

  //----------------------------------------------------------------------          
  void CellTreeTruth::processBeamOpFlash( const art::Event& event)
  {
    std::cout << "Process Beam OpFlash" << std::endl;
    art::Handle<std::vector<recob::OpFlash> > flash_handle;
    if(! event.getByLabel(fBeamOpFlashLabel, flash_handle)){
      cout << "WARNING: no label " << fBeamOpFlashLabel << endl;
      return;
    }
    std::vector<art::Ptr<recob::OpFlash> > flashes;
    art::fill_ptr_vector(flashes, flash_handle);
    ofBeam_nFlash = (int)flashes.size();

    int a=0;
    int nOpDet = fGeometry->NOpDets();

    for(auto const& flash: flashes){
      ofBeam_t.push_back(flash->Time());
      ofBeam_peTotal.push_back(flash->TotalPE());
      TH1F *h = new ((*fBeamPEperOpDet)[a]) TH1F("","",nOpDet,0,nOpDet);

      int mult = 0;
      for(int i=0; i<nOpDet; ++i){
        if(flash->PE(i) >= opBeamMultPEThresh){
          mult++;
        }
        h->SetBinContent(i, flash->PE(i));
      }
      ofBeam_multiplicity.push_back(mult);
      a++;
    }
  }

  //----------------------------------------------------------------------
void CellTreeTruth::processCosmicOpFlash( const art::Event& event)
  {
    std::cout << "Process Cosmic OpFlash" << std::endl;
    art::Handle<std::vector<recob::OpFlash> > flash_handle;
    if(! event.getByLabel(fCosmicOpFlashLabel, flash_handle)){
      cout << "WARNING: no label " << fCosmicOpFlashLabel << endl;
      return;
    }
    std::vector<art::Ptr<recob::OpFlash> > flashes;
    art::fill_ptr_vector(flashes, flash_handle);
    ofCosmic_nFlash = (int)flashes.size();

    int a=0;
    int nOpDet = fGeometry->NOpDets();

    for(auto const& flash: flashes){
      ofCosmic_t.push_back(flash->Time());
      ofCosmic_peTotal.push_back(flash->TotalPE());
      TH1F *h = new ((*fCosmicPEperOpDet)[a]) TH1F("","",nOpDet,0,nOpDet);

      int mult = 0;
      for(int i=0; i<nOpDet; ++i){
        if(flash->PE(i) >= opCosmicMultPEThresh){
          mult++;
        }
        h->SetBinContent(i, flash->PE(i));
      }
      ofCosmic_multiplicity.push_back(mult);
      a++;
    }
  }


  DEFINE_ART_MODULE(CellTreeTruth)
} // namespace wc
#endif
