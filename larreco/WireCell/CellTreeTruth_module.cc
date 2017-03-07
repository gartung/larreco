#ifndef CELLTREETRUTH_MODULE
#define CELLTREETRUTH_MODULE

#ifdef __CINT__
#pramga link C++ class std::vector<TVector3>+;
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

// LArSoft includes
#include "lardataobj/MCBase/MCStep.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"

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

#define MAX_TRACKS 30000

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

    void processMCTruth(const art::Event& evt);
    void processMCNeutrino(const art::Event& evt);
    void processMCParticle(const art::Event& evt);
    void processMCTrack(const art::Event& evt);
    void processMCShower(const art::Event& evt);
    void processSimChannel(const art::Event& evt);

  private:
    bool fSaveMCTruth;
    bool fSaveMCNeutrino;
    bool fSaveMCParticle;
    bool fSaveMCTrack;
    bool fSaveMCShower;
    bool fSaveSimChannel;
    std::string fMCTruthLabel;
    std::string fMCNeutrinoLabel;
    std::string fMCParticleLabel;
    std::string fMCTrackLabel;
    std::string fMCShowerLabel;
    std::string fSimChannelLabel;
    std::string fOutFileName;

    TFile *fOutFile;
    TTree *fEventTree;
    TTree *fTrueTree;
    TTree *fNeutrinoTree;
    TTree *fParticleTree;
    TTree *fTrackTree;
    TTree *fShowerTree;
    TTree *fSimChTree;

    int fRun;
    int fSubRun;
    int fEvent;
    double fEventTime;

    // MC TRUTH
    int Nmctruth;

    // MC Neutrino;
    int Nmcneutrino;
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
    int Nmcparticle;
    vector<int> mcparticle_id;
    vector<int> mcparticle_statusCode;
    vector<int> mcparticle_pdg;
    vector<int> mcparticle_mother;
    TObjArray *fmcparticle_polarization;
    vector<std::string> mcparticle_process;
    vector<std::string> mcparticle_endProcess;
    vector<int> mcparticle_ndaughters;
    vector<int> mcparticle_daughterId;
    vector<int> mcparticle_ntrajpts;
    float mcparticle_startXYZT[MAX_TRACKS][4];
    float mcparticle_endXYZT[MAX_TRACKS][4];
    float mcparticle_startMomentum[MAX_TRACKS][4];
    float mcparticle_endMomentum[MAX_TRACKS][4];
    TObjArray *fmcparticle_position;
    vector<double> mcparticle_vx;                                      
    vector<double> mcparticle_vy;                              
    vector<double> mcparticle_vz;                              
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
    TObjArray *fmcparticle_Gvtx;                  
    vector<double> mcparticle_Gvx;
    vector<double> mcparticle_Gvy;
    vector<double> mcparticle_Gvz;
    vector<double> mcparticle_Gvt;
    vector<int> mcparticle_firstDaughter;
    vector<int> mcparticle_lastDaughter;
    vector<int> mcparticle_rescatter;
    // trajectory                                              
    vector<double> mcparticle_weight;


    // MC TRACK
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
    TObjArray *fmctrack_startPosition;
    TObjArray *fmctrack_motherStartPosition;
    TObjArray *fmctrack_ancestorStartPosition;
    TObjArray *fmctrack_endPosition;
    TObjArray *fmctrack_motherEndPosition;
    TObjArray *fmctrack_ancestorEndPosition;
    TObjArray *fmctrack_startMomentum;
    TObjArray *fmctrack_motherStartMomentum;
    TObjArray *fmctrack_ancestorStartMomentum;
    TObjArray *fmctrack_endMomentum;
    TObjArray *fmctrack_motherEndMomentum;
    TObjArray *fmctrack_ancestorEndMomentum;
    vector<vector<vector<double> > >mctrack_dQdx;
    vector<vector<double> > mctrack_dEdx;

    // MC SHOWER
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
    TObjArray *fmcshower_start;
    TObjArray *fmcshower_motherStart;
    TObjArray *fmcshower_ancestorStart;
    TObjArray *fmcshower_end;
    TObjArray *fmcshower_motherEnd;
    TObjArray *fmcshower_ancestorEnd;
    TObjArray *fmcshower_detprofilePos;
    TObjArray *fmcshower_detprofileMom;
    vector<vector<unsigned int> > mcshower_daughterTrackID;
    vector<vector<double> > mcshower_charge;
    vector<vector<double> > mcshower_dQdx;
    vector<double> mcshower_dEdx;
    TObjArray *fmcshower_startDir;

    // SIM CHANNEL
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
  }; // class

  //-------------------------------------------------------------------
  CellTreeTruth::CellTreeTruth(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
  {
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
    fMCTruthLabel = p.get<std::string>("MCTruthLabel");
    fMCNeutrinoLabel = p.get<std::string>("MCNeutrinoLabel");
    fMCParticleLabel = p.get<std::string>("MCParticleLabel");
    fMCTrackLabel = p.get<std::string>("MCTrackLabel");
    fMCShowerLabel = p.get<std::string>("MCShowerLabel");
    fSimChannelLabel = p.get<std::string>("SimChannelLabel");

    fSaveMCTruth = p.get<bool>("saveMCTruth");
    fSaveMCNeutrino = p.get<bool>("saveMCNeutrino");
    fSaveMCParticle = p.get<bool>("saveMCParticle");
    fSaveMCTrack = p.get<bool>("saveMCTrack");
    fSaveMCShower = p.get<bool>("saveMCShower");
    fSaveSimChannel = p.get<bool>("saveSimChannel");
   
    fOutFileName = p.get<std::string>("outFile");
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::initOutput()
  {
    TDirectory* tmpDir = gDirectory;
    fOutFile = new TFile(fOutFileName.c_str(), "recreate");
    TDirectory* subDir = fOutFile->mkdir("Event");
    subDir->cd();
    fEventTree = new TTree("Sim","Event Tree from Simulation");
    fEventTree->Branch("runNo", &fRun);
    fEventTree->Branch("subRunNo", &fSubRun);
    fEventTree->Branch("eventNo", &fEvent);    
    fEventTree->Branch("eventTime", &fEventTime);

    // MC TRUTH
    fTrueTree = new TTree("MCTruth","MC Truth Tree");
    fTrueTree->Branch("Nmctruth",                  &Nmctruth);

    // MC Neutrino
    fNeutrinoTree = new TTree("MCNeutrino","MC Neutrino Tree from Output of LArG4");
    fNeutrinoTree->Branch("Nmcneutrino",               &Nmcneutrino);
    fNeutrinoTree->Branch("mcneutrino_ccnc",           &mcneutrino_ccnc);
    fNeutrinoTree->Branch("mcneutrino_mode",           &mcneutrino_mode);
    fNeutrinoTree->Branch("mcneutrino_interactionType",&mcneutrino_interactionType);
    fNeutrinoTree->Branch("mcneutrino_target",         &mcneutrino_target);
    fNeutrinoTree->Branch("mcneutrino_hitNuc",         &mcneutrino_hitNuc);
    fNeutrinoTree->Branch("mcneutrino_hitQuark",       &mcneutrino_hitQuark);
    fNeutrinoTree->Branch("mcneutrino_w",              &mcneutrino_w);
    fNeutrinoTree->Branch("mcneutrino_x",              &mcneutrino_x);
    fNeutrinoTree->Branch("mcneutrino_y",              &mcneutrino_y);
    fNeutrinoTree->Branch("mcneutrino_qSqr",           &mcneutrino_qSqr);
    fNeutrinoTree->Branch("mcneutrino_pt",             &mcneutrino_pt);
    fNeutrinoTree->Branch("mcneutrino_theta",          &mcneutrino_theta);

    // MC PARTICLE
    fParticleTree = new TTree("MCParticle","MC Particle Tree from Output of LArG4");
    fParticleTree->Branch("Nmcparticle",               &Nmcparticle);
    fParticleTree->Branch("mcparticle_id",             &mcparticle_id);
    fParticleTree->Branch("mcparticle_statusCode",     &mcparticle_statusCode);
    fParticleTree->Branch("mcparticle_pdg",            &mcparticle_pdg);
    fParticleTree->Branch("mcparticle_mother",         &mcparticle_mother);
    fmcparticle_polarization = new TObjArray();
    fmcparticle_polarization->SetOwner(kTRUE);
    fParticleTree->Branch("mcparticle_polarization",   &fmcparticle_polarization); 
    fParticleTree->Branch("mcparitcle_process",        &mcparticle_process);
    fParticleTree->Branch("mcparticle_endProcess",     &mcparticle_endProcess);           
    fParticleTree->Branch("mcparticle_ndaughters",     &mcparticle_ndaughters);
    fParticleTree->Branch("mcparticle_daughterId",     &mcparticle_daughterId);
    fParticleTree->Branch("mcparticle_ntrajpts",       &mcparticle_ntrajpts);
    fParticleTree->Branch("mcparticle_startXYZT",      &mcparticle_startXYZT, "mcparticle_startXYZT[Nmcparticle][4]");
    fParticleTree->Branch("mcparticle_endXYZT",        &mcparticle_endXYZT, "mcparticle_endXYZT[Nmcparticle][4]");
    fParticleTree->Branch("mcparticle_startMomentum",  &mcparticle_startMomentum, "mcparticle_startMomentum[Nmcparticle][4]");
    fParticleTree->Branch("mcparticle_endMomentum",    &mcparticle_endMomentum, "mcparticle_endMomentum[Nmcparticle][4]");
    fmcparticle_position = new TObjArray();
    fmcparticle_position->SetOwner(kTRUE);
    fParticleTree->Branch("mcparticle_position",       &fmcparticle_position);
    fParticleTree->Branch("mcparticle_vx",             &mcparticle_vx);
    fParticleTree->Branch("mcparticle_vy",             &mcparticle_vy);
    fParticleTree->Branch("mcparticle_vz",             &mcparticle_vz);
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
    fParticleTree->Branch("mcparitlce_endPX",          &mcparticle_endPX);
    fParticleTree->Branch("mcparticle_endPY",          &mcparticle_endPY);
    fParticleTree->Branch("mcparticle_endPZ",          &mcparticle_endPZ);
    fParticleTree->Branch("mcparticle_endE",           &mcparticle_endE);
    fmcparticle_Gvtx = new TObjArray();
    fmcparticle_Gvtx->SetOwner(kTRUE);
    fParticleTree->Branch("mcparticle_Gvtx",           &fmcparticle_Gvtx);
    fParticleTree->Branch("mcparticle_Gvx",            &mcparticle_Gvx);
    fParticleTree->Branch("mcparticle_Gvy",            &mcparticle_Gvy);
    fParticleTree->Branch("mcparticle_Gvz",            &mcparticle_Gvz);
    fParticleTree->Branch("mcparticle_Gvt",            &mcparticle_Gvt);
    fParticleTree->Branch("mcparticle_firsttDaughter", &mcparticle_firstDaughter);
    fParticleTree->Branch("mcparticle_lastDaughter",   &mcparticle_lastDaughter);
    fParticleTree->Branch("mcparticle_rescatter",      &mcparticle_rescatter);
    // trajectory
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
    fmctrack_startPosition = new TObjArray();
    fmctrack_startPosition->SetOwner(kTRUE);
    fTrackTree->Branch("mctrack_startPosition",     &fmctrack_startPosition);
    fmctrack_motherStartPosition = new TObjArray();
    fmctrack_motherStartPosition->SetOwner(kTRUE);
    fTrackTree->Branch("mctrack_motherStartPosition", &fmctrack_motherStartPosition);
    fmctrack_ancestorStartPosition = new TObjArray();
    fmctrack_ancestorStartPosition->SetOwner(kTRUE);
    fTrackTree->Branch("mctrack_ancestorStartPosition", &fmctrack_ancestorStartPosition);
    fmctrack_endPosition = new TObjArray();
    fmctrack_endPosition->SetOwner(kTRUE);
    fTrackTree->Branch("mctrack_endPosition",       &fmctrack_endPosition);
    fmctrack_motherEndPosition = new TObjArray();
    fmctrack_motherEndPosition->SetOwner(kTRUE);
    fTrackTree->Branch("mctrack_motherEndPosition", &fmctrack_motherEndPosition);
    fmctrack_ancestorEndPosition = new TObjArray();
    fmctrack_ancestorEndPosition->SetOwner(kTRUE);
    fTrackTree->Branch("mctrack_ancestorEndPosition", &fmctrack_ancestorEndPosition);
    fmctrack_startMomentum = new TObjArray();
    fmctrack_startMomentum->SetOwner(kTRUE);
    fTrackTree->Branch("mctrack_startMomentum",     &fmctrack_startMomentum);
    fmctrack_motherStartMomentum = new TObjArray();
    fmctrack_motherStartMomentum->SetOwner(kTRUE);
    fTrackTree->Branch("mctrack_motherStartMomentum", &fmctrack_motherStartMomentum);
    fmctrack_ancestorStartMomentum = new TObjArray();
    fmctrack_ancestorStartMomentum->SetOwner(kTRUE);
    fTrackTree->Branch("mctrack_ancestorStartMomentum", &fmctrack_ancestorStartMomentum);
    fmctrack_endMomentum = new TObjArray();
    fmctrack_endMomentum->SetOwner(kTRUE);  
    fTrackTree->Branch("mctrack_endMomentum",       &fmctrack_endMomentum);
    fmctrack_motherEndMomentum = new TObjArray();
    fmctrack_motherEndMomentum->SetOwner(kTRUE);
    fTrackTree->Branch("mctrack_motherEndMomentum", &fmctrack_motherEndMomentum);
    fmctrack_ancestorEndMomentum = new TObjArray();
    fmctrack_ancestorEndMomentum->SetOwner(kTRUE);
    fTrackTree->Branch("mctrack_ancestorEndMomentum", &fmctrack_ancestorEndMomentum);
    fTrackTree->Branch("mctrack_dQdx",             &mctrack_dQdx);
    fTrackTree->Branch("mctrack_dEdx",             &mctrack_dEdx);

    // MC SHOWER
    fShowerTree = new TTree("MCShower","MC Shower Tree from Output of LArG4");
    fShowerTree->Branch("Nmcshower",               &Nmcshower);
    fShowerTree->Branch("mcshower_pdg",            &mcshower_pdg);
    fShowerTree->Branch("mcshower_motherPdg",      &mcshower_motherPdg);
    fShowerTree->Branch("mcshower_ancestorPdg",    &mcshower_ancestorPdg);
    fShowerTree->Branch("mcshower_id",             &mcshower_id);
    fShowerTree->Branch("mcshower_motherId",       &mcshower_motherId);
    fShowerTree->Branch("mcshower_ancestorId",     &mcshower_ancestorId);
    fShowerTree->Branch("mcshower_process",        &mcshower_process);
    fShowerTree->Branch("mcshower_motherProcess",  &mcshower_motherProcess);
    fShowerTree->Branch("mcshower_ancestorProcess", &mcshower_ancestorProcess);
    fmcshower_start = new TObjArray();
    fmcshower_start->SetOwner(kTRUE);
    fShowerTree->Branch("mcshower_start",          &fmcshower_start);
    fmcshower_motherStart = new TObjArray();
    fmcshower_motherStart->SetOwner(kTRUE);
    fShowerTree->Branch("mcshower_motherStart",    &fmcshower_motherStart);
    fmcshower_ancestorStart = new TObjArray();
    fmcshower_ancestorStart->SetOwner(kTRUE);
    fShowerTree->Branch("mcshower_ancestorStart",  &fmcshower_ancestorStart);
    fmcshower_end = new TObjArray();
    fmcshower_end->SetOwner(kTRUE);
    fShowerTree->Branch("mcshower_end",            &fmcshower_end);
    fmcshower_motherEnd = new TObjArray();
    fmcshower_motherEnd->SetOwner(kTRUE);
    fShowerTree->Branch("mcshower_motherEnd",      &fmcshower_motherEnd);
    fmcshower_ancestorEnd = new TObjArray();
    fmcshower_ancestorEnd->SetOwner(kTRUE);
    fShowerTree->Branch("mcshower_ancestorEnd",    &fmcshower_ancestorEnd);
    fmcshower_detprofilePos = new TObjArray();
    fmcshower_detprofilePos->SetOwner(kTRUE);
    fShowerTree->Branch("mcshower_detprofilePos",  &fmcshower_detprofilePos);
    fmcshower_detprofileMom = new TObjArray();
    fmcshower_detprofileMom->SetOwner(kTRUE);
    fShowerTree->Branch("mcshower_detprofileMom",  &fmcshower_detprofileMom);
    fShowerTree->Branch("mcshower_daughterTrackID", &mcshower_daughterTrackID);
    fShowerTree->Branch("mcshower_charge",         &mcshower_charge);
    fShowerTree->Branch("mcshower_dQdx",           &mcshower_dQdx);
    fShowerTree->Branch("mcshower_dEdx",           &mcshower_dEdx);
    fmcshower_startDir = new TObjArray();
    fmcshower_startDir->SetOwner(kTRUE);
    fShowerTree->Branch("mcshower_startDir",       &fmcshower_startDir);

    // SIM CHANNEL
    fSimChTree = new TTree("SimChannel","SimChannel Tree from Output of LArG4");
    fSimChTree->Branch("Nsimchannel",               &Nsimchannel);
    fSimChTree->Branch("simchannel_channelID",      &simchannel_channelID);
    fSimChTree->Branch("simchannel_tdc",            &simchannel_tdc);
    fSimChTree->Branch("simchannel_nEnergyDeposits", &simchannel_nEnergyDeposits);
    fSimChTree->Branch("simchannel_id",             &simchannel_id);
    fSimChTree->Branch("simchannel_nElectrons",     &simchannel_nElectrons);
    fSimChTree->Branch("simchannel_energy",         &simchannel_energy);
    fSimChTree->Branch("simchannel_x",              &simchannel_x);
    fSimChTree->Branch("simchannel_y",              &simchannel_y);
    fSimChTree->Branch("simchannel_z",              &simchannel_z);

    gDirectory = tmpDir;
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::beginJob()
  {
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::endJob()
  {
    TDirectory* tmpDir = gDirectory;
    fOutFile->cd("/Event");
    fEventTree->Write();
    fTrueTree->Write();
    fNeutrinoTree->Write();
    fParticleTree->Write();
    fTrackTree->Write();
    fShowerTree->Write(); cout << "Wrote fShowerTree to file inside endJob() " << endl;
    fSimChTree->Write();
    gDirectory = tmpDir;
    fOutFile->Close();
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::beginRun(const art::Run&)
  {
    mf::LogInfo("CellTreeTruth") << "begin run";
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::analyze(const art::Event& event)
  {
    reset();

    fRun = event.run();
    fSubRun = event.subRun();
    fEvent = event.id().event();
    art::Timestamp ts = event.time();
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    fEventTime = tts.AsDouble();

    if(fSaveMCTruth == true) { processMCTruth(event); fTrueTree->Fill(); }
    if(fSaveMCNeutrino == true) { processMCNeutrino(event); fNeutrinoTree->Fill(); }
    if(fSaveMCParticle == true) { processMCParticle(event); fParticleTree->Fill(); }
    if(fSaveMCTrack == true) { processMCTrack(event); fTrackTree->Fill(); }
    if(fSaveMCShower == true) { processMCShower(event); fShowerTree->Fill(); }
    if(fSaveSimChannel == true) { processSimChannel(event); fSimChTree->Fill(); }
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::reset()
  {
    if(fSaveMCNeutrino == true){
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
    mcparticle_id.clear();
    mcparticle_statusCode.clear();
    mcparticle_pdg.clear();
    mcparticle_mother.clear();
    fmcparticle_polarization->Clear();
    mcparticle_process.clear();
    mcparticle_endProcess.clear();
    mcparticle_ndaughters.clear();
    mcparticle_daughterId.clear();
    mcparticle_ntrajpts.clear();
    Nmcparticle = 0;
    for(int i=0; i<MAX_TRACKS; i++){
      for(int j=0; j<4; j++){
	mcparticle_startXYZT[i][j] = 0;
	mcparticle_endXYZT[i][j] = 0;
	mcparticle_startMomentum[i][j] = 0;
	mcparticle_endMomentum[i][j] = 0;
      }
    }
    fmcparticle_position->Clear();
    mcparticle_vx.clear();
    mcparticle_vy.clear();
    mcparticle_vz.clear();
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
    fmcparticle_Gvtx->Clear();
    mcparticle_Gvx.clear();
    mcparticle_Gvy.clear();
    mcparticle_Gvz.clear();
    mcparticle_Gvt.clear();
    mcparticle_firstDaughter.clear();
    mcparticle_lastDaughter.clear();
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
    fmctrack_startPosition->Clear();
    fmctrack_motherStartPosition->Clear();
    fmctrack_ancestorStartPosition->Clear();
    fmctrack_endPosition->Clear();
    fmctrack_motherEndPosition->Clear();
    fmctrack_ancestorEndPosition->Clear();
    fmctrack_startMomentum->Clear();
    fmctrack_motherStartMomentum->Clear();
    fmctrack_ancestorStartMomentum->Clear();
    fmctrack_endMomentum->Clear();
    fmctrack_motherEndMomentum->Clear();
    fmctrack_ancestorEndMomentum->Clear();
    mctrack_dQdx.clear();
    mctrack_dEdx.clear();
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
    fmcshower_start->Clear();
    fmcshower_motherStart->Clear();
    fmcshower_ancestorStart->Clear();
    fmcshower_end->Clear();
    fmcshower_motherEnd->Clear();
    fmcshower_ancestorEnd->Clear();
    fmcshower_detprofilePos->Clear();
    fmcshower_detprofileMom->Clear();
    mcshower_daughterTrackID.clear();
    mcshower_charge.clear();
    mcshower_dQdx.clear();
    mcshower_dEdx.clear();
    fmcshower_startDir->Clear();
    }

    if(fSaveSimChannel == true){
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
    }
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::processMCTruth(const art::Event& event)
  {
    art::Handle<std::vector<simb::MCTruth> > truthHandle;
    if(! event.getByLabel(fMCTruthLabel, truthHandle)){
      cout << "WARNING: no label " << fMCTruthLabel << endl;
    }
    std::vector<art::Ptr<simb::MCTruth> > truth;
    art::fill_ptr_vector(truth, truthHandle); 

    Nmctruth = (int)truth.size();
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::processMCNeutrino(const art::Event& event)
  {
    art::Handle<std::vector<simb::MCNeutrino> > neutrinoHandle;
    if(! event.getByLabel(fMCNeutrinoLabel, neutrinoHandle)){
      cout << "WARNING: no label " << fMCNeutrinoLabel << endl;
    }
    std::vector<art::Ptr<simb::MCNeutrino> > neutrinos;
    art::fill_ptr_vector(neutrinos, neutrinoHandle);

    Nmcneutrino = (int)neutrinos.size();
    for(auto const& neutrino : neutrinos){
      mcneutrino_ccnc.push_back(neutrino->CCNC());
      mcneutrino_mode.push_back(neutrino->Mode());
      mcneutrino_interactionType.push_back(neutrino->InteractionType());
      mcneutrino_target.push_back(neutrino->Target());
      mcneutrino_hitNuc.push_back(neutrino->HitNuc());
      mcneutrino_hitQuark.push_back(neutrino->HitQuark());
      mcneutrino_w.push_back(neutrino->W());
      mcneutrino_x.push_back(neutrino->X());
      mcneutrino_y.push_back(neutrino->Y());
      mcneutrino_qSqr.push_back(neutrino->QSqr());
      mcneutrino_pt.push_back(neutrino->Pt());
      mcneutrino_theta.push_back(neutrino->Theta());
    }
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::processMCParticle(const art::Event& event)
  {
   
    art::Handle<std::vector<simb::MCParticle> > particleHandle;
    if(! event.getByLabel(fMCParticleLabel, particleHandle)){
      cout << "WARNING: no label " << fMCParticleLabel << endl;
      return;
    }
   
    std::vector<art::Ptr<simb::MCParticle> > particles;
    art::fill_ptr_vector(particles, particleHandle);
   
    //    Nmcparticle = (int)particles.size();
    int i=0; // track index in saved MCParticles
    for(auto const& particle : particles){
      mcparticle_id.push_back(particle->TrackId());
      mcparticle_statusCode.push_back(particle->StatusCode());
      mcparticle_pdg.push_back(particle->PdgCode());
      mcparticle_mother.push_back(particle->Mother());
      TClonesArray *Lpolarization = new TClonesArray("TVector3",0);
      new ((*Lpolarization)[0]) TVector3(particle->Polarization());
      fmcparticle_polarization->Add(Lpolarization);
      mcparticle_process.push_back(particle->Process());
      mcparticle_endProcess.push_back(particle->EndProcess());
      mcparticle_ndaughters.push_back(particle->NumberDaughters());
      for(int a=0; a<particle->NumberDaughters(); a++){
	mcparticle_daughterId.push_back(particle->Daughter(a));
      }
      mcparticle_ntrajpts.push_back(particle->NumberTrajectoryPoints());
      for(int b=0; b<(int)particle->NumberTrajectoryPoints(); b++){
	mcparticle_vx.push_back(particle->Vx(b));
	mcparticle_vy.push_back(particle->Vy(b));
	mcparticle_vz.push_back(particle->Vz(b));
	mcparticle_t.push_back(particle->T(b));
	mcparticle_px.push_back(particle->Px(b));
	mcparticle_py.push_back(particle->Py(b));
	mcparticle_pz.push_back(particle->Pz(b));
	mcparticle_e.push_back(particle->E(b));
	mcparticle_p.push_back(particle->P(b));
	mcparticle_pt.push_back(particle->Pt(b));
      }

      size_t numberTrajectoryPoints = particle->NumberTrajectoryPoints();
      int last = numberTrajectoryPoints - 1;
      const TLorentzVector& positionStart = particle->Position(0);
      const TLorentzVector& positionEnd = particle->Position(last);
      const TLorentzVector& momentumStart = particle->Momentum(0);
      const TLorentzVector& momentumEnd = particle->Momentum(last);
      positionStart.GetXYZT(mcparticle_startXYZT[i]);
      positionEnd.GetXYZT(mcparticle_endXYZT[i]);
      momentumStart.GetXYZT(mcparticle_startMomentum[i]);
      momentumEnd.GetXYZT(mcparticle_endMomentum[i]);

      TClonesArray *Lposition = new TClonesArray("TLorentzVector",numberTrajectoryPoints);
      for(unsigned int j=0; j<numberTrajectoryPoints; j++){
	new ((*Lposition)[j]) TLorentzVector(particle->Position(j));
      }
      fmcparticle_position->Add(Lposition);
      i++;
      if(i == MAX_TRACKS){
	cout << "WARNING: # tracks exceeds MAX_TRACKS " << MAX_TRACKS << endl;
	break;
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

      TClonesArray *LGvtx = new TClonesArray("TLorentzVector",1);
      new ((*LGvtx)[0]) TLorentzVector(particle->GetGvtx());
      fmcparticle_Gvtx->Add(LGvtx);
      mcparticle_Gvx.push_back(particle->Gvx());
      mcparticle_Gvy.push_back(particle->Gvy());
      mcparticle_Gvz.push_back(particle->Gvz());
      mcparticle_Gvt.push_back(particle->Gvt());
      //mcparticle_firstDaughter.push_back(particle->FirstDaughter());
      //mcparticle_lastDaughter.push_back(particle->LastDaughter()); // -->(first/last daughter) failure
      mcparticle_rescatter.push_back(particle->Rescatter());
      // trajectory
      mcparticle_weight.push_back(particle->Weight());
    } 
    Nmcparticle = i;
  }

  //-------------------------------------------------------------------
  void CellTreeTruth::processMCTrack(const art::Event& event)
  {
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
      TClonesArray *LstartP = new TClonesArray("TLorentzVector",1);
      new ((*LstartP)[0]) TLorentzVector(track->Start().Position());
      fmctrack_startPosition->Add(LstartP);
      TClonesArray *LendP = new TClonesArray("TLorentzVector",1);
      new ((*LendP)[0]) TLorentzVector(track->End().Position());
      fmctrack_endPosition->Add(LendP);
      TClonesArray *LstartM = new TClonesArray("TLorentzVector",1);
      new ((*LstartM)[0]) TLorentzVector(track->Start().Momentum());
      fmctrack_startMomentum->Add(LstartM);
      TClonesArray *LendM = new TClonesArray("TLorentzVector",1);
      new ((*LendM)[0]) TLorentzVector(track->End().Momentum());
      fmctrack_endMomentum->Add(LendM);
      mctrack_dQdx.push_back(track->dQdx());
      mctrack_dEdx.push_back(track->dEdx());
      mctrack_motherPdg.push_back(track->MotherPdgCode());
      mctrack_motherId.push_back(track->MotherTrackID());
      mctrack_motherProcess.push_back(track->MotherProcess());

      TClonesArray *LmomstartP = new TClonesArray("TLorentzVector",1);
      new ((*LmomstartP)[0]) TLorentzVector(track->MotherStart().Position());
      fmctrack_motherStartPosition->Add(LmomstartP);
      TClonesArray *LmomendP = new TClonesArray("TLorentzVector",1);
      new ((*LmomendP)[0]) TLorentzVector(track->MotherEnd().Position());
      fmctrack_motherEndPosition->Add(LmomendP);
      TClonesArray *LmomstartM = new TClonesArray("TLorentzVector",1);
      new ((*LmomstartM)[0]) TLorentzVector(track->MotherStart().Momentum());
      fmctrack_motherStartMomentum->Add(LmomstartM);
      TClonesArray *LmomendM = new TClonesArray("TLorentzVector",1);
      new ((*LmomendM)[0]) TLorentzVector(track->MotherEnd().Momentum());
      fmctrack_motherEndMomentum->Add(LmomendM);

      mctrack_ancestorPdg.push_back(track->AncestorPdgCode());
      mctrack_ancestorId.push_back(track->AncestorTrackID());
      mctrack_ancestorProcess.push_back(track->AncestorProcess());
      TClonesArray *LancstartP = new TClonesArray("TLorentzVector",1);
      new ((*LancstartP)[0]) TLorentzVector(track->AncestorStart().Position());
      fmctrack_ancestorStartPosition->Add(LancstartP);
      TClonesArray *LancendP = new TClonesArray("TLorentzVector",1);
      new ((*LancendP)[0]) TLorentzVector(track->AncestorEnd().Position());
      fmctrack_ancestorEndPosition->Add(LancendP);
      TClonesArray *LancstartM = new TClonesArray("TLorentzVector",1);
      new ((*LancstartM)[0]) TLorentzVector(track->AncestorStart().Momentum());
      fmctrack_ancestorStartMomentum->Add(LancstartM);
      TClonesArray *LancendM = new TClonesArray("TLorentzVector",1);
      new ((*LancendM)[0]) TLorentzVector(track->AncestorEnd().Momentum());
      fmctrack_ancestorEndMomentum->Add(LancendM);         
      
    }
  }
  
  //-------------------------------------------------------------------
  void CellTreeTruth::processMCShower(const art::Event& event)
  {
    art::Handle< std::vector<sim::MCShower> > showerHandle;
    if(! event.getByLabel(fMCShowerLabel, showerHandle)){
      cout << "WARNING: no label " << fMCShowerLabel << endl;
      return;
    }
    std::vector<art::Ptr<sim::MCShower> > showers;
    art::fill_ptr_vector(showers, showerHandle);
    
    Nmcshower = (int)showers.size();
    for(auto const& shower : showers){
      // origin
      mcshower_pdg.push_back(shower->PdgCode());
      mcshower_id.push_back(shower->TrackID());
      mcshower_process.push_back(shower->Process());
      TClonesArray *Lstart = new TClonesArray("TLorentzVector",1);
      new ((*Lstart)[0]) TLorentzVector(shower->Start().Position());
      fmcshower_start->Add(Lstart);
      TClonesArray *Lend = new TClonesArray("TLorentzVector",1);
      new ((*Lend)[0]) TLorentzVector(shower->End().Position());
      fmcshower_end->Add(Lend);
      mcshower_dQdx.push_back(shower->dQdx());
      mcshower_dEdx.push_back(shower->dEdx());
      TClonesArray *Lsd = new TClonesArray("TVector3",1);
      new ((*Lsd)[0]) TVector3(shower->StartDir());
      fmcshower_startDir->Add(Lsd);
      TClonesArray *LdpPos = new TClonesArray("TLorentzVector", 1);
      new ((*LdpPos)[0]) TLorentzVector(shower->DetProfile().Position());
      fmcshower_detprofilePos->Add(LdpPos);
      TClonesArray *LdpMom = new TClonesArray("TLorentzVector", 1);
      new ((*LdpMom)[0]) TLorentzVector(shower->DetProfile().Momentum());
      fmcshower_detprofileMom->Add(LdpMom);
      mcshower_motherPdg.push_back(shower->MotherPdgCode());
      mcshower_motherId.push_back(shower->MotherTrackID());
      mcshower_motherProcess.push_back(shower->MotherProcess());
      TClonesArray *Lmomstart = new TClonesArray("TLorentzVector",1);
      new ((*Lmomstart)[0]) TLorentzVector(shower->MotherStart().Position());
      fmcshower_motherStart->Add(Lmomstart);
      TClonesArray *Lmomend = new TClonesArray("TLorentzVector",1);
      new ((*Lmomend)[0]) TLorentzVector(shower->MotherEnd().Position());
      fmcshower_motherEnd->Add(Lmomend);
      mcshower_ancestorPdg.push_back(shower->AncestorPdgCode());
      mcshower_ancestorId.push_back(shower->AncestorTrackID());
      mcshower_ancestorProcess.push_back(shower->AncestorProcess());
      TClonesArray *Lancstart = new TClonesArray("TLorentzVector",1);
      new ((*Lancstart)[0]) TLorentzVector(shower->AncestorStart().Position());
      fmcshower_ancestorStart->Add(Lancstart);
      TClonesArray *Lancend = new TClonesArray("TLorentzVector",1);
      new ((*Lancend)[0]) TLorentzVector(shower->AncestorEnd().Position());
      fmcshower_ancestorEnd->Add(Lancend);
    }
  }  

  //---------------------------------------------------------------------
  void CellTreeTruth::processSimChannel(const art::Event& event)
  {
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

  DEFINE_ART_MODULE(CellTreeTruth)
} // namespace wc
#endif
