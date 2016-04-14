#ifndef NDecayTrackingEff_Module
#define NDecayTrackingEff_Module

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/RecoBase/Track.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"
#include "larcore/SimpleTypesAndConstants/geo_types.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/RawData/ExternalTrigger.h"
// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

#define MAX_TRACKS 1000
using namespace std;

//========================================================================

namespace DUNE{

class NDecayTrackingEff : public art::EDAnalyzer {
public:

    explicit NDecayTrackingEff(fhicl::ParameterSet const& pset);
    virtual ~NDecayTrackingEff();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);

    void processEff(const art::Event& evt, bool &isFiducial);
    void truthMatcher( std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac);
    double truthLength( const simb::MCParticle *MCparticle );
    bool insideFV(double vertex[4]);
    void doEfficiencies();
    void reset();

private:

    // the parameters we'll read from the .fcl
    std::string fMCTruthModuleLabel;
    std::string fTrackModuleLabel;
    int		fLeptonPDGcode;
    bool	fSaveMCTree; 

    TTree *fEventTree;
    TTree *fHitsTree;

    TH1D *h_Pmu_den;
    TH1D *h_Pmu_num;
    TH1D *h_Pkaon_den;
    TH1D *h_Pkaon_num;

    TH1D *h_Efrac_lepton;     
    TH1D *h_Efrac_kaon;     
    TH1D *h_trackRes_lepton;
    TH1D *h_trackRes_kaon;

    TEfficiency* h_Eff_Pmu = 0;
    TEfficiency* h_Eff_Pkaon = 0;

    // Event 
    int Event;
    int Run;
    int SubRun;

    //MC truth
    int    MC_incoming_PDG;
    int    MC_lepton_PDG;
    int    MC_isCC;
    int    MC_channel;
    int    MC_target;
    double MC_Q2;
    double MC_W;
    double MC_vertex[4];
    double MC_incoming_P[4];
    double MC_lepton_startMomentum[4]; 
    double MC_lepton_endMomentum[4]; 
    double MC_lepton_startXYZT[4];
    double MC_lepton_endXYZT[4];
    double MC_leptonP; 
    int    MC_leptonID;
    int    MC_kaonID;
    double MC_kaonP;
    double MC_kaon_startMomentum[4]; 
    double MC_kaon_endMomentum[4]; 
    double MC_kaon_startXYZT[4];
    double MC_kaon_endXYZT[4];
 
    int    MC_LeptonTrack;
    int    MC_kaonTrack;
 
    int    MC_Ntrack;  
    int    MC_id[MAX_TRACKS];  
    int    MC_pdg[MAX_TRACKS]; 
    int    MC_mother[MAX_TRACKS];  
    double MC_startXYZT[MAX_TRACKS][4]; 
    double MC_endXYZT[MAX_TRACKS][4];  
    double MC_startMomentum[MAX_TRACKS][4]; 
    double MC_endMomentum[MAX_TRACKS][4];  

    double Reco_EfracLepton;
    double Reco_LengthRes; 
    double Reco_LengthReskaon;

    int    n_recoTrack;

    float fFidVolCutX;
    float fFidVolCutY;
    float fFidVolCutZ;

    float fFidVolXmin;
    float fFidVolXmax;
    float fFidVolYmin;
    float fFidVolYmax;
    float fFidVolZmin;
    float fFidVolZmax;

    detinfo::DetectorProperties const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    detinfo::DetectorClocks const *ts = lar::providerFrom<detinfo::DetectorClocksService>();
    double XDriftVelocity = detprop->DriftVelocity()*1e-3; //cm/ns
    double WindowSize     = detprop->NumberTimeSamples() * ts->TPCClock().TickPeriod() * 1e3;

    art::ServiceHandle<geo::Geometry> geom;
 }; // class NDecayTrackingEff


//========================================================================
NDecayTrackingEff::NDecayTrackingEff(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    reconfigure(parameterSet);
}
//========================================================================
NDecayTrackingEff::~NDecayTrackingEff(){
  //destructor
}
//========================================================================
void NDecayTrackingEff::reconfigure(fhicl::ParameterSet const& p){

    fTrackModuleLabel    = p.get<std::string>("TrackModuleLabel");
    fLeptonPDGcode	 = p.get<int>("LeptonPDGcode");
    fSaveMCTree		 = p.get<bool>("SaveMCTree");
    fFidVolCutX          = p.get<float>("FidVolCutX");
    fFidVolCutY          = p.get<float>("FidVolCutY");
    fFidVolCutZ          = p.get<float>("FidVolCutZ");
}
//========================================================================
void NDecayTrackingEff::beginJob(){
  std::cout<<"job begin..."<<std::endl;
  // Get geometry.
  auto const* geo = lar::providerFrom<geo::Geometry>();
  // Define histogram boundaries (cm).
  // For now only draw cryostat=0.
  double minx = 1e9;
  double maxx = -1e9;
  double miny = 1e9;
  double maxy = -1e9;
  double minz = 1e9;
  double maxz = -1e9;
  for (size_t i = 0; i<geo->NTPC(); ++i){
    double local[3] = {0.,0.,0.};
    double world[3] = {0.,0.,0.};
    const geo::TPCGeo &tpc = geo->TPC(i);
    tpc.LocalToWorld(local,world);
    if (minx>world[0]-geo->DetHalfWidth(i))
      minx = world[0]-geo->DetHalfWidth(i);
    if (maxx<world[0]+geo->DetHalfWidth(i))
      maxx = world[0]+geo->DetHalfWidth(i);
    if (miny>world[1]-geo->DetHalfHeight(i))
      miny = world[1]-geo->DetHalfHeight(i);
    if (maxy<world[1]+geo->DetHalfHeight(i))
      maxy = world[1]+geo->DetHalfHeight(i);
    if (minz>world[2]-geo->DetLength(i)/2.)
      minz = world[2]-geo->DetLength(i)/2.;
    if (maxz<world[2]+geo->DetLength(i)/2.)
      maxz = world[2]+geo->DetLength(i)/2.;
  }

  fFidVolXmin = minx + fFidVolCutX;
  fFidVolXmax = maxx - fFidVolCutX;
  fFidVolYmin = miny + fFidVolCutY;
  fFidVolYmax = maxy - fFidVolCutY;
  fFidVolZmin = minz + fFidVolCutZ;
  fFidVolZmax = maxz - fFidVolCutZ;

  std::cout<<"Fiducial volume:"<<"\n"
	   <<fFidVolXmin<<"\t< x <\t"<<fFidVolXmax<<"\n"
	   <<fFidVolYmin<<"\t< y <\t"<<fFidVolYmax<<"\n"
	   <<fFidVolZmin<<"\t< z <\t"<<fFidVolZmax<<"\n";

  art::ServiceHandle<art::TFileService> tfs;

  double Pbins[21] ={0,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0};
  

  h_Pmu_den = tfs->make<TH1D>("h_Pmu_den","Muon Momentum; Muon Momentum (GeV); Tracking Efficiency",20,Pbins);
  h_Pmu_num = tfs->make<TH1D>("h_Pmu_num","Muon Momentum; Muon Momentum (GeV); Tracking Efficiency",20,Pbins);
  h_Pkaon_den = tfs->make<TH1D>("h_Pkaon_den","Kaon; Kaon Momentum (GeV); Tracking Efficiency", 20, Pbins);
  h_Pkaon_num = tfs->make<TH1D>("h_Pkaon_num","Kaon; Kaon Momentum (GeV); Tracking Efficiency", 20, Pbins);
  h_Pmu_den->Sumw2();
  h_Pmu_num->Sumw2();
  h_Pkaon_den->Sumw2();
  h_Pkaon_num->Sumw2();

  h_Efrac_lepton = tfs->make<TH1D>("h_Efrac_lepton","Efrac Lepton; Track Energy fraction;",60,0,1.2);
  h_Efrac_kaon = tfs->make<TH1D>("h_Efrac_kaon","Efrac Kaon; Track Energy fraction;",60,0,1.2);
  h_trackRes_lepton = tfs->make<TH1D>("h_trackRes_lepton", "Muon Residual; Truth length - Reco length (cm);",200,-100,100);
  h_trackRes_kaon = tfs->make<TH1D>("h_trackRes_kaon","Kaonn Residual; Truth length - Reco length (cm);",200,-100,100);
  h_Efrac_lepton->Sumw2();
  h_Efrac_kaon->Sumw2();
  h_trackRes_lepton->Sumw2();
  h_trackRes_kaon->Sumw2();

  if( fSaveMCTree ){
    fEventTree = tfs->make<TTree>("Event", "Event Tree from Sim & Reco");
    fEventTree->Branch("eventNo", &Event);
    fEventTree->Branch("runNo", &Run);
    fEventTree->Branch("subRunNo", &SubRun);
    fEventTree->Branch("mc_incoming_PDG", &MC_incoming_PDG);
    fEventTree->Branch("mc_lepton_PDG", &MC_lepton_PDG);
    fEventTree->Branch("mc_isCC", &MC_isCC);
    fEventTree->Branch("mc_target", &MC_target);
    fEventTree->Branch("mc_channel", &MC_channel);
    fEventTree->Branch("mc_Q2", &MC_Q2);
    fEventTree->Branch("mc_W", &MC_W);
    fEventTree->Branch("mc_vertex", &MC_vertex, "mc_vertex[4]/D");
    fEventTree->Branch("mc_incoming_P", &MC_incoming_P, "mc_incoming_P[4]/D");
    fEventTree->Branch("mc_lepton_startMomentum", &MC_lepton_startMomentum, "mc_lepton_startMomentum[4]/D");
    fEventTree->Branch("mc_lepton_endMomentum", &MC_lepton_endMomentum, "mc_lepton_endMomentum[4]/D");
    fEventTree->Branch("mc_lepton_startXYZT", &MC_lepton_startXYZT, "mc_lepton_startXYZT[4]/D");
    fEventTree->Branch("mc_lepton_endXYZT", &MC_lepton_endXYZT, "mc_lepton_endXYZT[4]/D");
    fEventTree->Branch("mc_kaon_startMomentum", &MC_kaon_startMomentum, "mc_kaon_startMomentum[4]/D");
    fEventTree->Branch("mc_kaon_endMomentum", &MC_kaon_endMomentum, "mc_kaon_endMomentum[4]/D");
    fEventTree->Branch("mc_kaon_startXYZT", &MC_kaon_startXYZT, "mc_kaon_startXYZT[4]/D");
    fEventTree->Branch("mc_kaon_endXYZT", &MC_kaon_endXYZT, "mc_kaon_endXYZT[4]/D");
    fEventTree->Branch("mc_leptonP", &MC_leptonP, "mc_leptonP/D");
    fEventTree->Branch("mc_leptonID", &MC_leptonID);
    fEventTree->Branch("mc_kaonID", &MC_kaonID);
    fEventTree->Branch("mc_leptonTrack", &MC_LeptonTrack);
    fEventTree->Branch("mc_kaonTrack", &MC_kaonTrack);
    fEventTree->Branch("mc_Ntrack", &MC_Ntrack);  // number of particles 
    fEventTree->Branch("mc_id", &MC_id, "mc_id[mc_Ntrack]/I");  
    fEventTree->Branch("mc_pdg", &MC_pdg, "mc_pdg[mc_Ntrack]/I"); 
    fEventTree->Branch("mc_mother", &MC_mother, "mc_mother[mc_Ntrack]/I"); 
    fEventTree->Branch("mc_startXYZT", &MC_startXYZT, "mc_startXYZT[mc_Ntrack][4]/D");  
    fEventTree->Branch("mc_endXYZT", &MC_endXYZT, "mc_endXYZT[mc_Ntrack][4]/D"); 
    fEventTree->Branch("mc_startMomentum", &MC_startMomentum, "mc_startMomentum[mc_Ntrack][4]/D");  
    fEventTree->Branch("mc_endMomentum", &MC_endMomentum, "mc_endMomentum[mc_Ntrack][4]/D"); 
    
    fEventTree->Branch("reco_efracLepton", &Reco_EfracLepton);
    fEventTree->Branch("reco_lengthReskaon", &Reco_LengthReskaon);
    fEventTree->Branch("reco_MC_LeptonTrack", &MC_LeptonTrack);
    fEventTree->Branch("reco_MC_kaonTrack", &MC_kaonTrack);
    
  }


}
//========================================================================
void NDecayTrackingEff::endJob(){
     
  doEfficiencies();

}
//========================================================================
void NDecayTrackingEff::beginRun(const art::Run& /*run*/){
  mf::LogInfo("NDecayTrackingEff")<<"begin run..."<<std::endl;
}
//========================================================================
void NDecayTrackingEff::analyze( const art::Event& event ){
    if (event.isRealData()) return;
    reset();

    Event  = event.id().event(); 
    Run    = event.run();
    SubRun = event.subRun();
    bool isFiducial = false;
    processEff(event, isFiducial);
    if( fSaveMCTree ){
      if(isFiducial && MC_isCC == 1) fEventTree->Fill();
    }
}
//========================================================================
void NDecayTrackingEff::processEff( const art::Event& event, bool &isFiducial){

    //!save FS particles
    simb::MCParticle *MClepton = NULL; 
    double leptonLength =0.0;
    simb::MCParticle *MCkaon = NULL;
    double kaonLength =0.0;  
    art::ServiceHandle<cheat::BackTracker> bt;
    const sim::ParticleList& plist = bt->ParticleList();
    simb::MCParticle *particle=0;
    int i=0; // particle index
    MC_Ntrack = plist.size();

    for( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
       particle = ipar->second;
       MC_id[i] = particle->TrackId();
       MC_pdg[i] = particle->PdgCode();
       MC_mother[i] = particle->Mother();
       const TLorentzVector& positionStart = particle->Position(0);
       const TLorentzVector& positionEnd   = particle->EndPosition();
       const TLorentzVector& momentumStart = particle->Momentum(0);
       const TLorentzVector& momentumEnd   = particle->EndMomentum();
       positionStart.GetXYZT(MC_startXYZT[i]);
       positionEnd.GetXYZT(MC_endXYZT[i]);
       momentumStart.GetXYZT(MC_startMomentum[i]);
       momentumEnd.GetXYZT(MC_endMomentum[i]);
       if( particle->Mother() == 0 && particle->PdgCode() == 321 ){   //save primary Kaon
         const TLorentzVector& kaon_momentum =particle->Momentum(0); 
         const TLorentzVector& kaon_position =particle->Position(0); 
         const TLorentzVector& kaon_positionEnd   = particle->EndPosition();
         const TLorentzVector& kaon_momentumEnd   = particle->EndMomentum();
         kaon_momentum.GetXYZT(MC_lepton_startMomentum);
         kaon_position.GetXYZT(MC_lepton_startXYZT);
         kaon_positionEnd.GetXYZT(MC_lepton_endXYZT);
         kaon_momentumEnd.GetXYZT(MC_lepton_endMomentum);
         kaon_position.GetXYZT(MC_vertex); //Primary vertex
         MC_kaonID = particle->TrackId();          
         MC_kaonP = sqrt(pow(particle->Momentum().Px(),2)+pow(particle->Momentum().Py(),2)+pow(particle->Momentum().Pz(),2));
         kaonLength = truthLength(particle); 
         MCkaon = particle;
       }
       else if( particle->Mother() ==2 && particle->Process() =="Decay" && particle->PdgCode() == fLeptonPDGcode ){  //primary lepton
         const TLorentzVector& lepton_momentum =particle->Momentum(0); 
         const TLorentzVector& lepton_position =particle->Position(0); 
         const TLorentzVector& lepton_positionEnd   = particle->EndPosition();
         const TLorentzVector& lepton_momentumEnd   = particle->EndMomentum();
         lepton_momentum.GetXYZT(MC_lepton_startMomentum);
         lepton_position.GetXYZT(MC_lepton_startXYZT);
         lepton_positionEnd.GetXYZT(MC_lepton_endXYZT);
         lepton_momentumEnd.GetXYZT(MC_lepton_endMomentum);
         MC_leptonID = particle->TrackId();
         MC_leptonP = sqrt(pow(MC_lepton_startMomentum[0],2)+pow(MC_lepton_startMomentum[1],2)+pow(MC_lepton_startMomentum[2],2));
         leptonLength = truthLength(particle); 
         MClepton = particle;
       }
       i++; //paticle index
    } 
    isFiducial =insideFV( MC_vertex );
    if( !isFiducial ) return;
    //save events within the fiducial volume  
       if( MClepton ){
         h_Pmu_den->Fill(MC_leptonP);
       }
       if( MCkaon ){
         h_Pkaon_den->Fill(MC_kaonP);
       }
   
    //========================================================================
    //========================================================================
    // Reco  stuff
    //========================================================================
    //========================================================================
    art::Handle< std::vector<recob::Track> > trackListHandle;
    if(! event.getByLabel(fTrackModuleLabel, trackListHandle)) return;
    std::vector<art::Ptr<recob::Track> > tracklist;
    art::fill_ptr_vector(tracklist, trackListHandle);
    n_recoTrack = tracklist.size();

    art::FindManyP<recob::Hit> track_hits(trackListHandle, event, fTrackModuleLabel);
    if( n_recoTrack == 0 ){
      LOG_DEBUG("NDecayTrackingEff")<<"There are no reco tracks... bye";
      return; 
    }
    LOG_DEBUG("NDecayTrackingEff")<<"Found this many reco tracks "<<n_recoTrack;

    double Efrac_lepton =0.0;
    double Efrac_kaon =0.0;
    double trackLength_lepton =0.0;
    double trackLength_kaon =0.0;
    const simb::MCParticle *MClepton_reco = NULL; 
    const simb::MCParticle *MCkaon_reco = NULL;
    for(int i=0; i<n_recoTrack; i++) {
       art::Ptr<recob::Track> track = tracklist[i];
       //const TVector3 tmp_track_vtx = track->Vertex();
       //double track_vtx[4] ={tmp_track_vtx[0], tmp_track_vtx[1], tmp_track_vtx[2], -999};
       //bool track_isFiducial = insideFV( track_vtx );
       //if( !track_isFiducial ) continue;
       std::vector<art::Ptr<recob::Hit>> all_trackHits = track_hits.at(i);  
       double tmpEfrac = 0;
       const simb::MCParticle *particle;
       truthMatcher( all_trackHits, particle, tmpEfrac );
       if (!particle) continue;
       if(  (particle->PdgCode() == fLeptonPDGcode) && (particle->TrackId() == MC_leptonID) ){
         //save the best track ... based on Efrac if there is more than one track 
         if( tmpEfrac > Efrac_lepton ){
           Efrac_lepton = tmpEfrac;
           trackLength_lepton = track->Length(); 
           MClepton_reco = particle;
         }
       }
       else if( (particle->PdgCode() == 321) && (particle->TrackId() == MC_kaonID) ){
         //save the best track ... based on Efrac if there is more than one track 
         if( tmpEfrac > Efrac_kaon ){
           Efrac_kaon = tmpEfrac;
           trackLength_kaon = track->Length();
           MCkaon_reco = particle;
         }
       }
    }

    if( MClepton_reco && MClepton  ){
        MC_LeptonTrack = 1;
        h_Pmu_num->Fill(MC_leptonP);
        h_Efrac_lepton->Fill(Efrac_lepton);
        h_trackRes_lepton->Fill(leptonLength-trackLength_lepton);
    }
    if( MCkaon_reco && MCkaon ){
        MC_kaonTrack = 1;
        h_Pkaon_num->Fill(MC_kaonP);     
        h_Efrac_kaon->Fill(Efrac_kaon);
        h_trackRes_kaon->Fill(kaonLength-trackLength_kaon);
    }


}
//========================================================================
void NDecayTrackingEff::truthMatcher( std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac){

    //std::cout<<"truthMatcher..."<<std::endl;
    art::ServiceHandle<cheat::BackTracker> bt;
    std::map<int,double> trkID_E;
    for(size_t j = 0; j < track_hits.size(); ++j){
       art::Ptr<recob::Hit> hit = track_hits[j];
       std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackID(hit);
       for(size_t k = 0; k < TrackIDs.size(); k++){
          trkID_E[TrackIDs[k].trackID] += TrackIDs[k].energy;
       }            
    }
    double max_E = -999.0;
    double total_E = 0.0;
    int TrackID = -999;
    double partial_E =0.0; // amount of energy deposited by the particle that deposited more energy... tomato potato... blabla
    //!if the collection of hits have more than one particle associate save the particle w/ the highest energy deposition 
    //!since we are looking for muons/pions/protons this should be enough 
    if( !trkID_E.size() ) {
      MCparticle = 0;
      return; //Ghost track???
    }
    for(std::map<int,double>::iterator ii = trkID_E.begin(); ii!=trkID_E.end(); ++ii){
       total_E += ii->second;
       if((ii->second)>max_E){
         partial_E = ii->second;
         max_E = ii->second;
         TrackID = ii->first;
       }
    } 

    MCparticle = bt->TrackIDToParticle(TrackID);
    Efrac = partial_E/total_E;
}

//========================================================================
 double NDecayTrackingEff::truthLength( const simb::MCParticle *MCparticle ){
   //calculate the truth length considering only the part that is inside the TPC
   int numberTrajectoryPoints = MCparticle->NumberTrajectoryPoints();
   double TPCLengthHits[numberTrajectoryPoints];
   int FirstHit=0, LastHit=0;
   double TPCLength = 0.0;
   bool BeenInVolume = false;

   for(int MCHit=0; MCHit < numberTrajectoryPoints; ++MCHit) {
      const TLorentzVector& tmpPosition= MCparticle->Position(MCHit);
      double const tmpPosArray[]={tmpPosition[0],tmpPosition[1],tmpPosition[2]};
      if (MCHit!=0) TPCLengthHits[MCHit] = sqrt( pow( (MCparticle->Vx(MCHit-1)-MCparticle->Vx(MCHit)),2)+ pow( (MCparticle->Vy(MCHit-1)-MCparticle->Vy(MCHit)),2)+ pow( (MCparticle->Vz(MCHit-1)-MCparticle->Vz(MCHit)),2));
      geo::TPCID tpcid = geom->FindTPCAtPosition(tmpPosArray);
      if(tpcid.isValid) { 
	// -- Check if hit is within drift window...
	geo::CryostatGeo const& cryo = geom->Cryostat(tpcid.Cryostat);
	geo::TPCGeo      const& tpc  = cryo.TPC(tpcid.TPC);
        double XPlanePosition      = tpc.PlaneLocation(0)[0];
        double DriftTimeCorrection = fabs( tmpPosition[0] - XPlanePosition ) / XDriftVelocity;
 	double TimeAtPlane         = MCparticle->T() + DriftTimeCorrection;
        if( TimeAtPlane < detprop->TriggerOffset() || TimeAtPlane > detprop->TriggerOffset() + WindowSize ) continue;
	LastHit = MCHit;
        if( !BeenInVolume ) {
          BeenInVolume = true;
          FirstHit = MCHit;
        }
      }
    }
    for (int Hit = FirstHit+1; Hit <= LastHit; ++Hit ) TPCLength += TPCLengthHits[Hit];

    return TPCLength;
}
//========================================================================
bool NDecayTrackingEff::insideFV( double vertex[4]){ 

     double x = vertex[0];
     double y = vertex[1];
     double z = vertex[2];

     if (x>fFidVolXmin && x<fFidVolXmax&&
	 y>fFidVolYmin && y<fFidVolYmax&&
	 z>fFidVolZmin && z<fFidVolZmax)
       return true;
     else
       return false;
}
//========================================================================
void NDecayTrackingEff::doEfficiencies(){

   art::ServiceHandle<art::TFileService> tfs;

   if(TEfficiency::CheckConsistency(*h_Pmu_num,*h_Pmu_den)){ 
     h_Eff_Pmu = tfs->make<TEfficiency>(*h_Pmu_num,*h_Pmu_den);
     TGraphAsymmErrors *grEff_Pmu = h_Eff_Pmu->CreateGraph();
     grEff_Pmu->Write("grEff_Pmu");
     h_Eff_Pmu->Write("h_Eff_Pmu");
   }
   if(TEfficiency::CheckConsistency(*h_Pkaon_num,*h_Pkaon_den)){
     h_Eff_Pkaon = tfs->make<TEfficiency>(*h_Pkaon_num,*h_Pkaon_den);
     TGraphAsymmErrors *grEff_Pkaon = h_Eff_Pkaon->CreateGraph();
     grEff_Pkaon->Write("grEff_Pkaon");
     h_Eff_Pkaon->Write("h_Eff_Pkaonn");
   }

}
//========================================================================
void NDecayTrackingEff::reset(){

   MC_incoming_PDG = -999;
   MC_lepton_PDG =-999;
   MC_isCC =-999;
   MC_channel =-999;
   MC_target =-999;
   MC_Q2 =-999.0;
   MC_W =-999.0;
   MC_leptonP = -999.0;
   MC_leptonID = -999;
   MC_kaonID = -999;
   MC_kaonP = -999.0;
   MC_LeptonTrack = -999;
   MC_kaonTrack =-999;
 
   for(int i = 0; i<4; i++){
      MC_vertex[i] = 0;
      MC_incoming_P[i] =  0;
      MC_lepton_startMomentum[i] = 0;
      MC_lepton_endMomentum[i] = 0;
      MC_lepton_startXYZT[i] = 0;
      MC_lepton_endXYZT[i] = 0;
      MC_kaon_startMomentum[i] = 0;
      MC_kaon_endMomentum[i] = 0;
      MC_kaon_startXYZT[i] = 0;
      MC_kaon_endXYZT[i] = 0;
 
   }

   Reco_EfracLepton = -999.0; 
   Reco_LengthRes = -999.0;
   Reco_LengthReskaon = -999.0;
  
   for(int i=0; i<MAX_TRACKS; i++) {
       MC_id[i] = 0;
       MC_pdg[i] = 0;
       MC_mother[i] = 0;
       for(int j=0; j<4; j++) {
          MC_startXYZT[i][j]      = 0;
          MC_endXYZT[i][j]        = 0;
          MC_startMomentum[i][j] = 0;
          MC_endMomentum[i][j]   = 0;
       }
    }

}
//========================================================================
DEFINE_ART_MODULE(NDecayTrackingEff)

} 

#endif // NDecayTrackingEff_Module
