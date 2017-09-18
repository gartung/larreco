/**
 *  @file   TrackRecoAna_module.cc
 *
 *  @brief  Ana analysis module aimed at comparing PFParticles to MC Truth
 *          with a specific emphasis on comparisons to single particles
 *
 *          Steal a lot from Tracy's PFParticleMcAna_module.cc
 */

#ifndef TrackRecoAna_module
#define TrackRecoAna_module

// LArSoft includes
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCHitCollection.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"

#include "larreco/RecoAna/MCParticleToPFParticle.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/cpu_timer.h"
#include "cetlib/exception.h"

// ROOT
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// C++ Includes
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>
#include <memory>

namespace TrackRecoAna
{
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition

class TrackRecoAna : public art::EDAnalyzer
{
public:

    // Standard constructor and destructor for an ART module.
    explicit TrackRecoAna( fhicl::ParameterSet const& pset );
    virtual ~TrackRecoAna();

    // This method is called once, at the start of the job. In this
    // example, it will define the histograms and n-tuples we'll write.
    void beginJob();

    // This method is called once, at the start of each run. It's a
    // good place to read databases or files that may have
    // run-dependent information.
    void beginRun( const art::Run& run );

    // This method reads in any parameters from the .fcl files. This
    // method is called 'reconfigure' because it might be called in the
    // middle of a job; e.g., if the user changes parameter values in an
    // interactive event display.
    void reconfigure( fhicl::ParameterSet const& pset );

    // The analysis routine, called once per event.
    void analyze ( const art::Event& evt );

private:

    // This method is meant to be called at the start of the event
    void PrepareEvent( const art::Event& evt );

    // The parameters we'll read from the .fcl file.
    std::string fSimulationProducerLabel;    //> The name of the producer that tracked simulated particles through the detector
    std::string fMCHitCollectionModuleLabel; //> The name of the producer of MCHits
    std::string fPFParticleProducerLabel;    //> The name of the produder of the PFParticle hierarchy
    std::string fHitProducerLabel;           //> The name of the producer that created hits
    std::string fClusterProducerLabel;       //> The name of the producer that created clusters
    // std::string fTrackProducerLabel;         //> The name of the producer that created the tracks

    // Maximum number of entries per row in our tuple
    int         fMaxEntries;

    // The variables that will go into the n-tuple.
    // Start with basic run information
    Int_t    fEvent;
    Int_t    fRun;
    Int_t    fSubRun;

    // Event-wise information
    Int_t    fNHits;
    Int_t    fNMCHits;
    Int_t    fNUniqueHits;
    // Int_t    fNNoiseHits;
    // Int_t    fNegTrackIds;
    Int_t    fNMCParticles;
    Int_t    fNPFParticles;
    Int_t    fNMatchedPFParticles;
    Double_t fRecoPVx;
    Double_t fRecoPVy;
    Double_t fRecoPVz;

    // MCParticle-wise information
    std::vector< Int_t >       fPDGCode;
    std::vector< Int_t >       fNMCHitPerMCParticle;
    std::vector< Int_t >       fNUniqueHitsPerMCParticle;
    std::vector< Int_t >       fNMatchedPFParticlesPerMCParticle;
    std::vector< Int_t >       fNMatchedHitsPerMCParticle;
    std::vector< Int_t >       fNHitsBestMatchedPFParticle;
    std::vector< Int_t >       fBestMatchedPFParticlePDGCode;
    std::vector< Int_t >       fBestMatchedPFParticleNHits;
    std::vector< Int_t >       fBestMatchedPFParticleID;
    std::vector< Double_t >    fPx;
    std::vector< Double_t >    fPy;
    std::vector< Double_t >    fPz;
    std::vector< Double_t >    fE;
    std::vector< Double_t >    fEndPx;
    std::vector< Double_t >    fEndPy;
    std::vector< Double_t >    fEndPz;
    std::vector< Double_t >    fEndE;
    std::vector< Double_t >    fVx;
    std::vector< Double_t >    fVy;
    std::vector< Double_t >    fVz;
    std::vector< Double_t >    fT;
    std::vector< Double_t >    fEndx;
    std::vector< Double_t >    fEndy;
    std::vector< Double_t >    fEndz;
    std::vector< Double_t >    fEndT;
    std::vector< Double_t >    fBestMatchedPFParticleVx;
    std::vector< Double_t >    fBestMatchedPFParticleVy;
    std::vector< Double_t >    fBestMatchedPFParticleVz;


    std::vector< std::string > fProcNames;
    // std::vector< std::string > fParentProcNames;

    // Other variables that will be shared between different methods.
    double                                fElectronsToGeV;       // conversion factor

    // std::make_unique< lar::recoana::MCParticleToPFParticle >  fMCParticleToPFParticle;

    // Define our output TTree here
    TTree*    fAnaTree;


}; // class TrackRecoAna

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class implementation

//-----------------------------------------------------------------------
// Constructor
TrackRecoAna::TrackRecoAna( fhicl::ParameterSet const& parameterSet )
    : EDAnalyzer( parameterSet ), fAnaTree( 0 )
{
    // Read in the parameters from the .fcl file.
    this->reconfigure( parameterSet );

    // fMCParticleToPFParticle = std::make_unique< lar::recoana::MCParticleToPFParticle >( fPFParticleProducerLabel, fClusterProducerLabel );
}

//-----------------------------------------------------------------------
// Destructor
TrackRecoAna::~TrackRecoAna()
{}

//-----------------------------------------------------------------------
void TrackRecoAna::beginJob()
{
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us.
    art::ServiceHandle< art::TFileService > tfs;

    // Create our TTree and populate with basic information
    fAnaTree = tfs->make<TTree>("RecoAnalysis", "RecoAna");
    fAnaTree->Branch("Run",                  &fRun,                   "Run/I");
    fAnaTree->Branch("Subrun",               &fSubRun,                "Subrun/I");
    fAnaTree->Branch("Event",                &fEvent,                 "Event/I");
    fAnaTree->Branch("NHits",                &fNHits,                 "NHits/I");
    fAnaTree->Branch("NUniqueHits",          &fNUniqueHits,           "NUniqueHits/I");
    // fAnaTree->Branch("NoiseHits",            &fNoiseHits,             "NoiseHits/I");
    // fAnaTree->Branch("NegativeTrackIDs",     &fNegTrackIds,           "NegativeTrackIDs/I");
    fAnaTree->Branch("NMCParticles",         &fNMCParticles,          "NMCParticles/I");
    fAnaTree->Branch("NPFParticles",         &fNPFParticles,          "NPFParticles/I");
    fAnaTree->Branch("NMatchedPFParticles",  &fNMatchedPFParticles,   "NMatchedPFParticles/I");
    fAnaTree->Branch("RecoPVx",              &fRecoPVx,               "RecoPVx/D");
    fAnaTree->Branch("RecoPVy",              &fRecoPVy,               "RecoPVy/D");
    fAnaTree->Branch("RecoPVz",              &fRecoPVz,               "RecoPVz/D");

    fRun = -1;
    fSubRun = -1;
    fEvent = -1;
    fNHits = 0;
    fNUniqueHits = 0;
    fNMCParticles = 0;
    fNPFParticles = 0;
    fNMatchedPFParticles = 0;
    fRecoPVx = 0.;
    fRecoPVy = 0.;
    fRecoPVz = 0.;

    fPDGCode.resize( fMaxEntries, 0 );
    fNUniqueHitsPerMCParticle.resize( fMaxEntries, 0 );
    fNMatchedPFParticlesPerMCParticle.resize( fMaxEntries, 0 );
    fNMatchedHitsPerMCParticle.resize( fMaxEntries, 0 );
    fNHitsBestMatchedPFParticle.resize( fMaxEntries, 0 );
    fBestMatchedPFParticlePDGCode.resize( fMaxEntries, 0 );
    fBestMatchedPFParticleNHits.resize( fMaxEntries, 0 );
    fBestMatchedPFParticleID.resize( fMaxEntries, -1 );
    fPx.resize( fMaxEntries, 0. );
    fPy.resize( fMaxEntries, 0. );
    fPz.resize( fMaxEntries, 0. );
    fE.resize( fMaxEntries, 0. );
    fEndPx.resize( fMaxEntries, 0. );
    fEndPy.resize( fMaxEntries, 0. );
    fEndPz.resize( fMaxEntries, 0. );
    fEndE.resize( fMaxEntries, 0. );
    fVx.resize( fMaxEntries, 0. );
    fVy.resize( fMaxEntries, 0. );
    fVz.resize( fMaxEntries, 0. );
    fT.resize( fMaxEntries, 0. );
    fEndx.resize( fMaxEntries, 0. );
    fEndy.resize( fMaxEntries, 0. );
    fEndz.resize( fMaxEntries, 0. );
    fEndT.resize( fMaxEntries, 0. );
    fBestMatchedPFParticleVx.resize( fMaxEntries, 0. );
    fBestMatchedPFParticleVy.resize( fMaxEntries, 0. );
    fBestMatchedPFParticleVz.resize( fMaxEntries, 0. );

    fAnaTree->Branch("PDGCode",              fPDGCode.data(),         "PDGCode[NMCParticles]/I");
    fAnaTree->Branch("NUniqueHitsPerMCParticle", fNUniqueHitsPerMCParticle.data(), "NUniqueHitsPerMCParticle[NMCParticles]/I");
    fAnaTree->Branch("NMatchedPFParticlesPerMCParticle", fNMatchedPFParticlesPerMCParticle.data(), "NMatchedPFParticlesPerMCParticle[NMCParticles]/I");
    fAnaTree->Branch("NMatchedHitsPerMCParticle", fNMatchedHitsPerMCParticle.data(), "NMatchedHitsPerMCParticle[NMCParticles]/I");
    fAnaTree->Branch("NHitsBestMatchedPFParticle", fNHitsBestMatchedPFParticle.data(), "NHitsBestMatchedPFParticle[NMCParticles]/I");
    fAnaTree->Branch("BestMatchedPFParticlePDGCode", fBestMatchedPFParticlePDGCode.data(), "BestMatchedPFParticlePDGCode[NMCParticles]/I");
    fAnaTree->Branch("BestMatchedPFParticleNHits", fBestMatchedPFParticleNHits.data(), "BestMatchedPFParticleNHits[NMCParticles]/I");
    fAnaTree->Branch("BestMatchedPFParticleID", fBestMatchedPFParticleID.data(), "BestMatchedPFParticleID[NMCParticles]/I");
    fAnaTree->Branch("Px",  fPx.data(), "Px[NMCParticles]/D");
    fAnaTree->Branch("Py",  fPy.data(), "Py[NMCParticles]/D");
    fAnaTree->Branch("Pz",  fPz.data(), "Pz[NMCParticles]/D");
    fAnaTree->Branch("E",   fE.data(),  "E[NMCParticles]/D");
    fAnaTree->Branch("EndPx",  fEndPx.data(), "EndPx[NMCParticles]/D");
    fAnaTree->Branch("EndPy",  fEndPy.data(), "EndPy[NMCParticles]/D");
    fAnaTree->Branch("EndPz",  fEndPz.data(), "EndPz[NMCParticles]/D");
    fAnaTree->Branch("EndT",   fEndT.data(),  "EndT[NMCParticles]/D");
    fAnaTree->Branch("Vx",  fVx.data(), "Vx[NMCParticles]/D");
    fAnaTree->Branch("Vy",  fVy.data(), "Vy[NMCParticles]/D");
    fAnaTree->Branch("Vz",  fVz.data(), "Vz[NMCParticles]/D");
    fAnaTree->Branch("T",   fT.data(),  "T[NMCParticles]/D");
    fAnaTree->Branch("Endx",  fEndx.data(), "Endx[NMCParticles]/D");
    fAnaTree->Branch("Endy",  fEndy.data(), "Endy[NMCParticles]/D");
    fAnaTree->Branch("Endz",  fEndz.data(), "Endz[NMCParticles]/D");
    fAnaTree->Branch("EndT",  fEndT.data(), "EndT[NMCParticles]/D");
    fAnaTree->Branch("BestMatchedPFParticleVx",  fBestMatchedPFParticleVx.data(), "BestMatchedPFParticleVx[NMCParticles]/D");
    fAnaTree->Branch("BestMatchedPFParticleVy",  fBestMatchedPFParticleVy.data(), "BestMatchedPFParticleVy[NMCParticles]/D");
    fAnaTree->Branch("BestMatchedPFParticleVz",  fBestMatchedPFParticleVz.data(), "BestMatchedPFParticleVz[NMCParticles]/D");


    fProcNames.resize( fMaxEntries, "processname  " );
    // fParentProcNames.resize( fMaxEntries, "processname  " );

    fAnaTree->Branch("MCParticleProc",    &fProcNames);
    // fAnaTree->Branch("MCParticleParentProc",  &fParentProcNames);

}

//-----------------------------------------------------------------------
void TrackRecoAna::beginRun( const art::Run& /*run*/ )
{
    // How to convert from number of electrons to GeV.  The ultimate
    // source of this conversion factor is
    // ${LARSIM_DIR}/include/SimpleTypesAndConstants/PhysicalConstants.h.
//    art::ServiceHandle<sim::LArG4Parameters> larParameters;
//    fElectronsToGeV = 1./larParameters->GeVToElectrons();
}

//-----------------------------------------------------------------------
void TrackRecoAna::reconfigure( fhicl::ParameterSet const& p )
{
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fSimulationProducerLabel    = p.get< std::string >( "SimulationLabel" );
    fMCHitCollectionModuleLabel = p.get< std::string >( "MCHitFinderLabel" );
    fPFParticleProducerLabel    = p.get< std::string >( "PFParticleLabel" );
    fHitProducerLabel           = p.get< std::string >( "HitLabel" );
    fClusterProducerLabel       = p.get< std::string >( "ClusterLabel" );
    // fTrackProducerLabel         = p.get< std::string >( "TrackProducerLabel" );
    fMaxEntries                 = p.get< int         >( "MaxEntries", 300 );

    return;
}

//-----------------------------------------------------------------------
void TrackRecoAna::PrepareEvent( const art::Event &evt )
{
    fRun                   = evt.run();
    fEvent                 = evt.id().event();
    fSubRun                = evt.subRun();
    fNHits                 = 0;
    fNUniqueHits           = 0;
    // fNNoiseHits            = 0;
    fNMCParticles          = 0;
    fNPFParticles          = 0;
    fNMatchedPFParticles   = 0;
    fRecoPVx               = 0.;
    fRecoPVy               = 0.;
    fRecoPVz               = 0.;

    fPDGCode.assign( fMaxEntries, 0 );
    fNUniqueHitsPerMCParticle.assign( fMaxEntries, 0 );
    fNMatchedPFParticlesPerMCParticle.assign( fMaxEntries, 0 );
    fNMatchedHitsPerMCParticle.assign( fMaxEntries, 0 );
    fNHitsBestMatchedPFParticle.assign( fMaxEntries, 0 );
    fBestMatchedPFParticlePDGCode.assign( fMaxEntries, 0 );
    fBestMatchedPFParticleNHits.assign( fMaxEntries, 0 );
    fBestMatchedPFParticleID.assign( fMaxEntries, -1 );
    fPx.resize( fMaxEntries, 0. );
    fPy.resize( fMaxEntries, 0. );
    fPz.resize( fMaxEntries, 0. );
    fE.resize( fMaxEntries, 0. );
    fEndPx.resize( fMaxEntries, 0. );
    fEndPy.resize( fMaxEntries, 0. );
    fEndPz.resize( fMaxEntries, 0. );
    fEndE.resize( fMaxEntries, 0. );
    fVx.resize( fMaxEntries, 0. );
    fVy.resize( fMaxEntries, 0. );
    fVz.resize( fMaxEntries, 0. );
    fT.resize( fMaxEntries, 0. );
    fEndx.resize( fMaxEntries, 0. );
    fEndy.resize( fMaxEntries, 0. );
    fEndz.resize( fMaxEntries, 0. );
    fEndT.resize( fMaxEntries, 0. );
    fBestMatchedPFParticleVx.resize( fMaxEntries, 0. );
    fBestMatchedPFParticleVy.resize( fMaxEntries, 0. );
    fBestMatchedPFParticleVz.resize( fMaxEntries, 0. );
    fProcNames.assign( fMaxEntries, "processname  " );
    // fParentProcNames.assign( fMaxEntries, "processname  " );
}

//-----------------------------------------------------------------------
void TrackRecoAna::analyze( const art::Event& event )
{
    this->PrepareEvent( event );

    // The first step is to attempt to recover the collection of MCHits that
    // we will need for doing our PFParticle to MC matching
    art::Handle< std::vector< sim::MCHitCollection > > mcHitCollectionHandle;
    event.getByLabel( fMCHitCollectionModuleLabel, mcHitCollectionHandle );

    if ( !mcHitCollectionHandle.isValid() ) {
        mf::LogWarning("TrackRecoAna") << "===>> NO MCHitColllection found.";
        fAnaTree->Fill();
        return;
    }

    // Recover this into a local stl version
    const std::vector< sim::MCHitCollection > &MCHits = *mcHitCollectionHandle;

    // Recover the PFParticles, the main products for our next major loop
    art::Handle< std::vector< recob::PFParticle > > pfParticleHandle;
    event.getByLabel( fPFParticleProducerLabel, pfParticleHandle );

    // no point continuing if no PFParticles
    if ( !pfParticleHandle.isValid() ) {
        mf::LogWarning("TrackRecoAna") << "===>> NO PFParticle collection found.";
        fAnaTree->Fill();
        return;
    }

    if ( pfParticleHandle->empty() ) {
        mf::LogWarning("TrackRecoAna") << "===>> PFParticle collection is empty.";
        fAnaTree->Fill();
        return;
    }

    // Recover the MC particles
    art::Handle< std::vector< simb::MCParticle > > particleHandle;
    event.getByLabel( fSimulationProducerLabel, particleHandle );

    if ( !particleHandle.isValid() ) {
        mf::LogWarning("TrackRecoAna") << "===>> NO MCParticle collection found.";
        fAnaTree->Fill();
        return;
    }

    // Let's recover the MCTruth objects
    // art::FindOneP< simb::MCTruth > mcTruthAssns( particleHandle, event, fSimulationProducerLabel );

    // The MCParticle objects are not necessarily in any particular
    // order. Since we may have to search the list of particles, let's
    // put them into a sorted map that will make searching fast and
    // easy. To save both space and time, the map will not contain a
    // copy of the MCParticle, but a pointer to it.
    std::map< int, const simb::MCParticle* > particleMap;

    // Before starting to loop through the particles, we are going to want to
    // build a mapping between hits and track id's
    // Start by recovering info from the event store
    art::Handle< std::vector< recob::Hit > > hitHandle;
    event.getByLabel( fHitProducerLabel, hitHandle );

    if ( !hitHandle.isValid() ) {
        mf::LogWarning("TrackRecoAna") << "===>> No hit collection found.";
        fAnaTree->Fill();
        return;
    }

    auto MCParticleToPFParticle = std::make_unique< lar::recoana::MCParticleToPFParticle >( fPFParticleProducerLabel, fClusterProducerLabel );
    // our ultimate goal here are the following two maps:
    // 1) a map between a recob::Hit and the track ID (meaning MCParticle)
    // 2) a map between the trackID (MCParticle) and the hits it created
    lar::recoana::MCParticleToPFParticle::HitToParticleMap HitToParticle;
    lar::recoana::MCParticleToPFParticle::ParticleToHitMap ParticleToHit;

    // Recover the list of hits in stl vector form
    const std::vector< recob::Hit >& Hits = *hitHandle;

    MCParticleToPFParticle->setup( *( lar::providerFrom< detinfo::DetectorClocksService >() ),
                                   *( lar::providerFrom<geo::Geometry>() ) );

    // Build out our MCParticle <---> reco Hit maps
    MCParticleToPFParticle->MakeHitParticleMaps( MCHits, Hits, HitToParticle, ParticleToHit );

    mf::LogDebug("TrackRecoAna") << "Size of Hits         : " << Hits.size();
    mf::LogDebug("TrackRecoAna") << "Size of MCParticles  : " << particleHandle->size();
    mf::LogDebug("TrackRecoAna") << "Size of PFParticles  : " << pfParticleHandle->size();
    mf::LogDebug("TrackRecoAna") << "Size of HitToParticle: " << HitToParticle.size();
    mf::LogDebug("TrackRecoAna") << "Size of ParticleToHit: " << ParticleToHit.size();

    // Now define the maps relating pfparticles to tracks
    lar::recoana::MCParticleToPFParticle::TrackIDToPFParticleVecMap TrackIDToPFParticle;
    lar::recoana::MCParticleToPFParticle::PFParticleToTrackHit2DMap PFParticleToTrackHit;

    // Something to keep track of number of hits associated to a cluster
    std::map< art::Ptr< recob::PFParticle >, int > PFParticleToHitCnt;

    lar::recoana::MCParticleToPFParticle::TrackIDToMatchedPFParticleHitCntMap TrackIDToMatchedPFParticleHitCnt;

    // Call the method for building out the PFParticle data structures
    MCParticleToPFParticle->MakePFParticleMaps( event, pfParticleHandle, HitToParticle, 
                                                PFParticleToTrackHit, TrackIDToPFParticle, PFParticleToHitCnt,
                                                TrackIDToMatchedPFParticleHitCnt );

    // Recover the association from the PFParticle to the vertex
    art::FindManyP< recob::Vertex > pfParticleVertexAssns( pfParticleHandle, event, fPFParticleProducerLabel );

    // Loop over MC particles
    for ( size_t iMCPart = 0; iMCPart < particleHandle->size(); ++iMCPart ) {

        art::Ptr< simb::MCParticle > particle( particleHandle, iMCPart );

        // Process of the MC particle
        auto Process = particle->Process();

        // At the moment we only care of primary particles for the events generated by a particle gun
        if ( Process != "primary" ) continue;
        

        fProcNames[fNMCParticles] = Process;
        // Recover the track ID, for historical reasons we call it "best"
        int ParticleTrackID = particle->TrackId();
        mf::LogDebug("TrackRecoAna") << "Process[" << fNMCParticles << "]: " << fProcNames[fNMCParticles] 
                                     << ": TrackID: " << ParticleTrackID << ": Mother: " << particle->Mother();

        // Did this mc particle leave hits in the TPC?
        lar::recoana::MCParticleToPFParticle::ParticleToHitMap::iterator ParticleToHitItr = ParticleToHit.find( ParticleTrackID );
        lar::recoana::MCParticleToPFParticle::TrackIDToMatchedPFParticleHitCntMap::iterator TrackIDToMatchedPFParticleHitCntItr = TrackIDToMatchedPFParticleHitCnt.find( ParticleTrackID );

        // Let's get the total number of reconstructed hits that are created by this MCParticle
        // Count number of hits in each view
        if ( ParticleToHitItr != ParticleToHit.end() ) {
            fNMatchedHitsPerMCParticle[fNMCParticles] = ParticleToHitItr->second.size();
            mf::LogDebug("TrackRecoAna") << "NMatchedHitsPerMCParticle[" << fNMCParticles << "]: " << fNMatchedHitsPerMCParticle[fNMCParticles];

            for ( const auto& hit : ParticleToHitItr->second ) {
                if ( HitToParticle[hit].size() < 2 ) fNUniqueHitsPerMCParticle[fNMCParticles]++;
            }
        }


        // How many PFParticles does the MCParticle match to?
        lar::recoana::MCParticleToPFParticle::TrackIDToPFParticleVecMap::iterator TrackIDToPFParticleItr = TrackIDToPFParticle.find( ParticleTrackID );
        if ( TrackIDToPFParticleItr != TrackIDToPFParticle.end() ) {
            fNMatchedPFParticlesPerMCParticle[fNMCParticles] = TrackIDToPFParticleItr->second.size();

            // Look for the best matched PFParticle by counting hits
            // Since we have sorted the PFParticles in the TrackIDToPFParticle map, the first PFParticle is the one
            if ( TrackIDToPFParticleItr->second.empty() ) continue;
            const auto& pfpart = TrackIDToPFParticleItr->second.at(0);
            fNHitsBestMatchedPFParticle[fNMCParticles] = PFParticleToHitCnt[pfpart];
            fBestMatchedPFParticlePDGCode[fNMCParticles] = pfpart->PdgCode();
            if ( TrackIDToMatchedPFParticleHitCntItr != TrackIDToMatchedPFParticleHitCnt.end() ) {
                if ( TrackIDToMatchedPFParticleHitCntItr->second.size() > 0 )
                    fBestMatchedPFParticleNHits[fNMCParticles] = TrackIDToMatchedPFParticleHitCntItr->second.at(0);
                else fBestMatchedPFParticleNHits[fNMCParticles] = 0;
            }
            fBestMatchedPFParticleID[fNMCParticles] = pfpart.key();

            // Fill the reconstructed vertex associated to the PFParticle
            double xyz[3];
            std::vector< art::Ptr< recob::Vertex > > const& vertices = pfParticleVertexAssns.at( pfpart.key() );
            if ( vertices.empty() ) continue;
            art::Ptr< recob::Vertex > vertex = vertices.at(0);
            vertex->XYZ( xyz );
            fBestMatchedPFParticleVx[fNMCParticles] = xyz[0];
            fBestMatchedPFParticleVy[fNMCParticles] = xyz[1];
            fBestMatchedPFParticleVz[fNMCParticles] = xyz[2];
        }

        // It is pointless to go through the rest of the loop if there are no hits to analyze
        // if ( nTrueMCHits < 1 ) continue;

        // Recover the particle's identity
        int TrackPDGCode = particle->PdgCode();
        fPDGCode[fNMCParticles] = TrackPDGCode;        
        mf::LogDebug("TrackRecoAna") << "fPDGCode[" << fNMCParticles << "]: " << fPDGCode[fNMCParticles];

        fPx[fNMCParticles] = particle->Px();
        fPy[fNMCParticles] = particle->Py();
        fPz[fNMCParticles] = particle->Pz();
        fE[fNMCParticles]  = particle->E();
        fEndPx[fNMCParticles] = particle->EndPx();
        fEndPy[fNMCParticles] = particle->EndPy();
        fEndPz[fNMCParticles] = particle->EndPz();
        fEndE[fNMCParticles]  = particle->EndE();
        fVx[fNMCParticles] = particle->Vx();
        fVy[fNMCParticles] = particle->Vy();
        fVz[fNMCParticles] = particle->Vz();
        fT[fNMCParticles] = particle->T();
        fEndx[fNMCParticles] = particle->EndX();
        fEndy[fNMCParticles] = particle->EndY();
        fEndz[fNMCParticles] = particle->EndZ();
        fEndT[fNMCParticles] = particle->EndT();

        fNHits += fNMatchedHitsPerMCParticle[fNMCParticles];
        fNUniqueHits += fNUniqueHitsPerMCParticle[fNMCParticles];
        // fNNoiseHits +=;
        fNMatchedPFParticles += fNMatchedPFParticlesPerMCParticle[fNMCParticles];
        fNMCParticles++;

        if ( fNMCParticles > fMaxEntries ) {
            mf::LogWarning("TrackRecoAna") << "===>> Exceeding particle with hits count: " << fNMCParticles;
            break;
        }
    } // end of loop over MC particles

    fNPFParticles = pfParticleHandle->size();

    // Looking for the primary vertex.
    // Currently we look for a neutrino PFParticle which has the most daughter PFParticles
    int nDaughtersMax = 0;
    size_t iNuPFParticle = 0;
    for ( size_t iPFPart = 0; iPFPart < pfParticleHandle->size(); ++iPFPart ) {
        recob::PFParticle const& pfpart = (*pfParticleHandle)[iPFPart];
        int PDGCode = pfpart.PdgCode();
        if ( PDGCode != 12 && PDGCode != 14 ) continue;
        mf::LogDebug("TrackRecoAna") << "PDGCode of PFParticle " << iPFPart << ": " << PDGCode;
        int nDaughters = pfpart.NumDaughters();
        if ( nDaughters > nDaughtersMax ) {
            nDaughtersMax = nDaughters;
            iNuPFParticle = iPFPart;
        }
    }
    std::vector< art::Ptr< recob::Vertex > > const& PVtxs = pfParticleVertexAssns.at( iNuPFParticle );
    if ( !PVtxs.empty() ) {
        art::Ptr< recob::Vertex > PVtx = PVtxs.at(0);
        double pvxyz[3];
        PVtx->XYZ( pvxyz );
        fRecoPVx = pvxyz[0];
        fRecoPVy = pvxyz[1];
        fRecoPVz = pvxyz[2];
    }
    fAnaTree->Fill();

    return;
}
// ----------------------------------------------------------------------------------------------------
DEFINE_ART_MODULE(TrackRecoAna)

} // namespace TrackRecoAna

#endif
