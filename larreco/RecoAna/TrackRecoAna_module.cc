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

    void PrepareParticle( int numColumns );

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

    // MCParticle-wise information
    std::vector< Int_t >       fPDGCode;
    std::vector< Int_t >       fNMCHitPerMCParticle;
    std::vector< Int_t >       fNUniqueHitsPerMCParticle;
    std::vector< Int_t >       fNMatchedPFParticlesPerMCParticle;
    std::vector< Int_t >       fNMatchedHitsPerMCParticle;
    std::vector< Int_t >       fNHitsBestMatchedPFParticle;
    std::vector< Int_t >       fBestMatchedPFParticlePDGCode;

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

    fRun = -1;
    fSubRun = -1;
    fEvent = -1;
    fNHits = 0;
    fNUniqueHits = 0;
    fNMCParticles = 0;
    fNPFParticles = 0;
    fNMatchedPFParticles = 0;

    fPDGCode.resize( fMaxEntries, 0 );
    fNUniqueHitsPerMCParticle.resize( fMaxEntries, 0 );
    fNMatchedPFParticlesPerMCParticle.resize( fMaxEntries, 0 );
    fNMatchedHitsPerMCParticle.resize( fMaxEntries, 0 );
    fNHitsBestMatchedPFParticle.resize( fMaxEntries, 0 );
    fBestMatchedPFParticlePDGCode.resize( fMaxEntries, 0 );

    fAnaTree->Branch("PDGCode",              fPDGCode.data(),         "PDGCode[NMCParticles]/I");
    fAnaTree->Branch("NUniqueHitsPerMCParticle", fNUniqueHitsPerMCParticle.data(), "NUniqueHitsPerMCParticle[NMCParticles]/I");
    fAnaTree->Branch("NMatchedPFParticlesPerMCParticle", fNMatchedPFParticlesPerMCParticle.data(), "NMatchedPFParticlesPerMCParticle[NMCParticles]/I");
    fAnaTree->Branch("NMatchedHitsPerMCParticle", fNMatchedHitsPerMCParticle.data(), "NMatchedHitsPerMCParticle[NMCParticles]/I");
    fAnaTree->Branch("NHitsBestMatchedPFParticle", fNHitsBestMatchedPFParticle.data(), "NHitsBestMatchedPFParticle[NMCParticles]/I");
    fAnaTree->Branch("BestMatchedPFParticlePDGCode", fBestMatchedPFParticlePDGCode.data(), "BestMatchedPFParticlePDGCode[NMCParticles]/I");

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
}

//-----------------------------------------------------------------------
void TrackRecoAna::PrepareParticle( int numColumns )
{
    // Set the number of entries for each vector
    // But no more than the maximum requested when we started
    size_t maxEntries = numColumns < fMaxEntries ? numColumns : fMaxEntries;

    fPDGCode.assign( maxEntries, 0 );
    fNUniqueHitsPerMCParticle.assign( maxEntries, 0 );
    fNMatchedPFParticlesPerMCParticle.assign( maxEntries, 0 );
    fNMatchedHitsPerMCParticle.assign( maxEntries, 0 );
    fNHitsBestMatchedPFParticle.assign( maxEntries, 0 );
    fBestMatchedPFParticlePDGCode.assign( maxEntries, 0 );
    fProcNames.assign( maxEntries, "processname  " );
    // fParentProcNames.assign( maxEntries, "processname  " );
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

    // Call the method for building out the PFParticle data structures
    MCParticleToPFParticle->MakePFParticleMaps( event, pfParticleHandle, HitToParticle, 
                                                PFParticleToTrackHit, TrackIDToPFParticle, PFParticleToHitCnt );

    this->PrepareParticle( ParticleToHit.size() );

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
            const auto& pfpart = TrackIDToPFParticleItr->second.at(0);
            fNHitsBestMatchedPFParticle[fNMCParticles] = PFParticleToTrackHit[pfpart].size();
            fBestMatchedPFParticlePDGCode[fNMCParticles] = pfpart->PdgCode();
        }

        // It is pointless to go through the rest of the loop if there are no hits to analyze
        // if ( nTrueMCHits < 1 ) continue;

        // Recover the particle's identity
        int TrackPDGCode = particle->PdgCode();
        fPDGCode[fNMCParticles] = TrackPDGCode;        
        mf::LogDebug("TrackRecoAna") << "fPDGCode[" << fNMCParticles << "]: " << fPDGCode[fNMCParticles];

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
    fAnaTree->Fill();

    return;
}
// ----------------------------------------------------------------------------------------------------
DEFINE_ART_MODULE(TrackRecoAna)

} // namespace TrackRecoAna

#endif
