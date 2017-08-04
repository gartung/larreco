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
    void PrepareEvent( const art::Event& evt, int numColumns );

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
    Int_t    fNMathcedPFParticles;

    // MCParticle-wise information
    std::vector< Int_t >       fPDGCode;
    std::vector< Int_t >       fNMCHitPerMCParticle;
    std::vector< Int_t >       fNUniqueHitsPerMCParticle;
    std::vector< Int_t >       fNMatchedPFParticlesPerMCParticle;
    std::vector< Int_t >       fNMatchedHitsPerMCParticle;
    std::vector< Int_t >       fNHitsBestMatchedPFParticle;
    std::vector< Int_t >       fBestMatchedPFParticlePDGCode;

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
    fAnaTree->Branch("NMathcedPFParticles",  &fNMathcedPFParticles,   "NMathcedPFParticles/I");

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
    fSimulationProducerLabel    = p.get< std::string >("SimulationLabel");
    fMCHitCollectionModuleLabel = p.get< std::string >("McHitFinderLabel");
    fPFParticleProducerLabel    = p.get< std::string >("PFParticleLabel");
    fHitProducerLabel           = p.get< std::string >("HitLabel");
    fClusterProducerLabel       = p.get< std::string >("ClusterLabel");
    // fTrackProducerLabel         = p.get< std::string >("TrackProducerLabel");
    fMaxEntries                 = p.get< int         >("MaxEntries", 300);

    return;
}

//-----------------------------------------------------------------------
void TrackRecoAna::PrepareEvent( const art::Event &evt, int numColumns )
{
    fRun                   = evt.run();
    fEvent                 = evt.id().event();
    fSubRun                = evt.subRun();
    fNHits                 = 0;
    fNUniqueHits           = 0;
    // fNNoiseHits            = 0;
    fNMCParticles          = 0;
    fNPFParticles          = 0;
    fNMathcedPFParticles   = 0;

    // Set the number of entries for each vector
    // But no more than the maximum requested when we started
    size_t maxEntries = numColumns < fMaxEntries ? numColumns : fMaxEntries;

    fPDGCode.assign( maxEntries, 0 );
    fNUniqueHitsPerMCParticle.assign( maxEntries, 0 );
    fNMatchedPFParticlesPerMCParticle.assign( maxEntries, 0 );
    fNMatchedHitsPerMCParticle.assign( maxEntries, 0 );
    fNHitsBestMatchedPFParticle.assign( maxEntries, 0 );
    fBestMatchedPFParticlePDGCode.assign( maxEntries, 0 );
}

//-----------------------------------------------------------------------
void TrackRecoAna::analyze( const art::Event& event )
{

    // The first step is to attempt to recover the collection of MCHits that
    // we will need for doing our PFParticle to MC matching
    art::Handle< std::vector< sim::MCHitCollection > > mcHitCollectionHandle;
    event.getByLabel( fMCHitCollectionModuleLabel, mcHitCollectionHandle );

    if ( !mcHitCollectionHandle.isValid() ) {
        mf::LogDebug("TrackRecoAna") << "===>> NO McHitColllection found.";
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
        mf::LogDebug("TrackRecoAna") << "===>> NO PFParticle collection found.";
        fAnaTree->Fill();
        return;
    }

    // Recover the MC particles
    art::Handle< std::vector< simb::MCParticle > > particleHandle;
    event.getByLabel( fSimulationProducerLabel, particleHandle );

    if ( !particleHandle.isValid() ) {
        mf::LogDebug("TrackRecoAna") << "===>> NO MCParticle collection found.";
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
        mf::LogDebug("TrackRecoAna") << "===>> No hit collection found.";
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

    // Now define the maps relating pfparticles to tracks
    lar::recoana::MCParticleToPFParticle::TrackIDToPFParticleVecMap TrackIDToPFParticle;
    lar::recoana::MCParticleToPFParticle::PFParticleToTrackHit2DMap PFParticleToTrackHit;

    // Something to keep track of number of hits associated to a cluster
    std::map< art::Ptr< recob::PFParticle >, int > PFParticleToHitCnt;

    // Call the method for building out the PFParticle data structures
    MCParticleToPFParticle->MakePFParticleMaps( event, pfParticleHandle, HitToParticle, 
                                                PFParticleToTrackHit, TrackIDToPFParticle, PFParticleToHitCnt );

    this->PrepareEvent( event, ParticleToHit.size() );


    // Loop over MC particles
    for ( size_t iMCPart = 0; iMCPart < particleHandle->size(); ++iMCPart ) {

        art::Ptr< simb::MCParticle > particle( particleHandle, iMCPart );

        // Recover the track ID, for historical reasons we call it "best"
        int ParticleTrackID = particle->TrackId();

        // Did this mc particle leave hits in the TPC?
        lar::recoana::MCParticleToPFParticle::ParticleToHitMap::iterator ParticleToHitItr = ParticleToHit.find( ParticleTrackID );

        // Let's get the total number of reconstructed hits that are created by this MCParticle
        // Count number of hits in each view
        if ( ParticleToHitItr != ParticleToHit.end() ) {
            fNMatchedHitsPerMCParticle[iMCPart] = ParticleToHitItr->second.size();

            for ( const auto& hit : ParticleToHitItr->second ) {
                if ( HitToParticle[hit].size() < 2 ) fNUniqueHitsPerMCParticle[iMCPart]++;
            }
        }


        // How many PFParticles does the MCParticle match to?
        lar::recoana::MCParticleToPFParticle::TrackIDToPFParticleVecMap::iterator TrackIDToPFParticleItr = TrackIDToPFParticle.find( ParticleTrackID );
        if ( TrackIDToPFParticleItr != TrackIDToPFParticle.end() ) {
            fNMatchedPFParticlesPerMCParticle[iMCPart] = TrackIDToPFParticleItr->second.size();

            // Look for the best matched PFParticle by counting hits
            // Since we have sorted the PFParticles in the TrackIDToPFParticle map, the first PFParticle is the one
            const auto& pfpart = TrackIDToPFParticleItr->second.at(0);
            fNHitsBestMatchedPFParticle[iMCPart] = PFParticleToTrackHit[pfpart].size();
            fBestMatchedPFParticlePDGCode[iMCPart] = pfpart->PdgCode();
        }

        // It is pointless to go through the rest of the loop if there are no hits to analyze
        // if ( nTrueMCHits < 1 ) continue;

        // Recover the particle's identity
        int TrackPDGCode = particle->PdgCode();
        fPDGCode[iMCPart] = TrackPDGCode;        

        fNHits += fNMatchedHitsPerMCParticle[iMCPart];
        fNUniqueHits += fNUniqueHitsPerMCParticle[iMCPart];
        // fNNoiseHits +=;
        fNMCParticles += 1;
        fNMathcedPFParticles += fNMatchedPFParticlesPerMCParticle[iMCPart];

    } // end of loop over MC particles

    fNPFParticles = pfParticleHandle->size();
    fAnaTree->Fill();

    return;
}
// ----------------------------------------------------------------------------------------------------
DEFINE_ART_MODULE(TrackRecoAna)

} // namespace TrackRecoAna

#endif
