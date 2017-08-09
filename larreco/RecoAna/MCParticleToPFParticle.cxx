/**
 * @file   MCParticleToPFParticle.cxx
 * @brief  Algorithms to match reconstructed hits to track IDs in MC particles (through MC hits),
 *         and to match track IDs to PFParticles
 * @author Yun-Tse Tsai  (stealing codes from Tracy Usher)
 * @date   August 3rd, 2017
 *
 */

#ifndef LARRECO_ANA_MCPARTICLETOPFPARTICLE_CXX
#define LARRECO_ANA_MCPARTICLETOPFPARTICLE_CXX

// LArSoft
#include "larreco/RecoAna/MCParticleToPFParticle.h"

// Framework includes
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include <iterator>

//----------------------------------------------------------------------------
// Build map of reconstructed 2D hits <--> MC particles (track IDs)
//----------------------------------------------------------------------------
void lar::recoana::MCParticleToPFParticle::MakeHitParticleMaps( const std::vector< sim::MCHitCollection >& MCHits,
                                                                const std::vector< recob::Hit >&           Hits,
                                                                HitToParticleMap&                          HitToParticle,
                                                                ParticleToHitMap&                          ParticleToHit )
{

    int nBigHits = 0;
    int nBigMatchedHits = 0;
    // Ok, so this loop obviously takes the MC information and builds two maps
    // 1) a map from a Hit2D object to the track ID's that made it
    // 2) the reverse map, going from track ID to Hit2D object
    for ( const recob::Hit& hit : Hits ) {
        if ( hit.PeakAmplitude() > 5 ) nBigHits++;
        const geo::WireID& wireId = hit.WireID();

        unsigned int channel = geom->PlaneWireToChannel( wireId.Plane, wireId.Wire, wireId.TPC, wireId.Cryostat );

        const std::vector< sim::MCHit >& MCHitsChannel = MCHits[channel];

        int start_tdc    = detc->TPCTick2TDC( hit.StartTick() );
        int end_tdc      = detc->TPCTick2TDC( hit.EndTick()   );
        int hitStart_tdc = detc->TPCTick2TDC( hit.PeakTime() - 3.*hit.SigmaPeakTime() );
        int hitEnd_tdc   = detc->TPCTick2TDC( hit.PeakTime() + 3.*hit.SigmaPeakTime() );

        start_tdc = std::max(start_tdc, hitStart_tdc);
        end_tdc   = std::min(end_tdc,   hitEnd_tdc  );

        sim::MCHit startTime;
        sim::MCHit endTime;

        startTime.SetTime(start_tdc, 0);
        endTime.SetTime(end_tdc, 0);

        std::vector< sim::MCHit >::const_iterator startItr = std::lower_bound( MCHitsChannel.begin(), MCHitsChannel.end(), startTime );
        std::vector< sim::MCHit >::const_iterator endItr   = std::upper_bound( startItr,              MCHitsChannel.end(), endTime );

        if ( startItr != MCHitsChannel.end() ) {
            double totalCharge = 0.;
            for ( auto tmpItr = startItr; tmpItr != endItr; ++tmpItr )  totalCharge += (*tmpItr).Charge();

            if ( !( totalCharge > 0. ) ) continue;
            if ( hit.PeakAmplitude() > 5 ) nBigMatchedHits++;

            while ( startItr != endItr ) {
                int trackID = (*startItr).PartTrackId();
                // mf::LogDebug("MCParticleToPFParticle::MakeHitParticleMaps") << "TrackID: " << trackID;
                if ( trackID < 0 ) {
                    trackID = -trackID;
                    // fNegTrackIds++;
                }

                HitToParticle[&hit].insert( trackID );
                ParticleToHit[trackID].insert( &hit );
                startItr++;
            }
        } // else {
          //   fNoiseHits++;
        // } // if ( startItr != MCHitsChannel.end() )
    }

    // This is for debugging
    mf::LogDebug("MCParticleToPFParticle::MakeHitParticleMaps") << "# big hits          : " << nBigHits;
    mf::LogDebug("MCParticleToPFParticle::MakeHitParticleMaps") << "# big matched hits  : " << nBigMatchedHits;
    /****
    unsigned int theChannel = 0;
    bool isTheChannel = true;
    int nRecoHits = 0;
    for ( const recob::Hit& hit : Hits ) {
        const geo::WireID& wireId = hit.WireID();
        unsigned int channel = geom->PlaneWireToChannel( wireId.Plane, wireId.Wire, wireId.TPC, wireId.Cryostat );
        if ( isTheChannel && channel > 5500 ) {
            theChannel = channel;
            isTheChannel = false;
        }

        if ( channel != theChannel ) continue;

        nRecoHits++;
        mf::LogDebug("MCParticleToPFParticle::MakeHitParticleMaps") << "Reco'ed hit channel: " << channel;
        int start_tdc    = detc->TPCTick2TDC( hit.StartTick() );
        int end_tdc      = detc->TPCTick2TDC( hit.EndTick()   );
        int hitStart_tdc = detc->TPCTick2TDC( hit.PeakTime() - 3.*hit.SigmaPeakTime() );
        int hitEnd_tdc   = detc->TPCTick2TDC( hit.PeakTime() + 3.*hit.SigmaPeakTime() );

        start_tdc = std::max(start_tdc, hitStart_tdc);
        end_tdc   = std::min(end_tdc,   hitEnd_tdc  );
        mf::LogDebug("MCParticleToPFParticle::MakeHitParticleMaps") << "Reco'ed hit start tdc: " << start_tdc;
        mf::LogDebug("MCParticleToPFParticle::MakeHitParticleMaps") << "Reco'ed hit end tdc  : " << end_tdc;
        const auto matchedParticles = HitToParticle[&hit];
        // std::set< int >::iterator it;
        mf::LogDebug("MCParticleToPFParticle::MakeHitParticleMaps") << "# matched particle   : " << matchedParticles.size();
        // mf::LogDebug("MCParticleToPFParticle::MakeHitParticleMaps") << "Matched trackIDs     : ";
        // for ( it = matchedParticles.begin(); it != matchedParticles.end(); ++it )
        //     mf::LogDebug("MCParticleToPFParticle::MakeHitParticleMaps") << (*it);
    }

    const std::vector< sim::MCHit >& MCHitsChannel = MCHits[theChannel];
    int nMCHits = MCHitsChannel.size();
    for ( const sim::MCHit& mcHit : MCHitsChannel ) {
        mf::LogDebug("MCParticleToPFParticle::MakeHitParticleMaps") << "MC hit peak time     : " << mcHit.PeakTime();
        mf::LogDebug("MCParticleToPFParticle::MakeHitParticleMaps") << "MC hit trackIDs      : " << mcHit.PartTrackId();
    }

    mf::LogDebug("MCParticleToPFParticle::MakeHitParticleMaps") << "# reco'ed hits in channel " << theChannel << ": " << nRecoHits;
    mf::LogDebug("MCParticleToPFParticle::MakeHitParticleMaps") << "# mc      hits in channel " << theChannel << ": " << nMCHits;
    ***/
    return;
}

//----------------------------------------------------------------------------
// Build maps for PFParticles
//----------------------------------------------------------------------------
void lar::recoana::MCParticleToPFParticle::MakePFParticleMaps( const art::Event&                                      event,
                                                               const art::Handle< std::vector< recob::PFParticle > >& pfParticleHandle,
                                                               const HitToParticleMap&                                HitToParticle,
                                                               PFParticleToTrackHit2DMap&                             PFParticleToTrackHit2D,
                                                               TrackIDToPFParticleVecMap&                             TrackIDToPFParticles,
                                                               PFParticleHitCntMap&                                   PFParticleToHitCnt) {

    // Recover pfparticle to cluster associations
    art::FindManyP< recob::Cluster > pfParticleClusterAssns( pfParticleHandle, event, fPFParticleProducerLabel );

    // We need a handle to the collection of clusters in the data store so we can
    // handle associations to hits.
    art::Handle< std::vector< recob::Cluster > > clusterHandle;
    event.getByLabel( fClusterProducerLabel, clusterHandle );

    // Now recover the collection of associations to hits
    art::FindManyP< recob::Hit > clusterHitAssns( clusterHandle, event, fClusterProducerLabel );

    // Commence looping over PFParticles
    for ( size_t iPFPart = 0; iPFPart < pfParticleHandle->size(); ++iPFPart ) {
        // Get our PFParticle
        art::Ptr< recob::PFParticle > pfParticle( pfParticleHandle, iPFPart );

        // Get the list of clusters associated to this PFParticle
        std::vector< art::Ptr< recob::Cluster > > clusters = pfParticleClusterAssns.at( pfParticle.key() );

        // Number 2D hits this PFParticle
        int nPFParticleHits(0);

        // To categorize the PFParticle (associate to MCParticle), we will
        // create and fill an instance of a TrackIDToHitMap.
        // So, for a given PFParticle we will have a map of the track ID's (MCParticles)
        // and the hits that they contributed energy too
        PFParticleToTrackHit2D[pfParticle] = TrackIDToHit2DMap();
        TrackIDToHit2DMap& trackIdToHit2DMap  = PFParticleToTrackHit2D[pfParticle];

        // Loop over the associated clusters
        for ( const auto& cluster : clusters ) {
            std::vector< art::Ptr< recob::Hit > > hits;

            try {
                hits = clusterHitAssns.at( cluster->ID() );
            } catch(...) {
                continue;
            }

            // Keep track of this to hopefully save some time later
            PFParticleToHitCnt[pfParticle] += hits.size();

            // Something to count MCParticles contributing here
            std::set< int > trackIdCntSet;
            int             nMultiParticleHits(0);
            int             nNoParticleHits(0);

            nPFParticleHits += hits.size();

            // To fill out this map we now loop over the recovered Hit2D's and stuff into the map
            for ( const auto& hit : hits ) {
                // Given the Hit2D, recover the list of asscociated track ID's (MCParticles)
                HitToParticleMap::const_iterator HitToParticleMapItr = HitToParticle.find( hit.get() );

                if ( HitToParticleMapItr != HitToParticle.end() ) {
                    // More than one MCParticle can contribute energy to make a given hit
                    // So loop over the track ID's for this hit
                    for ( const auto& trackId : HitToParticleMapItr->second ) {
                        // mf::LogDebug("MCParticleToPFParticle::MakePFParticleMaps") << "Hit matched trackIDs: " << trackId;
                        trackIdToHit2DMap[trackId].insert( hit );
                        trackIdCntSet.insert( trackId );
                    }

                    if ( HitToParticleMapItr->second.size() > 1) nMultiParticleHits++;
                } else nNoParticleHits++;  // end of finding the hit in the HitToParticle map
            } // end loop over hits

            /***
            mf::LogDebug("MCParticleToPFParticle::MakePFParticleMaps") << "PFParticle " << iPFPart 
                            << ": PDGCode: " << pfParticle->PdgCode()
                            << ": Cluster: " << cluster->View() << ": nHits: " << hits.size()
                            << ": nPFParticleHits: " << nPFParticleHits
                            << ": nNoParticleHits: " << nNoParticleHits
                            << ": nMultiParticleHits: " << nMultiParticleHits
                            << ": ntrackIdCntSet: " << trackIdCntSet.size()
                            << ": ntrackIdToHit2DMap: " << trackIdToHit2DMap.size();
            ***/


        } // end loop over clusters

        // Fill the map of trackIDs --> PFParticles (through hits)
        if ( !trackIdToHit2DMap.empty() ) {
            // Now spin through the trackIdToHit2DMap to build the reverse map of above, taking track ID's to PFParticles...
            for ( const auto& trackItr : trackIdToHit2DMap ) {
                TrackIDToPFParticles[trackItr.first].push_back( pfParticle );
            }
        } else {
            mf::LogWarning("MCParticleToPFParticle::MakePFParticleMaps") << "===>> No PFParticle to MCParticle match made; the PFParticl has a PDGCode: "
                                                                         << pfParticle->PdgCode();
        }
    } // end of loop over the PFParticle collection


    // Sort the PFParticles in the track ID to PFParticle map
    for ( auto& trackItr : TrackIDToPFParticles ) {
        std::sort( trackItr.second.begin(), trackItr.second.end(), SortPFParticleVec( PFParticleToHitCnt ) );
    }

    return;
}

#endif
