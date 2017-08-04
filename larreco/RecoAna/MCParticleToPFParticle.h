/**
 * @file   MCParticleToPFParticle.h
 * @brief  Algorithms to match reconstructed hits to track IDs in MC particles (through MC hits),
 *         and to match track IDs to PFParticles
 * @author Yun-Tse Tsai  (stealing codes from Tracy Usher)
 * @date   August 3rd, 2017
 *
 */

#ifndef LARRECO_ANA_MCPARTICLETOPFPARTICLE_H
#define LARRECO_ANA_MCPARTICLETOPFPARTICLE_H

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/MCBase/MCHitCollection.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/cpu_timer.h"
#include "cetlib/exception.h"

// ROOT
#include "TVector3.h"

// C/C++ standard libraries
#include <vector>
#include <map>
#include <set>
#include <string>
#include <algorithm>

namespace lar {
    namespace recoana {

    class MCParticleToPFParticle {

        public:

        MCParticleToPFParticle( std::string pfparticleLabel, std::string clusterLabel )
            : fPFParticleProducerLabel( pfparticleLabel ), fClusterProducerLabel( clusterLabel )
            {}

        void setup( detinfo::DetectorClocks const& detclock, geo::GeometryCore const& geometry )
            { detc = &detclock; geom = &geometry; }


        // Useful typedefs to keep from getting lost in stl'ese
        typedef std::map< const recob::Hit*, std::set< int > > HitToParticleMap;    // Maps recob::Hits to track ID's (hence MCParticles)
        typedef std::map< int, std::set< const recob::Hit* > > ParticleToHitMap;    // Maps track IDs to recob Hits
        typedef std::vector< art::Ptr< recob::Hit > >          HitVector;
    
        // More useful typedefs, this time for relating PFParticles to Track IDs and hits
        typedef std::set< art::Ptr< recob::Hit > >                           Hit2DSet;                   // Need to count only unique hits
        typedef std::map< int, Hit2DSet >                                    TrackIDToHit2DMap;          // Maps Track ID's to lists of hits
        typedef std::vector< art::Ptr< recob::PFParticle > >                 PFParticleVec;              // Typedef a vector of PFParticles
        typedef std::map< int, PFParticleVec >                               TrackIDToPFParticleVecMap;  // Maps Track ID's to the PFParticles associated to them
        typedef std::map< art::Ptr< recob::PFParticle >, TrackIDToHit2DMap > PFParticleToTrackHit2DMap;  // Maps PFParticles to a map of associated Track ID's to their hits
        typedef std::map< art::Ptr< recob::PFParticle >, int >               PFParticleHitCntMap;        // Allows us to count reco hits per PFParticle directly


        // Define a function to fill the above data structures
        void MakeHitParticleMaps( const std::vector< sim::MCHitCollection >& MCHits,
                                  const std::vector< recob::Hit >&           Hits,
                                  HitToParticleMap&                          HitToParticle,
                                  ParticleToHitMap&                          ParticleToHit );
    
        // Define a function to fill all of the above
        void MakePFParticleMaps( const art::Event&                                      event,
                                 const art::Handle< std::vector< recob::PFParticle > >& pfParticleHandle,
                                 const HitToParticleMap&                                HitToParticle,
                                 PFParticleToTrackHit2DMap&                             PFParticleToTrackHit2D,
                                 TrackIDToPFParticleVecMap&                             TrackIDtoPFParticles,
                                 PFParticleHitCntMap&                                   PFParticleHitCnt );

        private:

        std::string fPFParticleProducerLabel;
        std::string fClusterProducerLabel;
        class SortPFParticleVec;

        /// Pointer to the geometry to be used
        geo::GeometryCore const* geom = nullptr;

        /// Pointer to the detector clock
        detinfo::DetectorClocks const* detc = nullptr;

    }; // class MCParticleToPFParticle

    // This will be used to sort PFParticles in the map below
    class MCParticleToPFParticle::SortPFParticleVec {
        /**
         * @brief This is used to sort "Hough Clusters" by the maximum entries in a bin
         */
        public:

        SortPFParticleVec( const MCParticleToPFParticle::PFParticleHitCntMap& pfPartCntMap ) : fPFPartCntMap( pfPartCntMap ) {}
    
        bool operator()( const art::Ptr< recob::PFParticle >& left, const art::Ptr< recob::PFParticle >& right ) {
            size_t numHitsLeft  = fPFPartCntMap.at( left );
            size_t numHitsRight = fPFPartCntMap.at( right );
    
            return numHitsLeft > numHitsRight;
        }

        private:
        const MCParticleToPFParticle::PFParticleHitCntMap& fPFPartCntMap;
    }; // class MCParticleToPFParticle::SortPFParticleVec


    } // namespace recoana
} // namespace lar

#endif
