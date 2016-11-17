////////////////////////////////////////////////////////////////////////////
// Class:       TrackIdentifier
// Module Type: producer
// File:        TrackIdentifier_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffield.ac.uk), November 2016
//
// Module to separate tracks from showers.
//
// Produces recob::Tracks (previously found using other techniques),
// retaining all associations, which are determined to be actual tracks.
// Also produces a set of hits, determined to be shower hits.  The
// intention is that shower finding algorithms (e.g. BlurredCluster &
// EMShower run over this output to create showers.
//
// This module is basically a wrapper around the TrackShowerSeparationAlg
// algorithm, used by BlurredCluster to idetify shower-like hits.  It
// can be used as a way of considering only track-like tracks in analyses.
// NB It may make more sense for this to be included in the tracking
//    algorithms themselves (e.g. PMTrack).  This should be considered.
////////////////////////////////////////////////////////////////////////////

// Framework includes:
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/RecoAlg/TrackShowerSepAlg.h"

// ROOT includes

namespace shower {
  class TrackIdentifier;
}

class shower::TrackIdentifier : public art::EDProducer {
public:

  TrackIdentifier(const fhicl::ParameterSet& pset);

  void produce(art::Event& evt);
  void reconfigure(const fhicl::ParameterSet& pset);

private:

  std::string fHitsModuleLabel, fTrackModuleLabel;

  TrackShowerSepAlg fTrackShowerSepAlg;

};

DEFINE_ART_MODULE(shower::TrackIdentifier)

shower::TrackIdentifier::TrackIdentifier(const fhicl::ParameterSet& pset) : fTrackShowerSepAlg(pset.get<fhicl::ParameterSet>("TrackShowerSepAlg")) {
  this->reconfigure(pset);
  produces<std::vector<recob::Track> >();
  produces<std::vector<recob::Hit> >();
  produces<art::Assns<recob::Track, recob::Hit> >();
  produces<art::Assns<recob::Track, recob::Cluster> >();
  produces<art::Assns<recob::Track, recob::SpacePoint> >();
  produces<art::Assns<recob::Track, recob::Vertex> >();
}

void shower::TrackIdentifier::reconfigure(const fhicl::ParameterSet& pset) {
  fHitsModuleLabel  = pset.get<std::string>("HitsModuleLabel");
  fTrackModuleLabel = pset.get<std::string>("TrackModuleLabel");
}

void shower::TrackIdentifier::produce(art::Event& evt) {

  // Output data products
  std::unique_ptr<std::vector<recob::Track> > outTracks(new std::vector<recob::Track>);
  std::unique_ptr<std::vector<recob::Hit> > outHits(new std::vector<recob::Hit>);
  std::unique_ptr<art::Assns<recob::Track, recob::Hit> > hitAssns(new art::Assns<recob::Track, recob::Hit>);
  std::unique_ptr<art::Assns<recob::Track, recob::Cluster> > clusterAssns(new art::Assns<recob::Track, recob::Cluster>);
  std::unique_ptr<art::Assns<recob::Track, recob::SpacePoint> > spacePointAssns(new art::Assns<recob::Track, recob::SpacePoint>);
  std::unique_ptr<art::Assns<recob::Track, recob::Vertex> > vertexAssns(new art::Assns<recob::Track, recob::Vertex>);

  // Input data products
  std::vector<art::Ptr<recob::Hit> > hits;
  art::Handle<std::vector<recob::Hit> > hitHandle;
  if (evt.getByLabel(fHitsModuleLabel, hitHandle))
    art::fill_ptr_vector(hits, hitHandle);

  std::vector<art::Ptr<recob::Track> > tracks;
  art::Handle<std::vector<recob::Track> > trackHandle;
  if (evt.getByLabel(fTrackModuleLabel, trackHandle))
    art::fill_ptr_vector(tracks, trackHandle);

  std::vector<art::Ptr<recob::SpacePoint> > spacePoints;
  art::Handle<std::vector<recob::SpacePoint> > spacePointHandle;
  if (evt.getByLabel(fTrackModuleLabel, spacePointHandle))
    art::fill_ptr_vector(spacePoints, spacePointHandle);

  art::FindManyP<recob::Hit>        fmht(trackHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Track>      fmth(hitHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::SpacePoint> fmspt(trackHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Track>      fmtsp(spacePointHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Cluster>    fmct(trackHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Vertex>     fmvt(trackHandle, evt, fTrackModuleLabel);

  // Run track shower separation
  fTrackShowerSepAlg.RunTrackShowerSep(evt.event(), hits, tracks, spacePoints, fmht, fmth, fmspt, fmtsp);

  // Get the output collections
  const std::vector<art::Ptr<recob::Track> > trackTracks = fTrackShowerSepAlg.TrackTracks();
  const std::vector<art::Ptr<recob::Hit> > showerHits = fTrackShowerSepAlg.ShowerHits();

  // Make output data products
  for (std::vector<art::Ptr<recob::Track> >::const_iterator trackIt = trackTracks.begin(); trackIt != trackTracks.end(); ++trackIt) {
    outTracks->push_back(*(*trackIt));
    std::vector<art::Ptr<recob::Hit> > trackHits = fmht.at(trackIt->key());
    std::vector<art::Ptr<recob::Cluster> > trackClusters = fmct.at(trackIt->key());
    std::vector<art::Ptr<recob::SpacePoint> > trackSpacePoints = fmspt.at(trackIt->key());
    std::vector<art::Ptr<recob::Vertex> > trackVertices = fmvt.at(trackIt->key());
    util::CreateAssn(*this, evt, *(outTracks.get()), trackHits,        *(hitAssns.get()));
    util::CreateAssn(*this, evt, *(outTracks.get()), trackClusters,    *(clusterAssns.get()));
    util::CreateAssn(*this, evt, *(outTracks.get()), trackSpacePoints, *(spacePointAssns.get()));
    util::CreateAssn(*this, evt, *(outTracks.get()), trackVertices,    *(vertexAssns.get()));
  }
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = showerHits.begin(); hitIt != showerHits.end(); ++hitIt)
    outHits->push_back(*(*hitIt));

  evt.put(std::move(outTracks));
  evt.put(std::move(outHits));
  evt.put(std::move(hitAssns));
  evt.put(std::move(clusterAssns));
  evt.put(std::move(spacePointAssns));
  evt.put(std::move(vertexAssns));

  return;

}
