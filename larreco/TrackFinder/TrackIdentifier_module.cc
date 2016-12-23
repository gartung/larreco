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
// EMShower) run over this output to create showers.
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
  std::string fTrackInstanceLabel, fShowerInstanceLabel;

  TrackShowerSepAlg fTrackShowerSepAlg;

};

DEFINE_ART_MODULE(shower::TrackIdentifier)

shower::TrackIdentifier::TrackIdentifier(const fhicl::ParameterSet& pset) : fTrackShowerSepAlg(pset.get<fhicl::ParameterSet>("TrackShowerSepAlg")) {
  this->reconfigure(pset);
  produces<std::vector<recob::Track> >(fTrackInstanceLabel);
  produces<std::vector<recob::Hit> >(fShowerInstanceLabel);
  produces<std::vector<recob::Vertex> >(fShowerInstanceLabel);
  produces<art::Assns<recob::Hit, recob::Track> >(fShowerInstanceLabel);
  produces<art::Assns<recob::Track, recob::Hit> >(fTrackInstanceLabel);
  produces<art::Assns<recob::Track, recob::Vertex> >(fTrackInstanceLabel);
  produces<art::Assns<recob::Track, recob::SpacePoint> >(fTrackInstanceLabel);
  produces<art::Assns<recob::SpacePoint, recob::Hit> >(fTrackInstanceLabel);
}

void shower::TrackIdentifier::reconfigure(const fhicl::ParameterSet& pset) {
  fHitsModuleLabel     = pset.get<std::string>("HitsModuleLabel");
  fTrackModuleLabel    = pset.get<std::string>("TrackModuleLabel");
  fTrackInstanceLabel  = pset.get<std::string>("TrackInstanceLabel","track");
  fShowerInstanceLabel = pset.get<std::string>("ShowerInstanceLabel","shower");
}

void shower::TrackIdentifier::produce(art::Event& evt) {

  // Output data products
  // Currently making all associations which PMTrack does itself, except for PFParticles and track meta data
  std::unique_ptr<std::vector<recob::Track> > outTracks(new std::vector<recob::Track>);
  std::unique_ptr<std::vector<recob::Hit> > outHits(new std::vector<recob::Hit>);
  std::unique_ptr<std::vector<recob::Vertex> > outVertices(new std::vector<recob::Vertex>);
  std::unique_ptr<art::Assns<recob::Hit, recob::Track> > hitTrackAssns(new art::Assns<recob::Hit, recob::Track>);
  std::unique_ptr<art::Assns<recob::Track, recob::Hit> > trackHitAssns(new art::Assns<recob::Track, recob::Hit>);
  std::unique_ptr<art::Assns<recob::Track, recob::Vertex> > vertexAssns(new art::Assns<recob::Track, recob::Vertex>);
  std::unique_ptr<art::Assns<recob::Track, recob::SpacePoint> > spacePointAssns(new art::Assns<recob::Track, recob::SpacePoint>);
  std::unique_ptr<art::Assns<recob::SpacePoint, recob::Hit> > spHitAssns(new art::Assns<recob::SpacePoint, recob::Hit>);

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
  art::FindManyP<recob::Vertex>     fmvt(trackHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Hit>        fmhsp(spacePointHandle, evt, fTrackModuleLabel);

  // Run track shower separation
  fTrackShowerSepAlg.RunTrackShowerSep(evt.event(), hits, tracks, spacePoints, fmht, fmth, fmspt, fmtsp);

  // Get the output collections
  const std::vector<art::Ptr<recob::Track> > trackTracks = fTrackShowerSepAlg.TrackTracks();
  const std::vector<art::Ptr<recob::Hit> > showerHits = fTrackShowerSepAlg.ShowerHits();
  const std::vector<TVector3> showerStarts = fTrackShowerSepAlg.ShowerStarts();

  // Make output data products

  // Track
  art::ProductID trackID = getProductID<std::vector<recob::Track> >(evt, fTrackInstanceLabel);
  const art::EDProductGetter* trackGetter = evt.productGetter(trackID);
  for (std::vector<art::Ptr<recob::Track> >::const_iterator trackIt = trackTracks.begin(); trackIt != trackTracks.end(); ++trackIt) {
    outTracks->push_back(*(*trackIt));
    art::Ptr<recob::Track> trackPtr(trackID, std::distance(trackTracks.begin(), trackIt), trackGetter);
    const std::vector<art::Ptr<recob::Hit> > trackHits = fmht.at(trackIt->key());
    const std::vector<art::Ptr<recob::SpacePoint> > trackSpacePoints = fmspt.at(trackIt->key());
    const std::vector<art::Ptr<recob::Vertex> > trackVertices = fmvt.at(trackIt->key());
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator trackHitIt = trackHits.begin(); trackHitIt != trackHits.end(); ++trackHitIt)
      trackHitAssns->addSingle(trackPtr, *trackHitIt);
    for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator trackSpIt = trackSpacePoints.begin(); trackSpIt != trackSpacePoints.end(); ++trackSpIt)
      spacePointAssns->addSingle(trackPtr, *trackSpIt);
    for (std::vector<art::Ptr<recob::Vertex> >::const_iterator trackVertexIt = trackVertices.begin(); trackVertexIt != trackVertices.end(); ++trackVertexIt)
      vertexAssns->addSingle(trackPtr, *trackVertexIt);
  }

  // Track space points
  for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spIt = spacePoints.begin(); spIt != spacePoints.end(); ++spIt) {
    const std::vector<art::Ptr<recob::Hit> > spHits = fmhsp.at(spIt->key());
    util::CreateAssn(*this, evt, *spIt, spHits, *(spHitAssns.get()));
  }

  // Shower
  art::ProductID showerHitID = getProductID<std::vector<recob::Hit> >(evt, fShowerInstanceLabel);
  const art::EDProductGetter* showerHitGetter = evt.productGetter(showerHitID);
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = showerHits.begin(); hitIt != showerHits.end(); ++hitIt) {
    outHits->push_back(*(*hitIt));
    art::Ptr<recob::Hit> hitPtr(showerHitID, std::distance(showerHits.begin(), hitIt), showerHitGetter);
    const std::vector<art::Ptr<recob::Track> > hitTracks = fmth.at(hitIt->key());
    for (std::vector<art::Ptr<recob::Track> >::const_iterator trackIt = hitTracks.begin(); trackIt != hitTracks.end(); ++trackIt)
      hitTrackAssns->addSingle(hitPtr, *trackIt);
  }

  // Vertices
  for (std::vector<TVector3>::const_iterator vertexIt = showerStarts.begin(); vertexIt != showerStarts.end(); ++vertexIt) {
    double point[3] = {vertexIt->X(), vertexIt->Y(), vertexIt->Z()};
    outVertices->emplace_back(point, outVertices->size());
  }

  evt.put(std::move(outTracks), fTrackInstanceLabel);
  evt.put(std::move(outHits), fShowerInstanceLabel);
  evt.put(std::move(outVertices), fShowerInstanceLabel);
  evt.put(std::move(hitTrackAssns), fShowerInstanceLabel);
  evt.put(std::move(trackHitAssns), fTrackInstanceLabel);
  evt.put(std::move(spacePointAssns), fTrackInstanceLabel);
  evt.put(std::move(vertexAssns), fTrackInstanceLabel);
  evt.put(std::move(spHitAssns), fTrackInstanceLabel);

  return;

}
