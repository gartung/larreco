#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "lardata/Utilities/ForEachAssociatedGroup.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "cetlib/exception.h"
//
#include <memory>
//
#include "art/Utilities/make_tool.h"
#include "larreco/TrackFinder/TrackMaker.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
//
  /**
   * @file  larreco/TrackFinder/TrackProducerFromPFParticle_module.cc
   * @class TrackProducerFromPFParticle
   *
   * @brief Produce a reco::Track collection, as a result of the fit of an existing recob::PFParticle collection.
   *
   * This producer takes an input an existing recob::PFParticle collection and refits the associated tracks;
   * it can make track fits of the associated showers as well, but this is experimental - do it at your own risk.
   * The mandatory outputs are: the resulting recob::Track collection, the associated hits, and the association
   * between the input PFParticle and the output Track.
   * Optional outputs are recob::TrackFitHitInfo and recob::SpacePoint collections, plus the Assns of SpacePoints to Hits.
   * An option is provided to create SpacePoints from the TrajectoryPoints in the Track.
   * The fit is performed by an user-defined tool, which must inherit from larreco/TrackFinder/TrackMaker.
   *
   * Parameters: trackMaker (fhicl::ParameterSet for the trkmkr::TrackMaker tool used to do the fit), inputCollection (art::InputTag of the input recob::Track collection),
   * doTrackFitHitInfo (bool to decide whether to produce recob::TrackFitHitInfo's), doSpacePoints (bool to decide whether to produce recob::SpacePoint's),
   * spacePointsFromTrajP (bool to decide whether the produced recob::SpacePoint's are taken from the recob::tracking::TrajectoryPoint_t's of the fitted recob::Track),
   * trackFromPF (bool to decide whether to fit the recob::Track associated to the recob::PFParticle),
   * showerFromPF (bool to decide whether to fit the recob::Shower associated to the recob::PFParticle - this option is intended to mitigate possible problems due to tracks being mis-identified as showers),
   * clusterFromPF (bool to decide whether to fit the recob::Clusters associated to the recob::PFParticle).
   *
   * @author  G. Cerati (FNAL, MicroBooNE)
   * @date    2017
   * @version 1.0
   */
//
//
class TrackProducerFromPFParticle : public art::EDProducer {
public:
  explicit TrackProducerFromPFParticle(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  //
  // Plugins should not be copied or assigned.
  TrackProducerFromPFParticle(TrackProducerFromPFParticle const &) = delete;
  TrackProducerFromPFParticle(TrackProducerFromPFParticle &&) = delete;
  TrackProducerFromPFParticle & operator = (TrackProducerFromPFParticle const &) = delete;
  TrackProducerFromPFParticle & operator = (TrackProducerFromPFParticle &&) = delete;
  // Required functions.
  void produce(art::Event & e) override;
private:
  std::unique_ptr<trkmkr::TrackMaker> trackMaker_;
  art::InputTag pfpInputTag;
  art::InputTag trkInputTag;
  art::InputTag shwInputTag;
  art::InputTag clsInputTag;
  bool doTrackFitHitInfo_;
  bool doSpacePoints_;
  bool spacePointsFromTrajP_;
  bool trackFromPF_;
  bool showerFromPF_;
  bool clusterFromPF_;
};
//
TrackProducerFromPFParticle::TrackProducerFromPFParticle(fhicl::ParameterSet const & p)
  : trackMaker_{art::make_tool<trkmkr::TrackMaker>(p.get<fhicl::ParameterSet>("trackMaker"))}
  , pfpInputTag{p.get<art::InputTag>("inputCollection")}
  , doTrackFitHitInfo_{p.get<bool>("doTrackFitHitInfo")}
  , doSpacePoints_{p.get<bool>("doSpacePoints")}
  , spacePointsFromTrajP_{p.get<bool>("spacePointsFromTrajP")}
  , trackFromPF_{p.get<bool>("trackFromPF")}
  , showerFromPF_{p.get<bool>("showerFromPF")}
  , clusterFromPF_{p.get<bool>("clusterFromPF")}
{
  //
  if (p.has_key("trackInputTag")) trkInputTag = p.get<art::InputTag>("trackInputTag");
  else trkInputTag = pfpInputTag;
  if (p.has_key("showerInputTag")) shwInputTag = p.get<art::InputTag>("showerInputTag");
  else shwInputTag = pfpInputTag;
  if (p.has_key("clusterInputTag")) clsInputTag = p.get<art::InputTag>("clusterInputTag");
  else clsInputTag = shwInputTag;
  produces<std::vector<recob::Track> >();
  produces<art::Assns<recob::Track, recob::Hit> >();
  produces<art::Assns<recob::PFParticle, recob::Track> >();
  if (doTrackFitHitInfo_) produces<std::vector<std::vector<recob::TrackFitHitInfo> > >();
  if (doSpacePoints_) {
    produces<std::vector<recob::SpacePoint> >();
    produces<art::Assns<recob::Hit, recob::SpacePoint> >();
  }
}
//
void TrackProducerFromPFParticle::produce(art::Event & e)
{
  //std::cout << "TrackProducerFromPFParticle::produce" << std::endl;
  // Output collections
  auto outputTracks  = std::make_unique<std::vector<recob::Track> >();
  auto outputHits    = std::make_unique<art::Assns<recob::Track, recob::Hit> >();
  auto outputPfpTAssn = std::make_unique<art::Assns<recob::PFParticle, recob::Track> >();
  auto outputHitInfo = std::make_unique<std::vector<std::vector<recob::TrackFitHitInfo> > >();
  auto outputSpacePoints  = std::make_unique<std::vector<recob::SpacePoint> >();
  auto outputHitSpacePointAssn = std::make_unique<art::Assns<recob::Hit, recob::SpacePoint> >();
  //
  // PtrMakers for Assns
  art::PtrMaker<recob::Track> trackPtrMaker(e, *this);
  art::PtrMaker<recob::SpacePoint>* spacePointPtrMaker = nullptr;
  if (doSpacePoints_) spacePointPtrMaker = new art::PtrMaker<recob::SpacePoint>(e, *this);
  //
  // Input from event
  art::ValidHandle<std::vector<recob::PFParticle> > inputPfps = e.getValidHandle<std::vector<recob::PFParticle> >(pfpInputTag);
  const auto assocTracks = std::unique_ptr<art::FindManyP<recob::Track> >(new art::FindManyP<recob::Track>(inputPfps, e, trkInputTag));
  const auto assocShowers = std::unique_ptr<art::FindManyP<recob::Shower> >(new art::FindManyP<recob::Shower>(inputPfps, e, shwInputTag));
  const auto assocClusters = std::unique_ptr<art::FindManyP<recob::Cluster> >(new art::FindManyP<recob::Cluster>(inputPfps, e, clsInputTag));
  //
  const auto* tkHitsAssn = (trackFromPF_ ? e.getValidHandle<art::Assns<recob::Track, recob::Hit> >(trkInputTag).product() : nullptr);
  using TrackHitsGroupsPtr = decltype(util::associated_groups(*tkHitsAssn));
  std::unique_ptr<TrackHitsGroupsPtr> trackHitsGroups = nullptr;
  if (trackFromPF_) {
    trackHitsGroups = std::make_unique<TrackHitsGroupsPtr>(util::associated_groups(*tkHitsAssn));
  }
  //
  auto const& pfClustersAssn = *e.getValidHandle<art::Assns<recob::PFParticle, recob::Cluster> >(clsInputTag);
  const auto& pfpClusterGroups = util::associated_groups(pfClustersAssn);
  auto const& clHitsAssn = *e.getValidHandle<art::Assns<recob::Cluster, recob::Hit> >(clsInputTag);
  const auto& clusterHitsGroups = util::associated_groups(clHitsAssn);
  auto pfpsproxy = proxy::getCollection<std::vector<recob::PFParticle> > (e, pfpInputTag, proxy::withAssociated<recob::SpacePoint>());
  //
  // Initialize tool for this event
  trackMaker_->initEvent(e);
  //
  // Loop over pfps to fit
  for (unsigned int iPfp = 0; iPfp < inputPfps->size(); ++iPfp) {
    //
    const art::Ptr<recob::PFParticle> pfp(inputPfps, iPfp);
    //
    //std::cout << "pfp with key=" << pfp.key() << " iPfp=" << iPfp << " pdg=" << pfp->PdgCode() << " clsize=" << util::groupByIndex(pfpClusterGroups, pfp.key()).size() << std::endl;
    //
    if (trackFromPF_) {
      // Tracks associated to PFParticles
      const std::vector<art::Ptr<recob::Track> >& tracks = assocTracks->at(iPfp);
      // Loop over tracks to refit
      for (art::Ptr<recob::Track> const& track: tracks) {
	//
	// Get track and its hits
	std::vector<art::Ptr<recob::Hit> > inHits;
	decltype(auto) hitsRange = util::groupByIndex(*trackHitsGroups, track.key());
	for (art::Ptr<recob::Hit> const& hit: hitsRange) inHits.push_back(hit);
	//
	// Declare output objects
	recob::Track outTrack;
	std::vector<art::Ptr<recob::Hit> > outHits;
	trkmkr::OptionalOutputs optionals;
	if (doTrackFitHitInfo_) optionals.initTrackFitInfos();
	if (doSpacePoints_ && !spacePointsFromTrajP_) optionals.initSpacePoints();
	//
	// Invoke tool to fit track and fill output objects
	bool fitok = trackMaker_->makeTrack(track, inHits, outTrack, outHits, optionals);
	if (!fitok) continue;
	//
	// Check that the requirement Nhits == Npoints is satisfied
	// We also require the hits to the in the same order as the points (this cannot be enforced, can it?)
	if (outTrack.NumberTrajectoryPoints()!=outHits.size()) {
	  throw cet::exception("TrackProducerFromPFParticle") << "Produced recob::Track required to have 1-1 correspondance between hits and points.\n";
	}
	//
	// Fill output collections, including Assns
	outputTracks->emplace_back(std::move(outTrack));
	const art::Ptr<recob::Track> aptr = trackPtrMaker(outputTracks->size()-1);
	outputPfpTAssn->addSingle(pfp, aptr);
	unsigned int ip = 0;
	for (auto const& trhit: outHits) {
	  outputHits->addSingle(aptr, trhit);
	  //
	  if (doSpacePoints_ && spacePointsFromTrajP_ && outputTracks->back().HasValidPoint(ip)) {
	    auto& tp = outputTracks->back().Trajectory().LocationAtPoint(ip);
	    const double fXYZ[3] = {tp.X(),tp.Y(),tp.Z()};
	    const double fErrXYZ[6] = {0};
	    recob::SpacePoint sp(fXYZ, fErrXYZ, -1.);
	    outputSpacePoints->emplace_back(std::move(sp));
	    const art::Ptr<recob::SpacePoint> apsp = (*spacePointPtrMaker)(outputSpacePoints->size()-1);
	    outputHitSpacePointAssn->addSingle(trhit, apsp);
	  }
	  ip++;
	}
	if (doSpacePoints_ && !spacePointsFromTrajP_) {
	  auto osp = optionals.spacePointHitPairs();
	  for (auto it = osp.begin(); it!=osp.end(); ++it ) {
	    outputSpacePoints->emplace_back(std::move(it->first));
	    const art::Ptr<recob::SpacePoint> apsp = (*spacePointPtrMaker)(outputSpacePoints->size()-1);
	    outputHitSpacePointAssn->addSingle(it->second,apsp);
	  }
	}
	if (doTrackFitHitInfo_) {
	  outputHitInfo->emplace_back(std::move(optionals.trackFitHitInfos()));
	}
      }
    }
    //
    if (showerFromPF_) {
      //
      // Showers associated to PFParticles
      const std::vector<art::Ptr<recob::Shower> >& showers = assocShowers->at(iPfp);
      // if there is more than one shower the logic below to get the hits does not work! this works, at least for uboone
      if (showers.size()!=1) continue;
      //std::cout << "pfp pdg=" << pfp->PdgCode() << std::endl;
      //
      // Loop over showers to refit (should be only one)
      for (unsigned int iShower = 0; iShower < showers.size(); ++iShower) {
	//
	// Get hits for shower (through the chain pfp->clusters->hits)
	std::vector<art::Ptr<recob::Hit> > inHits;
	// Clusters associated to PFParticles
	const std::vector<art::Ptr<recob::Cluster> >& clusters = assocClusters->at(iPfp);
	//std::cout << "pfp associated to clusters with size=" << clusters.size() << std::endl;
	if (clusters.size()<2) continue;
	for (art::Ptr<recob::Cluster> const& cluster: clusters) {
	  std::cout << "cluster found in shower pfp" << std::endl;
	  decltype(auto) hitsRange = util::groupByIndex(clusterHitsGroups, cluster.key());
	  for (art::Ptr<recob::Hit> const& hit: hitsRange) inHits.push_back(hit);
	}
	//
	// Get the shower and convert/hack it into a trajectory so that the fit is initialized
	// (all positions and directions are the same, but only the first one matters in the fitter)
	art::Ptr<recob::Shower> shower = showers[iShower++];
	std::cout << "trying to fit shower with id=" << shower->ID() << std::endl;
	recob::tracking::Point_t pos(shower->ShowerStart().X(),shower->ShowerStart().Y(),shower->ShowerStart().Z());
	recob::tracking::Vector_t dir(shower->Direction().X(),shower->Direction().Y(),shower->Direction().Z());
	std::vector<recob::tracking::Point_t> p;
	std::vector<recob::tracking::Vector_t> d;
	std::cout << "pos=" << pos << " dir=" << dir << " nhits=" << inHits.size() << std::endl;
	for (unsigned int i=0; i<inHits.size(); ++i) {
	  p.push_back(pos);
	  d.push_back(dir);
	}
	recob::TrackTrajectory traj(std::move(p), std::move(d), recob::TrackTrajectory::Flags_t(p.size()), false);
	//
	// Declare output objects
	recob::Track outTrack;
	std::vector<art::Ptr<recob::Hit> > outHits;
	trkmkr::OptionalOutputs optionals;
	if (doTrackFitHitInfo_) optionals.initTrackFitInfos();
	if (doSpacePoints_ && !spacePointsFromTrajP_) optionals.initSpacePoints();
	//
	// Invoke tool to fit track and fill output objects
	std::cout << "before fit" << std::endl;
	bool fitok = trackMaker_->makeTrack(traj, shower->ID(), inHits, outTrack, outHits, optionals);
	std::cout << "fitok=" << fitok << " with nvalid=" << outTrack.CountValidPoints() << std::endl;
	if (!fitok) continue;
	//
	// Check that the requirement Nhits == Npoints is satisfied
	// We also require the hits to the in the same order as the points (this cannot be enforced, can it?)
	if (outTrack.NumberTrajectoryPoints()!=outHits.size()) {
	  throw cet::exception("TrackProducerFromPFParticle") << "Produced recob::Track required to have 1-1 correspondance between hits and points.\n";
	}
	//
	// Fill output collections, including Assns
	outputTracks->emplace_back(std::move(outTrack));
	const art::Ptr<recob::Track> aptr = trackPtrMaker(outputTracks->size()-1);
	outputPfpTAssn->addSingle(pfp, aptr);
	unsigned int ip = 0;
	for (auto const& trhit: outHits) {
	  outputHits->addSingle(aptr, trhit);
	  //
	  if (doSpacePoints_ && spacePointsFromTrajP_ && outputTracks->back().HasValidPoint(ip)) {
	    auto& tp = outputTracks->back().Trajectory().LocationAtPoint(ip);
	    const double fXYZ[3] = {tp.X(),tp.Y(),tp.Z()};
	    const double fErrXYZ[6] = {0};
	    recob::SpacePoint sp(fXYZ, fErrXYZ, -1.);
	    outputSpacePoints->emplace_back(std::move(sp));
	    const art::Ptr<recob::SpacePoint> apsp = (*spacePointPtrMaker)(outputSpacePoints->size()-1);
	    outputHitSpacePointAssn->addSingle(trhit, apsp);
	  }
	  ip++;
	}
	if (doSpacePoints_ && !spacePointsFromTrajP_) {
	  auto osp = optionals.spacePointHitPairs();
	  for (auto it = osp.begin(); it!=osp.end(); ++it ) {
	    outputSpacePoints->emplace_back(std::move(it->first));
	    const art::Ptr<recob::SpacePoint> apsp = (*spacePointPtrMaker)(outputSpacePoints->size()-1);
	    outputHitSpacePointAssn->addSingle(it->second,apsp);
	  }
	}
	if (doTrackFitHitInfo_) {
	  outputHitInfo->emplace_back(std::move(optionals.trackFitHitInfos()));
	}
      }
    }
    //
    if (clusterFromPF_) {
      //
      // Clusters associated to PFParticles
      const std::vector<art::Ptr<recob::Cluster> >& clusters = assocClusters->at(iPfp);
      std::cout << "pfp associated to clusters with size=" << clusters.size() << std::endl;
      if (clusters.size()<2) continue;
      //
      // Get hits for clusters
      std::vector<art::Ptr<recob::Hit> > inHits;
      for (art::Ptr<recob::Cluster> const& cluster: clusters) {
	decltype(auto) hitsRange = util::groupByIndex(clusterHitsGroups, cluster.key());
	for (art::Ptr<recob::Hit> const& hit: hitsRange) inHits.push_back(hit);
      }
      //
      // Get the space points and convert/hack the inputs into a trajectory so that the fit is initialized
      auto ppxy = pfpsproxy[pfp.key()];
      auto& sps = ppxy.get<recob::SpacePoint>();
      std::cout << "nHits=" << inHits.size() << " nSpacePoints=" << sps.size() << std::endl;
      if (sps.size()<2) continue;
      auto& sp0 = sps[0];
      auto& sp1 = sps[std::min(size_t(10),sps.size()-1)];
      const double* xyz0 = sp0->XYZ();
      recob::tracking::Point_t pos(xyz0[0],xyz0[1],xyz0[2]);
      const double* xyz1 = sp1->XYZ();
      recob::tracking::Vector_t dir(recob::tracking::Point_t(xyz1[0],xyz1[1],xyz1[2]) - pos);
      dir = dir.Unit();
      std::cout << "pos0=" << pos << std::endl;
      std::cout << "pos1=" << recob::tracking::Point_t(xyz1[0],xyz1[1],xyz1[2]) << std::endl;
      std::cout << "dir=" << dir << std::endl;
      std::vector<recob::tracking::Point_t> p;
      std::vector<recob::tracking::Vector_t> d;
      for (unsigned int i=0; i<inHits.size(); ++i) {
	p.push_back(pos);
	d.push_back(dir);
      }
      recob::TrackTrajectory traj(std::move(p), std::move(d), recob::TrackTrajectory::Flags_t(p.size()), false);
      //
      // Declare output objects
      recob::Track outTrack;
      std::vector<art::Ptr<recob::Hit> > outHits;
      trkmkr::OptionalOutputs optionals;
      if (doTrackFitHitInfo_) optionals.initTrackFitInfos();
      if (doSpacePoints_ && !spacePointsFromTrajP_) optionals.initSpacePoints();
      //
      // Invoke tool to fit track and fill output objects
      bool fitok = trackMaker_->makeTrack(traj, outputTracks->size(), inHits, outTrack, outHits, optionals);
      if (!fitok) continue;
      //
      // Check that the requirement Nhits == Npoints is satisfied
      // We also require the hits to the in the same order as the points (this cannot be enforced, can it?)
      if (outTrack.NumberTrajectoryPoints()!=outHits.size()) {
	throw cet::exception("TrackProducerFromPFParticle") << "Produced recob::Track required to have 1-1 correspondance between hits and points.\n";
      }
      //
      // Fill output collections, including Assns
      outputTracks->emplace_back(std::move(outTrack));
      const art::Ptr<recob::Track> aptr = trackPtrMaker(outputTracks->size()-1);
      outputPfpTAssn->addSingle(pfp, aptr);
      unsigned int ip = 0;
      for (auto const& trhit: outHits) {
	outputHits->addSingle(aptr, trhit);
	//
	if (doSpacePoints_ && spacePointsFromTrajP_ && outputTracks->back().HasValidPoint(ip)) {
	  auto& tp = outputTracks->back().Trajectory().LocationAtPoint(ip);
	  const double fXYZ[3] = {tp.X(),tp.Y(),tp.Z()};
	  const double fErrXYZ[6] = {0};
	  recob::SpacePoint sp(fXYZ, fErrXYZ, -1.);
	  outputSpacePoints->emplace_back(std::move(sp));
	  const art::Ptr<recob::SpacePoint> apsp = (*spacePointPtrMaker)(outputSpacePoints->size()-1);
	  outputHitSpacePointAssn->addSingle(trhit, apsp);
	}
	ip++;
      }
      if (doSpacePoints_ && !spacePointsFromTrajP_) {
	auto osp = optionals.spacePointHitPairs();
	for (auto it = osp.begin(); it!=osp.end(); ++it ) {
	  outputSpacePoints->emplace_back(std::move(it->first));
	  const art::Ptr<recob::SpacePoint> apsp = (*spacePointPtrMaker)(outputSpacePoints->size()-1);
	  outputHitSpacePointAssn->addSingle(it->second,apsp);
	}
      }
      if (doTrackFitHitInfo_) {
	outputHitInfo->emplace_back(std::move(optionals.trackFitHitInfos()));
      }
    }
    //
  }
  //
  // Put collections in the event
  e.put(std::move(outputTracks));
  e.put(std::move(outputHits));
  e.put(std::move(outputPfpTAssn));
  if (doTrackFitHitInfo_) {
    e.put(std::move(outputHitInfo));
  }
  if (doSpacePoints_) {
    e.put(std::move(outputSpacePoints));
    e.put(std::move(outputHitSpacePointAssn));
  }
  if (doSpacePoints_) delete spacePointPtrMaker;
}
//
DEFINE_ART_MODULE(TrackProducerFromPFParticle)
