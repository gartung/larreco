////////////////////////////////////////////////////////////////////////
// Class:       BlurredCluster
// Module Type: producer
// File:        BlurredCluster_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffield.ac.uk), May 2015
//
// Reconstructs showers by blurring the hit map image to introduce fake
// hits before clustering to make fuller and more complete clusters.
//
// See DUNE-DocDB 54 (public) for details.
////////////////////////////////////////////////////////////////////////

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

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/ClusterFinder/ClusterCreator.h"
#include "larreco/RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "larreco/RecoAlg/ClusterParamsImportWrapper.h"
#include "larreco/RecoAlg/BlurredClusterAlg.h"
#include "larreco/RecoAlg/MergeClusterAlg.h"
#include "larreco/RecoAlg/TrackShowerSepAlg.h"

// ROOT & C++ includes
#include <string>
#include <vector>
#include <map>

namespace cluster {
  class BlurredCluster;
}

class cluster::BlurredCluster: public art::EDProducer {
public:

  explicit BlurredCluster(fhicl::ParameterSet const& pset);
  virtual ~BlurredCluster();

  void produce(art::Event &evt);
  void reconfigure(fhicl::ParameterSet const &p);

private:

  int fEvent, fRun, fSubrun;

  // Input data products
  std::string fHitsModuleLabel, fTrackModuleLabel, fVertexModuleLabel, fPFParticleModuleLabel;

  // Run options
  bool fCreateDebugPDF, fMergeClusters, fGlobalTPCRecon, fShowerReconOnly, fUseVertices, fUseReblurring;

  // Create instances of algorithm classes to perform the clustering
  cluster::BlurredClusterAlg fBlurredClusterAlg;
  cluster::MergeClusterAlg fMergeClusterAlg;
  shower::TrackShowerSepAlg fTrackShowerSepAlg;

  // Output containers to place in event
  std::unique_ptr<std::vector<recob::Cluster> > clusters;
  std::unique_ptr<art::Assns<recob::Cluster,recob::Hit> > associations;

};

cluster::BlurredCluster::BlurredCluster(fhicl::ParameterSet const &pset) : fBlurredClusterAlg(pset.get<fhicl::ParameterSet>("BlurredClusterAlg")),
									   fMergeClusterAlg(pset.get<fhicl::ParameterSet>("MergeClusterAlg")),
									   fTrackShowerSepAlg(pset.get<fhicl::ParameterSet>("TrackShowerSepAlg")) {
  this->reconfigure(pset);
  produces<std::vector<recob::Cluster> >();
  produces<art::Assns<recob::Cluster,recob::Hit> >();
}

cluster::BlurredCluster::~BlurredCluster() { }

void cluster::BlurredCluster::reconfigure(fhicl::ParameterSet const& p) {
  fHitsModuleLabel       = p.get<std::string>("HitsModuleLabel");
  fTrackModuleLabel      = p.get<std::string>("TrackModuleLabel");
  fVertexModuleLabel     = p.get<std::string>("VertexModuleLabel");
  fPFParticleModuleLabel = p.get<std::string>("PFParticleModuleLabel");
  fCreateDebugPDF        = p.get<bool>       ("CreateDebugPDF");
  fMergeClusters         = p.get<bool>       ("MergeClusters");
  fGlobalTPCRecon        = p.get<bool>       ("GlobalTPCRecon");
  fShowerReconOnly       = p.get<bool>       ("ShowerReconOnly");
  fUseVertices           = p.get<bool>       ("UseVertices");
  fUseReblurring         = p.get<bool>       ("UseReblurring");
  fBlurredClusterAlg.reconfigure(p.get<fhicl::ParameterSet>("BlurredClusterAlg"));
  fMergeClusterAlg.reconfigure(p.get<fhicl::ParameterSet>("MergeClusterAlg"));
  fTrackShowerSepAlg.reconfigure(p.get<fhicl::ParameterSet>("TrackShowerSepAlg"));
}

void cluster::BlurredCluster::produce(art::Event &evt) {

  fEvent  = evt.event();
  fRun    = evt.run();
  fSubrun = evt.subRun();

  // Create debug pdf to illustrate the blurring process
  if (fCreateDebugPDF)
    fBlurredClusterAlg.CreateDebugPDF(fRun, fSubrun, fEvent);

  // Output containers -- collection of clusters and associations
  clusters.reset(new std::vector<recob::Cluster>);
  associations.reset(new art::Assns<recob::Cluster,recob::Hit>);

  // Compute the cluster characteristics
  // Just use default for now, but configuration will go here
  ClusterParamsImportWrapper<StandardClusterParamsAlg> ClusterParamAlgo;

  // Create geometry handle
  art::ServiceHandle<geo::Geometry> geom;

  // Get the hits from the event
  art::Handle<std::vector<recob::Hit> > hitCollection;
  std::vector<art::Ptr<recob::Hit> > hits, hitsToCluster;
  if (evt.getByLabel(fHitsModuleLabel,hitCollection))
    art::fill_ptr_vector(hits, hitCollection);

  // Get vertices from the event
  art::Handle<std::vector<recob::Vertex> > vertexCollection;
  std::vector<art::Ptr<recob::Vertex> > vertices;
  if (evt.getByLabel(fVertexModuleLabel, vertexCollection))
    art::fill_ptr_vector(vertices, vertexCollection);

  if (fUseVertices and !vertexCollection.isValid()) {
    std::cout << "BlurredCluster: UseVertices is selected but no vertices (module label " << fVertexModuleLabel << ") found in event" << std::endl;
    throw;
  }

  // Remove track hits if required
  if (fShowerReconOnly) {

    // // Get the tracks from the event
    // art::Handle<std::vector<recob::Track> > trackCollection;
    // std::vector<art::Ptr<recob::Track> > tracks;
    // if (evt.getByLabel(fTrackModuleLabel,trackCollection))
    //   art::fill_ptr_vector(tracks, trackCollection);

    // // Get the space points from the event
    // art::Handle<std::vector<recob::SpacePoint> > spacePointCollection;
    // std::vector<art::Ptr<recob::SpacePoint> > spacePoints;
    // if (evt.getByLabel(fTrackModuleLabel,spacePointCollection))
    //   art::fill_ptr_vector(spacePoints, spacePointCollection);

    // if (trackCollection.isValid()) {
    //   art::FindManyP<recob::Hit> fmht(trackCollection, evt, fTrackModuleLabel);
    //   art::FindManyP<recob::Track> fmth(hitCollection, evt, fTrackModuleLabel);
    //   art::FindManyP<recob::SpacePoint> fmspt(trackCollection, evt, fTrackModuleLabel);
    //   art::FindManyP<recob::Track> fmtsp(spacePointCollection, evt, fTrackModuleLabel);
    //   fTrackShowerSepAlg.RunTrackShowerSep(evt.event(), hits, tracks, spacePoints, fmht, fmth, fmspt, fmtsp);
    //   hitsToCluster = fTrackShowerSepAlg.ShowerHits();
    // }
    // else
    //   hitsToCluster = hits;

    // Get Pandora PFParticles from the event
    art::Handle<std::vector<recob::PFParticle> > pfParticleHandle;
    std::vector<art::Ptr<recob::PFParticle> > pfParticles;
    if (evt.getByLabel(fPFParticleModuleLabel, pfParticleHandle))
      art::fill_ptr_vector(pfParticles, pfParticleHandle);

    art::Handle<std::vector<recob::Cluster> > clusterHandle;;
    std::vector<art::Ptr<recob::Cluster> > clusters;
    if (evt.getByLabel(fPFParticleModuleLabel, clusterHandle))
      art::fill_ptr_vector(clusters, clusterHandle);

    art::FindManyP<recob::Cluster> fmcpfp(pfParticleHandle, evt, fPFParticleModuleLabel);
    art::FindManyP<recob::Hit> fmhc(clusterHandle, evt, fPFParticleModuleLabel);

    for (std::vector<art::Ptr<recob::PFParticle> >::const_iterator pfParticleIt = pfParticles.begin(); pfParticleIt != pfParticles.end(); ++pfParticleIt) {
      if ((*pfParticleIt)->PdgCode() != 11)
	continue;
      const std::vector<art::Ptr<recob::Cluster> > pfParticleClusters = fmcpfp.at(pfParticleIt->key());
      for (std::vector<art::Ptr<recob::Cluster> >::const_iterator pfpClusterIt = pfParticleClusters.begin(); pfpClusterIt != pfParticleClusters.end(); ++pfpClusterIt) {
	const std::vector<art::Ptr<recob::Hit> > pfpClusterHits = fmhc.at(pfpClusterIt->key());
        for (std::vector<art::Ptr<recob::Hit> >::const_iterator pfpClusterHitIt = pfpClusterHits.begin(); pfpClusterHitIt != pfpClusterHits.end(); ++pfpClusterHitIt)
	  hitsToCluster.push_back(*pfpClusterHitIt);
      }
    }

  }

  else
    hitsToCluster = hits;

  // Make a map between the planes and the hits on each
  std::map<geo::PlaneID,std::vector<art::Ptr<recob::Hit> > > planeToHits;
  for (std::vector<art::Ptr<recob::Hit> >::iterator hitToCluster = hitsToCluster.begin(); hitToCluster != hitsToCluster.end(); ++hitToCluster) {
    if (fGlobalTPCRecon)
      planeToHits[geo::PlaneID((*hitToCluster)->WireID().Cryostat,(*hitToCluster)->WireID().TPC%2,(*hitToCluster)->WireID().Plane)].push_back(*hitToCluster);
    else
      planeToHits[(*hitToCluster)->WireID().planeID()].push_back(*hitToCluster);
  }

  // Loop over views
  for (std::map<geo::PlaneID,std::vector<art::Ptr<recob::Hit> > >::const_iterator planeIt = planeToHits.begin(); planeIt != planeToHits.end(); ++planeIt) {

    //std::cout << "Clustering in plane " << planeIt->first.Plane << " in global TPC " << planeIt->first.TPC << std::endl;

    std::vector<art::PtrVector<recob::Hit> > finalClusters;

    // Implement the algorithm
    if (planeIt->second.size() >= fBlurredClusterAlg.GetMinSize()) {

      std::vector<std::vector<double> > image = fBlurredClusterAlg.ConvertRecobHitsToVector(planeIt->second);

      // Add vertices
      std::vector<TVector2> planeVertices;
      if (fUseVertices) {
	for (std::vector<art::Ptr<recob::Vertex> >::const_iterator vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt) {
	  double xyz[3]; (*vertexIt)->XYZ(xyz);
	  planeVertices.push_back(fBlurredClusterAlg.Convert3DPointToPlaneBins(TVector3(xyz), planeIt->first.Plane));
	}
      }

      // Find clusters in hit map
      std::vector<std::vector<int> > allClusterBins;
      if (fUseReblurring)
	allClusterBins = fBlurredClusterAlg.FindClusters(image, planeVertices, true);
      else {
	std::vector<std::vector<double> > blurred = fBlurredClusterAlg.GaussianBlur(image);
	allClusterBins = fBlurredClusterAlg.FindClusters(blurred, planeVertices);
      }
      mf::LogVerbatim("Blurred Clustering") << "Found " << allClusterBins.size() << " clusters" << std::endl;

      // Create output clusters from the vector of clusters made in FindClusters
      std::vector<art::PtrVector<recob::Hit> > planeClusters;
      fBlurredClusterAlg.ConvertBinsToClusters(image, allClusterBins, planeClusters);

      // Use the cluster merging algorithm
      if (fMergeClusters) {
	int numMergedClusters = fMergeClusterAlg.MergeClusters(planeClusters, finalClusters);
	mf::LogVerbatim("Blurred Clustering") << "After merging, there are " << numMergedClusters << " clusters" << std::endl;
      }
      else finalClusters = planeClusters;

      // Make the debug PDF
      if (fCreateDebugPDF) {
	std::stringstream name;
	name << "blurred_image";
	TH2F* imageHist = fBlurredClusterAlg.MakeHistogram(image, TString(name.str()));
	name << "_convolved";
	std::vector<std::vector<double> > blurred = fBlurredClusterAlg.GaussianBlur(image);
	TH2F* blurredHist = fBlurredClusterAlg.MakeHistogram(blurred, TString(name.str()));
      	fBlurredClusterAlg.SaveImage(imageHist, 1, planeIt->first.TPC, planeIt->first.Plane);
      	fBlurredClusterAlg.SaveImage(blurredHist, 2, planeIt->first.TPC, planeIt->first.Plane);
      	fBlurredClusterAlg.SaveImage(blurredHist, allClusterBins, 3, planeIt->first.TPC, planeIt->first.Plane);
      	fBlurredClusterAlg.SaveImage(imageHist, finalClusters, 4, planeIt->first.TPC, planeIt->first.Plane);
	imageHist->Delete();
	blurredHist->Delete();
      }

    } // End min hits check

    //fBlurredClusterAlg.fHitMap.clear();

    // Make the output cluster objects
    for (std::vector<art::PtrVector<recob::Hit> >::iterator clusIt = finalClusters.begin(); clusIt != finalClusters.end(); ++clusIt) {

      art::PtrVector<recob::Hit> clusterHits = *clusIt;
      if (clusterHits.size() > 0) {

	// Get the start and end wires of the cluster
	// unsigned int startWire = fBlurredClusterAlg.GlobalWire(clusterHits.front()->WireID());
	// unsigned int endWire = fBlurredClusterAlg.GlobalWire(clusterHits.back()->WireID());
	unsigned int startWire = clusterHits.front()->WireID().Wire;
	unsigned int endWire = clusterHits.back()->WireID().Wire;

	// Put cluster hits in the algorithm
	ClusterParamAlgo.ImportHits(clusterHits);

	// Create the recob::Cluster and place in the vector of clusters
	ClusterCreator cluster(
			       ClusterParamAlgo,                        // algo
			       float(startWire),                        // start_wire
			       0.,                                      // sigma_start_wire
			       clusterHits.front()->PeakTime(),         // start_tick
			       clusterHits.front()->SigmaPeakTime(),    // sigma_start_tick
			       float(endWire),                          // end_wire
			       0.,                                      // sigma_end_wire,
			       clusterHits.back()->PeakTime(),          // end_tick
			       clusterHits.back()->SigmaPeakTime(),     // sigma_end_tick
			       clusters->size(),                        // ID
			       clusterHits.front()->View(),             // view
			       clusterHits.front()->WireID().planeID(), // plane
			       recob::Cluster::Sentry                   // sentry
			       );

	clusters->emplace_back(cluster.move());

	// Associate the hits to this cluster
	util::CreateAssn(*this, evt, *(clusters.get()), clusterHits, *(associations.get()));

      } // End this cluster

    } // End loop over all clusters

  }

  evt.put(std::move(clusters));
  evt.put(std::move(associations));

  return;
    
}

DEFINE_ART_MODULE(cluster::BlurredCluster)
