////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgTrackMaker
// Module Type: producer
// File:        PMAlgTrackMaker_module.cc
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
//
// Creates 3D tracks using Projection Matching Algorithm, see RecoAlg/ProjectionMatchingAlg.h
// for basics of the PMA algorithm and its possible settings.
//
// Progress:
//    May-June 2015:   track finding and validation, growing tracks by iterative merging of matching
//                     clusters, no attempts to build multi-track structures, however cosmic tracking
//                     works fine as they are sets of independent tracks
//    June-July 2015:  merging track parts within a single tpc and stitching tracks across tpc's
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"

#include "MCCheater/BackTracker.h"
#include "Simulation/ParticleList.h"
#include "SimulationBase/MCParticle.h"

#include "RecoAlg/ProjectionMatchingAlg.h"
#include "RecoAlg/PMAlg/PmaTrack3D.h"
#include "RecoAlg/PMAlg/Utilities.h"

#include "TTree.h"

#include <memory>

namespace trkf {

typedef std::map< size_t, std::vector<double> > dedx_map;
typedef std::map< unsigned int, std::vector< art::Ptr<recob::Hit> > > view_hitmap;
typedef std::map< unsigned int, view_hitmap > tpc_view_hitmap;
typedef std::map< unsigned int, tpc_view_hitmap > cryo_tpc_view_hitmap;

typedef std::map< size_t, std::vector<pma::Track3D*> > tpc_track_map;

struct TrkCandidate {
	pma::Track3D* Track;
	std::vector< size_t > Clusters;
	double Mse, Validation;
	bool Good;
};

class PMAlgTrackMaker : public art::EDProducer {
public:
  explicit PMAlgTrackMaker(fhicl::ParameterSet const & p);

  PMAlgTrackMaker(PMAlgTrackMaker const &) = delete;
  PMAlgTrackMaker(PMAlgTrackMaker &&) = delete;
  PMAlgTrackMaker & operator = (PMAlgTrackMaker const &) = delete;
  PMAlgTrackMaker & operator = (PMAlgTrackMaker &&) = delete;

  void beginJob() override;

  void reconfigure(fhicl::ParameterSet const& p) override;

  void produce(art::Event & e) override;


private:
  // *** methods to create tracks from clusters ***

  // loop over all clusters and build as much as possible to find
  // the logic implemented here for sure is not exhaustive, it was
  // checked on long cosmic tracks and low energy stopping tracks
  // and seems to be a good example how to use the algorithm
  int fromMaxCluster(const art::Event& evt, std::vector< pma::Track3D* >& result);

  void fromMaxCluster_tpc(
    std::vector< pma::Track3D* >& result,
    art::Handle< std::vector<recob::Cluster> > clusters,
    const art::FindManyP< recob::Hit >& fbp,
    size_t minBuildSize, unsigned int tpc, unsigned int cryo);

  bool extendTrack(
	TrkCandidate& candidate,
	const std::vector< art::Ptr<recob::Hit> >& hits,
	unsigned int testView, bool add_nodes);

  int matchCluster(
    const TrkCandidate& trk,
	art::Handle< std::vector<recob::Cluster> > clusters,
	const art::FindManyP< recob::Hit >& fbp,
	size_t minSize, double fraction,
	unsigned int preferedView, unsigned int testView,
	unsigned int tpc, unsigned int cryo);

  // display what was used and what is left, just for development and debugging
  void listUsedClusters(art::Handle< std::vector<recob::Cluster> > clusters) const;
  // ------------------------------------------------------

  // build tracks from clusters associated by any other module - not yet implemented
  int fromExistingAssocs(const art::Event& evt, std::vector< pma::Track3D* >& result);
  // ------------------------------------------------------


  // ************* some common functionality **************
  cryo_tpc_view_hitmap c_t_v_hits;
  bool sortHits(const art::Event& evt);

  std::vector< size_t > used_clusters, initial_clusters;
  std::map< unsigned int, std::vector<size_t> > tried_clusters;
  bool has(const std::vector<size_t>& v, size_t idx) const
  {
  	for (auto c : v) if (c == idx) return true;
  	return false;
  }

  std::vector< TrkCandidate > fCandidates;

  int maxCluster(
    art::Handle< std::vector<recob::Cluster> > clusters,
    const art::FindManyP< recob::Hit >& fbp,
    size_t min_clu_size,
    geo::View_t view, unsigned int tpc, unsigned int cryo);
  int maxCluster(
	size_t first_idx,
    art::Handle< std::vector<recob::Cluster> > clusters,
    const art::FindManyP< recob::Hit >& fbp,
    float tmin, float tmax, size_t min_clu_size,
    geo::View_t view, unsigned int tpc, unsigned int cryo);

  bool areCoLinear(
	pma::Track3D* trk1, pma::Track3D* trk2,
	double& dist, double& cos, bool& reverseOrder,
	double distThr, double distThrMin,
	double distProjThr,
	double cosThr);
  void mergeCoLinear(std::vector< pma::Track3D* >& tracks);
  void mergeCoLinear(tpc_track_map& tracks);

  bool reassignHits(
	const std::vector< art::Ptr<recob::Hit> >& hits,
	std::vector< pma::Track3D* >& tracks, size_t trk_idx);
  bool reassignSingleViewEnds(std::vector< pma::Track3D* >& tracks);

  double validate(pma::Track3D& trk, unsigned int testView);
  recob::Track convertFrom(const pma::Track3D& src);

  bool isMcStopping(void) const; // to be moved to the testing module
  int getMcPdg(const pma::Track3D& trk) const; // to be moved to the testing module
  // ------------------------------------------------------

  art::ServiceHandle< geo::Geometry > fGeom;
  art::ServiceHandle<util::DetectorProperties> fDetProp;

  // ******************* tree output **********************
  int fEvNumber;        // event number
  int fTrkIndex;        // track index in the event, same for all dQ/dx points of the track
  int fPlaneIndex;      // wire plane index of the dQ/dx data point
  int fPlaneNPoints;    // number of data points in this plane
  int fIsStopping;      // tag tracks of stopping particles
  int fMcPdg;           // MC truth PDG matched for a track
  int fPidTag;          // Tag: 0=trk-like, 1=cascade-like
  double fdQdx;         // dQ/dx data point stored for each plane
  double fRange;        // residual range at dQ/dx data point, tracks are auto-flipped
  double fLength;       // track length
  double fHitsMse;      // MSE of hits: mean dist^2 of hit to 2D track projection
  double fSegAngMean;   // Mean segment-segment 3D angle.
  TTree* fTree_dQdx;    // dQ/dx info (saved optionally)
  TTree* fTree_trk;     // overall info

  // ******************** parameters **********************
  std::string fHitModuleLabel; // label for hits collection (used for trk validation)
  std::string fCluModuleLabel; // label for input cluster collection
  int fCluMatchingAlg;         // which algorithm for cluster association

  size_t fMinSeedSize1stPass;  // min. cluster size used to start building a track in the 1st pass
  size_t fMinSeedSize2ndPass;  // min. cluster size used to start building a track in the 2nd pass

  bool fFlipToBeam;            // set the track direction to increasing Z values
  bool fFlipDownward;          // set the track direction to decreasing Y values
  bool fAutoFlip_dQdx;         // set the track direction to increasing dQ/dx
  bool fSave_dQdx;             // for debugging purposes, off by default

  bool fMergeWithinTPC;          // merge witnin single TPC; finds tracks best matching by angle, with limits:
  double fMergeTransverseShift;  //   - max. transverse displacement [cm] between tracks
  double fMergeAngle;            //   - max. angle [degree] between tracks (nearest segments)

  bool fStitchBetweenTPCs;       // stitch between TPCs; finds tracks best matching by angle, with limits:
  double fStitchDistToWall;      //   - max. track endpoint distance [cm] to TPC boundary
  double fStitchTransverseShift; //   - max. transverse displacement [cm] between tracks
  double fStitchAngle;           //   - max. angle [degree] between tracks (nearest segments)

  pma::ProjectionMatchingAlg fProjectionMatchingAlg;
  double fMinTwoViewFraction;  // ProjectionMatchingAlg parameter used also in the module
};
// ------------------------------------------------------

PMAlgTrackMaker::PMAlgTrackMaker(fhicl::ParameterSet const & p) :
	fProjectionMatchingAlg(p.get< fhicl::ParameterSet >("ProjectionMatchingAlg"))
{
	this->reconfigure(p);
	produces< std::vector<recob::Track> >();
	produces< std::vector<recob::SpacePoint> >();
	produces< art::Assns<recob::Track, recob::Hit> >();
	produces< art::Assns<recob::Track, recob::SpacePoint> >();
	produces< art::Assns<recob::SpacePoint, recob::Hit> >();
}
// ------------------------------------------------------

void PMAlgTrackMaker::beginJob()
{
	art::ServiceHandle<art::TFileService> tfs;
	fTree_dQdx = tfs->make<TTree>("PMAlgTrackMaker_dQdx", "tracks dQ/dx info");
	fTree_dQdx->Branch("fEvNumber", &fEvNumber, "fEvNumber/I");
	fTree_dQdx->Branch("fTrkIndex", &fTrkIndex, "fTrkIndex/I");
	fTree_dQdx->Branch("fPlaneIndex", &fPlaneIndex, "fPlaneIndex/I");
	fTree_dQdx->Branch("fPlaneNPoints", &fPlaneNPoints, "fPlaneNPoints/I");
	fTree_dQdx->Branch("fIsStopping", &fIsStopping, "fIsStopping/I");
	fTree_dQdx->Branch("fMcPdg", &fMcPdg, "fMcPdg/I");
	fTree_dQdx->Branch("fdQdx", &fdQdx, "fdQdx/D");
	fTree_dQdx->Branch("fRange", &fRange, "fRange/D");
	fTree_dQdx->Branch("fLength", &fLength, "fLength/D");

	fTree_trk = tfs->make<TTree>("PMAlgTrackMaker_trk", "tracks overall info");
	fTree_trk->Branch("fEvNumber", &fEvNumber, "fEvNumber/I");
	fTree_trk->Branch("fTrkIndex", &fTrkIndex, "fTrkIndex/I");
	fTree_trk->Branch("fMcPdg", &fMcPdg, "fMcPdg/I");
	fTree_trk->Branch("fLength", &fLength, "fLength/D");
	fTree_trk->Branch("fHitsMse", &fHitsMse, "fHitsMse/D");
	fTree_trk->Branch("fSegAngMean", &fSegAngMean, "fSegAngMean/D");
	fTree_trk->Branch("fPidTag", &fPidTag, "fPidTag/I");
}

void PMAlgTrackMaker::reconfigure(fhicl::ParameterSet const& pset)
{
	fHitModuleLabel = pset.get< std::string >("HitModuleLabel");
	fCluModuleLabel = pset.get< std::string >("ClusterModuleLabel");
	fCluMatchingAlg = pset.get< int >("CluMatchingAlg");

	fMinSeedSize1stPass = pset.get< size_t >("MinSeedSize1stPass");
	fMinSeedSize2ndPass = pset.get< size_t >("MinSeedSize2ndPass");

	fFlipToBeam = pset.get< bool >("FlipToBeam");
	fFlipDownward = pset.get< bool >("FlipDownward");
	fAutoFlip_dQdx = pset.get< bool >("AutoFlip_dQdx");
	fSave_dQdx = pset.get< bool >("Save_dQdx");

	fMergeWithinTPC = pset.get< bool >("MergeWithinTPC");
	fMergeTransverseShift = pset.get< double >("MergeTransverseShift");
	fMergeAngle = pset.get< double >("MergeAngle");

	fStitchBetweenTPCs = pset.get< bool >("StitchBetweenTPCs");
	fStitchDistToWall = pset.get< double >("StitchDistToWall");
	fStitchTransverseShift = pset.get< double >("StitchTransverseShift");
	fStitchAngle = pset.get< double >("StitchAngle");

	fProjectionMatchingAlg.reconfigure(pset.get< fhicl::ParameterSet >("ProjectionMatchingAlg"));
	fMinTwoViewFraction = pset.get< double >("ProjectionMatchingAlg.MinTwoViewFraction");
}
// ------------------------------------------------------

recob::Track PMAlgTrackMaker::convertFrom(const pma::Track3D& src)
{
	std::vector< TVector3 > xyz, dircos;
	std::vector< std::vector<double> > dst_dQdx; // [view][dQ/dx]
	dst_dQdx.push_back(std::vector<double>()); // kU
	dst_dQdx.push_back(std::vector<double>()); // kV
	dst_dQdx.push_back(std::vector<double>()); // kZ

	unsigned int cryo = src.FrontCryo();
	unsigned int tpc = src.FrontTPC();

	std::map< unsigned int, dedx_map > src_dQdx;
	if (fGeom->TPC(tpc, cryo).HasPlane(geo::kU))
	{
		src_dQdx[geo::kU] = dedx_map();
		src.GetRawdEdxSequence(src_dQdx[geo::kU], geo::kU);
	}
	if (fGeom->TPC(tpc, cryo).HasPlane(geo::kV))
	{
		src_dQdx[geo::kV] = dedx_map();
		src.GetRawdEdxSequence(src_dQdx[geo::kV], geo::kV);
	}
	if (fGeom->TPC(tpc, cryo).HasPlane(geo::kZ))
	{
		src_dQdx[geo::kZ] = dedx_map();
		src.GetRawdEdxSequence(src_dQdx[geo::kZ], geo::kZ);
	}

	for (size_t i = 0; i < src.size(); i++)
	{
		xyz.push_back(src[i]->Point3D());

		if (i < src.size() - 1)
		{
			TVector3 dc(src[i + 1]->Point3D());
			dc -= src[i]->Point3D();
			dc *= 1.0 / dc.Mag();
			dircos.push_back(dc);
		}
		else dircos.push_back(dircos.back());

		double dQ = 0., dx = 0.;
		dst_dQdx[geo::kU].push_back(0.);
		dst_dQdx[geo::kV].push_back(0.);
		dst_dQdx[geo::kZ].push_back(0.);

		double dQdx;
		for (auto const& m : src_dQdx)
		{
			auto it = m.second.find(i);
			if (it != m.second.end())
			{
				dQ = it->second[5];
				dx = it->second[6];
				if (dx > 0.) dQdx = dQ/dx;
				else dQdx = 0.;

				dst_dQdx[m.first][i] = dQdx;

				break;
			}
		}
	}

	fLength = src.Length();
	fHitsMse = src.GetMse();
	fSegAngMean = src.GetMeanAng();
	fMcPdg = getMcPdg(src);

	fPidTag = 0;                 // 0 is track-like (long and/or very straight, well matching 2D hits); 0x10000 is EM shower-like trajectory
	if ( // (fLength < 80.0) &&  // tag only short tracks as EM shower-like
	    ((fHitsMse > 0.0001 * fLength) || (fSegAngMean < 3.0)))
		fPidTag = 0x10000;

	if (fSave_dQdx)
	{
		fIsStopping = (int)isMcStopping();
		for (unsigned int view = 0; view < fGeom->Nviews(); view++)
			if (fGeom->TPC(tpc, cryo).HasPlane(view))
			{
				fPlaneIndex = view;
				fPlaneNPoints = src_dQdx[view].size();
				for (auto const& data : src_dQdx[view])
				{
					double dQ = data.second[5];
					double dx = data.second[6];
					fRange = data.second[7];
					if (dx > 0.)
					{
						fdQdx = dQ/dx;
						fTree_dQdx->Fill();
					}
				}
			}
	}
	fTree_trk->Fill();

	if (xyz.size() != dircos.size())
	{
		mf::LogError("PMAlgTrackMaker") << "pma::Track3D to recob::Track conversion problem.";
	}
	return recob::Track(xyz, dircos, dst_dQdx, std::vector< double >(2, util::kBogusD), fTrkIndex + fPidTag);
}
// ------------------------------------------------------

bool PMAlgTrackMaker::isMcStopping(void) const
{
	art::ServiceHandle< cheat::BackTracker > bt;
	const sim::ParticleList& plist = bt->ParticleList();
	const simb::MCParticle* particle = plist.Primary(0);

	if (particle)
	{
//		std::cout << "...:SIM:... " << particle->EndProcess()
//			<< " n:" << particle->NumberDaughters()
//			<< " m:" << particle->Mass()
//			<< " E:" << particle->EndE() << std::endl;
		return (particle->EndE() - particle->Mass() < 0.001);
		//return (particle->NumberDaughters() == 0);
	}
	else return false;
}

int PMAlgTrackMaker::getMcPdg(const pma::Track3D& trk) const
{
	// return 0; // would be good in case of real data...

	art::ServiceHandle< cheat::BackTracker > bt;
	std::map< int, size_t > pdg_counts;
	for (size_t i = 0; i < trk.size(); i++)
		if (trk[i]->View2D() == geo::kZ)
	{
		art::Ptr< recob::Hit > hit = trk[i]->Hit2DPtr();
		std::vector< sim::TrackIDE > tids = bt->HitToTrackID(hit);
		for (auto const& tid : tids)
		{
			const simb::MCParticle* pi = bt->TrackIDToParticle(tid.trackID);
			int pdg = pi->PdgCode();
			auto it = pdg_counts.find(pdg);
			if (it != pdg_counts.end()) pdg_counts[pdg]++;
			else pdg_counts[pdg] = 1;
        }
	}

	int best_pdg = 0;
	size_t count, max_count = 0;
	for (auto const& p : pdg_counts)
	{
		count = p.second;
		if (count > max_count)
		{
			max_count = count; best_pdg = p.first;
		}
	}
	return best_pdg;
}
// ------------------------------------------------------

double PMAlgTrackMaker::validate(pma::Track3D& trk, unsigned int testView)
{
	if ((trk.FirstElement()->GetDistToWall() < -3.0) ||
	    (trk.LastElement()->GetDistToWall() < -3.0))
	{
		mf::LogVerbatim("PMAlgTrackMaker") << "first or last node too far out of its initial TPC";
		return 0.0;
	}

	if (testView != geo::kUnknown)
		mf::LogVerbatim("PMAlgTrackMaker") << "validation in plane: " << testView;
	else return 1.0;

	std::vector< art::Ptr<recob::Hit> >& hits = c_t_v_hits[trk.FrontCryo()][trk.FrontTPC()][testView];

	// always validate (needed for disambiguation postponed to 3D step):
	return fProjectionMatchingAlg.validate(trk, hits, testView);

	// in case of usual disambig on hit level one may validate only if there are at least a few hits:
/*	if (hits.size() > 10) return fProjectionMatchingAlg.validate(trk, hits, testView);
	else
	{
		mf::LogWarning("PMAlgTrackMaker") << "   too few hits (" << hits.size() << ")";
		return 1.0;
	}
*/
}
// ------------------------------------------------------

bool PMAlgTrackMaker::extendTrack(TrkCandidate& candidate,
	const std::vector< art::Ptr<recob::Hit> >& hits,
	unsigned int testView, bool add_nodes)
{
	double m_max = 2.0 * candidate.Mse; // max acceptable MSE value
	if (m_max < 0.05) m_max = 0.05;     // this is still good, low MSE value

	double v_min1 = 0.98 * candidate.Validation;
	double v_min2 = 0.9 * candidate.Validation;

	pma::Track3D* copy = fProjectionMatchingAlg.extendTrack(*(candidate.Track), hits, add_nodes);
	double m1 = copy->GetMse();
	double v1 = validate(*copy, testView);

	if (((m1 < candidate.Mse) && (v1 >= v_min2)) ||
	    ((m1 < 0.5) && (m1 <= m_max) && (v1 >= v_min1)))
	{
		mf::LogVerbatim("PMAlgTrackMaker")
			<< "  track EXTENDED, MSE = " << m1 << ", v = " << v1;
		delete candidate.Track;    // delete previous version of the track
		candidate.Track = copy;    // replace with the new track
		copy->SortHits();          // sort hits in the new track

		candidate.Mse = m1;        // save info
		candidate.Validation = v1;

		return true;
	}
	else
	{
		mf::LogVerbatim("PMAlgTrackMaker")
			<< "  track NOT extended, MSE = " << m1 << ", v = " << v1;
		delete copy;
		return false;
	}
}
// ------------------------------------------------------

bool PMAlgTrackMaker::areCoLinear(pma::Track3D* trk1, pma::Track3D* trk2,
	double& dist, double& cos, bool& reverseOrder,
	double distThr, double distThrMin,
	double distProjThr,
	double cosThr)
{
	double lmax;
	double l1 = trk1->Length();
	double l2 = trk2->Length();

	if (l1 > l2) lmax = l1;
	else lmax = l2;

	double d = lmax * distThr;
	if (d < distThrMin) d = distThrMin;

	unsigned int k = 0;
	double distFF = pma::Dist2(trk1->front()->Point3D(), trk2->front()->Point3D());
	dist = distFF;

	double distFB = pma::Dist2(trk1->front()->Point3D(), trk2->back()->Point3D());
	if (distFB < dist) { k = 1; dist = distFB; }

	double distBF = pma::Dist2(trk1->back()->Point3D(), trk2->front()->Point3D());
	if (distBF < dist) { k = 2; dist = distBF; }

	double distBB = pma::Dist2(trk1->back()->Point3D(), trk2->back()->Point3D());
	if (distBB < dist) { k = 3; dist = distBB; }

	dist = sqrt(dist);
	cos = 0.0;

	//std::cout << "  min dist:" << dist << " d:" << d << ", trk len:" << lmax << std::endl;
	if (dist < d)
	{
		pma::Track3D* tmp = 0;
		switch (k) // swap or flip to get trk1 end before trk2 start
		{
			case 0:	trk1->Flip(); break;
			case 1: tmp = trk1;	trk1 = trk2; trk2 = tmp; break;
			case 2: break;
			case 3: trk2->Flip(); break;
			default: mf::LogError("PMAlgTrackMaker") << "Should never happen.";
		}
		if (k == 1) reverseOrder = true;
		else reverseOrder = false;

		size_t nodeEndIdx = trk1->Nodes().size() - 1;

		double distProj1 = pma::Dist2(
			trk1->back()->Point3D(),
			pma::GetProjectionToSegment(trk1->back()->Point3D(),
				trk2->Nodes()[0]->Point3D(), trk2->Nodes()[1]->Point3D()));

		double distProj2 = pma::Dist2(
			trk2->front()->Point3D(),
			pma::GetProjectionToSegment(trk2->front()->Point3D(),
				trk1->Nodes()[nodeEndIdx - 1]->Point3D(), trk1->Nodes()[nodeEndIdx]->Point3D()));

		TVector3 dir1 = trk1->Nodes()[nodeEndIdx]->Point3D() - trk1->Nodes()[nodeEndIdx - 1]->Point3D();
		TVector3 dir2 = trk2->Nodes()[1]->Point3D() - trk2->Nodes()[0]->Point3D();

		cos = (dir1 * dir2) / (dir1.Mag() * dir2.Mag());

		//std::cout << "     cos:" << cos << " p1:" << distProj1 << " p2:" << distProj2 << std::endl;
		if ((cos > cosThr) && (distProj1 < distProjThr) && (distProj2 < distProjThr))
			return true;
	}
	return false;
}
// ------------------------------------------------------

void PMAlgTrackMaker::mergeCoLinear(std::vector< pma::Track3D* >& tracks)
{
	double distThr = 0.05;    // max gap as a fraction of the longer track length
	double distThrMin = 0.5;  // lower limit of max gap threshold [cm]

	double distProjThr = fMergeTransverseShift;
	double cosThr = cos(fMergeAngle);

	bool r;
	double d, c;
	size_t t = 0, u = 0;
	while (t < tracks.size())
	{
		pma::Track3D* trk1 = tracks[t];

		pma::Track3D* trk2 = 0;
		for (u = t + 1; u < tracks.size(); u++)
		{
			trk2 = tracks[u];

			if (areCoLinear(trk1, trk2, d, c, r, distThr, distThrMin, distProjThr, cosThr)) break;

			trk2 = 0;
		}

		if (trk2)
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "Merge track ("
				<< trk1->size() << ") with track (" << trk2->size() << ")";
			if (r)
			{
				fProjectionMatchingAlg.mergeTracks(*trk2, *trk1, true);
				tracks[t] = trk2;
				trk2 = trk1;
			}
			else
			{
				fProjectionMatchingAlg.mergeTracks(*trk1, *trk2, true);
			}

			tracks.erase(tracks.begin() + u);
			delete trk2;
		}
		else t++;
	}
}
// ------------------------------------------------------

void PMAlgTrackMaker::mergeCoLinear(tpc_track_map& tracks)
{
	double distThr = 0.25;    // max gap as a fraction of the longer track length
	double distThrMin = 0.5;  // lower limit of max gap threshold [cm]

	double distProjThr = fStitchTransverseShift;
	double cosThr = cos(fStitchAngle);

	double wallDistThr = fStitchDistToWall;
	double dfront1, dback1, dfront2, dback2;

	for (auto & tpc_entry1 : tracks)
	{
		unsigned int tpc1 = tpc_entry1.first;
		std::vector< pma::Track3D* >& tracks1 = tpc_entry1.second;

		size_t t = 0;
		while (t < tracks1.size())
		{
			bool r, reverse = false;
			double d, c, cmax = 0.0;
			pma::Track3D* best_trk2 = 0;
			unsigned int best_tpc = 0;
			size_t best_idx = 0;

			pma::Track3D* trk1 = tracks1[t];
			dfront1 = trk1->Nodes().front()->GetDistToWall();
			dback1 = trk1->Nodes().back()->GetDistToWall();
			if ((dfront1 < wallDistThr) || (dback1 < wallDistThr))
			{
				for (auto & tpc_entry2 : tracks)
				{
					unsigned int tpc2 = tpc_entry2.first;
					if (tpc2 == tpc1) continue;

					std::vector< pma::Track3D* >& tracks2 = tpc_entry2.second;

					for (size_t u = 0; u < tracks2.size(); u++)
					{
						pma::Track3D* trk2 = tracks2[u];
						dfront2 = trk2->Nodes().front()->GetDistToWall();
						dback2 = trk2->Nodes().back()->GetDistToWall();
						if ((dfront2 < wallDistThr) || (dback2 < wallDistThr))
						{
							if (areCoLinear(trk1, trk2, d, c, r, distThr, distThrMin, distProjThr, cosThr))
							{
								if (c > cmax)
								{
									cmax = c;
									best_trk2 = trk2;
									best_tpc = tpc2;
									best_idx = u;
									reverse = r;
								}
							}
						}
					}
				}
			}

			if (best_trk2)
			{
				mf::LogVerbatim("PMAlgTrackMaker") << "Merge track ("
					<< tpc1 << ":" << tracks1.size() << ":" << trk1->size() << ") with track ("
					<< best_tpc  << ":" << tracks[best_tpc].size() << ":" << best_trk2->size() << ")";
				if (reverse)
				{
					fProjectionMatchingAlg.mergeTracks(*best_trk2, *trk1, true);
					tracks1[t] = best_trk2;
					best_trk2 = trk1;
				}
				else
				{
					fProjectionMatchingAlg.mergeTracks(*trk1, *best_trk2, true);
				}

				tracks[best_tpc].erase(tracks[best_tpc].begin() + best_idx);
				delete best_trk2;
			}
			else t++;
		}
	}
}
// ------------------------------------------------------

bool PMAlgTrackMaker::reassignHits(
	const std::vector< art::Ptr<recob::Hit> >& hits,
	std::vector< pma::Track3D* >& tracks, size_t trk_idx)
{
	pma::Track3D* trk1 = tracks[trk_idx];
	pma::Track3D* best_trk = 0;

	size_t best_u = 0, n_max = 0;
	for (size_t u = 0; u < tracks.size(); u++)
		if (trk_idx != u)
	{
		pma::Track3D* trk2 = tracks[u];
		size_t n = fProjectionMatchingAlg.testHits(*trk2, hits);
		if (n > n_max) { n_max = n; best_u = u; best_trk = trk2; }
	}

	if (best_trk && (n_max >= hits.size() / 3)) // /2
	{
		mf::LogVerbatim("PMAlgTrackMaker") << "  Reassign " << n_max << " hits." << std::endl;

		trk1->RemoveHits(hits);
		trk1->CleanupTails();
		trk1->ShiftEndsToHits();

		pma::Track3D* ext = fProjectionMatchingAlg.extendTrack(*best_trk, hits,	false);
		ext->SortHits(); ext->ShiftEndsToHits();
		if (fProjectionMatchingAlg.isContained(*ext))
		{
			tracks[best_u] = ext;
			delete best_trk;
			return true;
		}
		else delete ext;
	}
	return false;
}

bool PMAlgTrackMaker::reassignSingleViewEnds(std::vector< pma::Track3D* >& tracks)
{
	bool result = false;
	for (size_t t = 0; t < tracks.size(); t++)
	{
		pma::Track3D* trk = tracks[t];
		if (trk->size() < 5) continue;

		trk->DisableSingleViewEnds();

		std::vector< art::Ptr<recob::Hit> > hits;

		size_t idx = 0;
		while ((idx < trk->size() - 1) && !((*trk)[idx]->IsEnabled()))
		{
			hits.push_back((*trk)[idx++]->Hit2DPtr());
		}

		double d2;
		if (idx > 0)
		{
			if ((idx < trk->size() - 1) &&
			    ((*trk)[idx]->View2D() == (*trk)[idx - 1]->View2D()))
			{
				double dprev = pma::Dist2((*trk)[idx]->Point3D(), (*trk)[idx - 1]->Point3D());
				double dnext = pma::Dist2((*trk)[idx]->Point3D(), (*trk)[idx + 1]->Point3D());
				if (dprev < dnext)
				{
					hits.push_back((*trk)[idx++]->Hit2DPtr());
				}
			}
			d2 = pma::Dist2((*trk)[idx]->Point3D(), (*trk)[idx - 1]->Point3D());
		}
		else d2 = 0.0;

		if ((hits.size() > 1) || (d2 > 1.0)) // min. 2 hits or single hit separated from the rest
		{
			result |= reassignHits(hits, tracks, t);
		}

		hits.clear();
		idx = trk->size() - 1;
		while ((idx > 0) && !((*trk)[idx]->IsEnabled()))
		{
			hits.push_back((*trk)[idx--]->Hit2DPtr());
		}

		if (idx < trk->size() - 1)
		{
			if ((idx > 0) &&
			    ((*trk)[idx]->View2D() == (*trk)[idx + 1]->View2D()))
			{
				double dprev = pma::Dist2((*trk)[idx]->Point3D(), (*trk)[idx + 1]->Point3D());
				double dnext = pma::Dist2((*trk)[idx]->Point3D(), (*trk)[idx - 1]->Point3D());
				if (dprev < dnext)
				{
					hits.push_back((*trk)[idx--]->Hit2DPtr());
				}
			}
			d2 = pma::Dist2((*trk)[idx]->Point3D(), (*trk)[idx + 1]->Point3D());
		}
		else d2 = 0.0;

		if ((hits.size() > 1) || (d2 > 1.0))  // min. 2 hits or single hit separated from the rest
		{
			result |= reassignHits(hits, tracks, t);
		}

		trk->SelectHits();
	}
	return result;
}
// ------------------------------------------------------

bool PMAlgTrackMaker::sortHits(const art::Event& evt)
{
	c_t_v_hits.clear();

	art::Handle< std::vector<recob::Hit> > hitListHandle;
	std::vector< art::Ptr<recob::Hit> > hitlist;
	if (evt.getByLabel(fHitModuleLabel, hitListHandle))
	{
		art::fill_ptr_vector(hitlist, hitListHandle);

		unsigned int cryo, tpc, view;
		for (auto const& h : hitlist)
		{
			cryo = h->WireID().Cryostat;
			tpc = h->WireID().TPC;
			view = h->WireID().Plane;

			c_t_v_hits[cryo][tpc][view].push_back(h);
		}
		return true;
	}
	else return false;
}
// ------------------------------------------------------

void PMAlgTrackMaker::produce(art::Event& evt)
{
	fEvNumber = evt.id().event();

	std::vector< pma::Track3D* > result;

	std::unique_ptr< std::vector< recob::Track > > tracks(new std::vector< recob::Track >);
	std::unique_ptr< std::vector< recob::SpacePoint > > allsp(new std::vector< recob::SpacePoint >);

	std::unique_ptr< art::Assns< recob::Track, recob::Hit > > trk2hit(new art::Assns< recob::Track, recob::Hit >);
	std::unique_ptr< art::Assns< recob::Track, recob::SpacePoint > > trk2sp(new art::Assns< recob::Track, recob::SpacePoint >);
	std::unique_ptr< art::Assns< recob::SpacePoint, recob::Hit > > sp2hit(new art::Assns< recob::SpacePoint, recob::Hit >);

	if (sortHits(evt))
	{
		int retCode = 0;
		switch (fCluMatchingAlg)
		{
			default:
			case 1: retCode = fromMaxCluster(evt, result); break;

			case 2: retCode = fromExistingAssocs(evt, result); break;
		}
		switch (retCode)
		{
			case -2: mf::LogError("Summary") << "problem"; break;
			case -1: mf::LogWarning("Summary") << "no input"; break;
			case  0: mf::LogVerbatim("Summary") << "no tracks done"; break;
			default:
				if (retCode < 0) mf::LogVerbatim("Summary") << "unknown result";
				else if (retCode == 1) mf::LogVerbatim("Summary") << retCode << " track ready";
				else mf::LogVerbatim("Summary") << retCode << " tracks ready";
				break;
		}

		if (result.size()) // ok, there is something to save
		{
			size_t spStart = 0, spEnd = 0;
			double sp_pos[3], sp_err[6];
			for (size_t i = 0; i < 6; i++) sp_err[i] = 1.0;

			double dQdxFlipThr = 0.0;
			if (fFlipToBeam) dQdxFlipThr = 0.4;

			fTrkIndex = 0;
			for (auto const& trk : result)
			{
				if (fFlipToBeam)    // flip the track to the beam direction
				{
					double z0 = trk->front()->Point3D().Z();
					double z1 = trk->back()->Point3D().Z();
					if (z0 > z1) trk->Flip();
				}
				if (fFlipDownward)  // flip the track to point downward
				{
					double y0 = trk->front()->Point3D().Y();
					double y1 = trk->back()->Point3D().Y();
					if (y0 < y1) trk->Flip();
				}
				if (fAutoFlip_dQdx) // flip the track by dQ/dx
					fProjectionMatchingAlg.autoFlip(*trk, pma::Track3D::kForward, dQdxFlipThr);
					/* test code: fProjectionMatchingAlg.autoFlip(*trk, pma::Track3D::kBackward, dQdxFlipThr); */

				tracks->push_back(convertFrom(*trk));
				fTrkIndex++;

				std::vector< art::Ptr< recob::Hit > > hits2d;
				art::PtrVector< recob::Hit > sp_hits;

				spStart = allsp->size();
				for (int h = trk->size() - 1; h >= 0; h--)
				{
					pma::Hit3D* h3d = (*trk)[h];
					hits2d.push_back(h3d->Hit2DPtr());

					if ((h == 0) ||
					      (sp_pos[0] != h3d->Point3D().X()) ||
					      (sp_pos[1] != h3d->Point3D().Y()) ||
					      (sp_pos[2] != h3d->Point3D().Z()))
					{
						if (sp_hits.size()) // hits assigned to the previous sp
						{
							util::CreateAssn(*this, evt, *allsp, sp_hits, *sp2hit);
							sp_hits.clear();
						}
						sp_pos[0] = h3d->Point3D().X();
						sp_pos[1] = h3d->Point3D().Y();
						sp_pos[2] = h3d->Point3D().Z();
						allsp->push_back(recob::SpacePoint(sp_pos, sp_err, 1.0));
					}
					sp_hits.push_back(h3d->Hit2DPtr());
				}
				if (sp_hits.size()) // hits assigned to the last sp
				{
					util::CreateAssn(*this, evt, *allsp, sp_hits, *sp2hit);
				}
				spEnd = allsp->size();

				if (hits2d.size())
				{
					util::CreateAssn(*this, evt, *tracks, *allsp, *trk2sp, spStart, spEnd);
					util::CreateAssn(*this, evt, *tracks, hits2d, *trk2hit);
				}
			}

			// data prods done, delete all pma::Track3D's
			for (size_t t = 0; t < result.size(); t++) delete result[t];
		}
	}
	else mf::LogError("PMAlgTrackMaker") << "Hits not found in the event.";

	evt.put(std::move(tracks));
	evt.put(std::move(allsp));

	evt.put(std::move(trk2hit));
	evt.put(std::move(trk2sp));
	evt.put(std::move(sp2hit));
}
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

int PMAlgTrackMaker::fromMaxCluster(const art::Event& evt, std::vector< pma::Track3D* >& result)
{
	art::Handle< std::vector<recob::Cluster> > cluListHandle;
	if (evt.getByLabel(fCluModuleLabel, cluListHandle))
	{
		initial_clusters.clear();
		tried_clusters.clear();
		used_clusters.clear();

		art::FindManyP< recob::Hit > fbp(cluListHandle, evt, fCluModuleLabel);

		tpc_track_map tracks; // track parts in tpc's

		// find reasonably large parts
		for (auto tpc_iter = fGeom->begin_TPC_id();
		          tpc_iter != fGeom->end_TPC_id();
		          tpc_iter++)
		{
			fromMaxCluster_tpc(tracks[tpc_iter->TPC], cluListHandle, fbp, fMinSeedSize1stPass, tpc_iter->TPC, tpc_iter->Cryostat);
		}

		// loop again to find small things
		for (auto tpc_iter = fGeom->begin_TPC_id();
		          tpc_iter != fGeom->end_TPC_id();
		          tpc_iter++)
		{
			fromMaxCluster_tpc(tracks[tpc_iter->TPC], cluListHandle, fbp, fMinSeedSize2ndPass, tpc_iter->TPC, tpc_iter->Cryostat);
		}

		// try correcting track ends:
		//   - single-view sections spuriously merged on 2D clusters level
		for (auto tpc_iter = fGeom->begin_TPC_id();
		          tpc_iter != fGeom->end_TPC_id();
		          tpc_iter++)
		{
			reassignSingleViewEnds(tracks[tpc_iter->TPC]);
		}

		// merge co-linear parts inside each tpc
		if (fMergeWithinTPC)
		{
			for (auto tpc_iter = fGeom->begin_TPC_id();
			          tpc_iter != fGeom->end_TPC_id();
			          tpc_iter++)
			{
				mf::LogVerbatim("PMAlgTrackMaker") << "Merge co-linear tracks within TPC " << tpc_iter->TPC << ".";
				mergeCoLinear(tracks[tpc_iter->TPC]);
			}
		}

		// merge co-linear parts between tpc's
		if (fStitchBetweenTPCs)
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "Stitch co-linear tracks between TPCs.";
			mergeCoLinear(tracks);
		}

		for (auto const & tpc_entry : tracks)
			for (auto const & trk : tpc_entry.second)
				result.push_back(trk);

		// used for development
		listUsedClusters(cluListHandle);
	}
	else
	{
		mf::LogWarning("PMAlgTrackMaker") << "no clusters";
		return -1;
	}

	return result.size();
}

int PMAlgTrackMaker::matchCluster(const TrkCandidate& trk,
	art::Handle< std::vector<recob::Cluster> > clusters,
	const art::FindManyP< recob::Hit >& fbp,
	size_t minSize, double fraction,
	unsigned int preferedView, unsigned int testView,
	unsigned int tpc, unsigned int cryo)
{
	double f, fmax = 0.0;
	unsigned int n, max = 0;
	int idx = -1;
	for (size_t i = 0; i < clusters->size(); ++i)
	{
		unsigned int view = (*clusters)[i].View();
		unsigned int nhits = (*clusters)[i].NHits();

		if (has(used_clusters, i) ||                             // don't try already used clusters
			has(trk.Clusters, i) ||                              // don't try clusters from this candidate
		    (view == testView) ||                                // don't use clusters from validation view
		    ((preferedView != geo::kUnknown)&&(view != preferedView)) || // only prefered view if specified
		    (nhits < minSize))                                   // skip small clusters
		    continue;

		std::vector< art::Ptr<recob::Hit> > v = fbp.at(i);
		n = fProjectionMatchingAlg.testHits(*(trk.Track), v);
		f = n / (double)v.size();
		if ((f > fraction) && (n > max))
		{
			max = n; fmax = f; idx = i;
		}
	}

	if (idx >= 0) mf::LogVerbatim("PMAlgTrackMaker") << "max matching hits: " << max << " (" << fmax << ")";
	else mf::LogVerbatim("PMAlgTrackMaker") << "no clusters to extend the track";

	return idx;
}


void PMAlgTrackMaker::fromMaxCluster_tpc(
	std::vector< pma::Track3D* >& result,
	art::Handle< std::vector<recob::Cluster> > clusters,
	const art::FindManyP< recob::Hit >& fbp,
	size_t minBuildSize, unsigned int tpc, unsigned int cryo)
{
	initial_clusters.clear();

	size_t minSizeCompl = minBuildSize / 8;  // smaller minimum required in complementary views
	if (minSizeCompl < 2) minSizeCompl = 2;  // but at least two hits!

	int max_first_idx = 0;
	while (max_first_idx >= 0) // loop over clusters, any view, starting from the largest
	{
		max_first_idx = maxCluster(clusters, fbp, minBuildSize, geo::kUnknown, tpc, cryo); // any view
		if (max_first_idx >= 0)
		{
			const recob::Cluster& clu_first = (*clusters)[max_first_idx];

			geo::View_t first_view = clu_first.View();
			geo::View_t sec_view_a, sec_view_b;
			switch (first_view)
			{
				case geo::kU: sec_view_a = geo::kZ; sec_view_b = geo::kV; break;
				case geo::kV: sec_view_a = geo::kZ; sec_view_b = geo::kU; break;
				case geo::kZ: sec_view_a = geo::kV; sec_view_b = geo::kU; break;
				default:
					mf::LogError("PMAlgTrackMaker") << "Not a 2D view.";
					return;
			}

			tried_clusters[geo::kU].clear();
			tried_clusters[geo::kV].clear();
			tried_clusters[geo::kZ].clear();

			tried_clusters[first_view].push_back(max_first_idx);
			initial_clusters.push_back(max_first_idx);

			unsigned int nFirstHits = clu_first.NHits();
			mf::LogVerbatim("PMAlgTrackMaker") << "use plane  *** " << first_view << " ***  size: " << nFirstHits;

			std::vector< art::Ptr<recob::Hit> > v_first = fbp.at(max_first_idx);

			float tmax = fDetProp->ConvertTicksToX(v_first.front()->PeakTime(), first_view, tpc, cryo);
			float t, tmin = tmax;
			for (size_t j = 1; j < v_first.size(); ++j)
			{
				t = fDetProp->ConvertTicksToX(v_first[j]->PeakTime(), first_view, tpc, cryo);
				if (t > tmax) { tmax = t; }
				if (t < tmin) { tmin = t; }
			}


			fCandidates.clear(); // possible solutions of the selected cluster and clusters in complementary views

			bool try_build = true;
			while (try_build) // loop over complementary views
			{
				TrkCandidate candidate;
				candidate.Clusters.push_back(max_first_idx);

				int idx, max_sec_a_idx, max_sec_b_idx;
				max_sec_a_idx = maxCluster(max_first_idx, clusters, fbp, tmin, tmax, minSizeCompl, sec_view_a, tpc, cryo);
				max_sec_b_idx = maxCluster(max_first_idx, clusters, fbp, tmin, tmax, minSizeCompl, sec_view_b, tpc, cryo);

				unsigned int nSecHitsA = 0, nSecHitsB = 0;
				if (max_sec_a_idx >= 0) nSecHitsA = (*clusters)[max_sec_a_idx].NHits();
				if (max_sec_b_idx >= 0) nSecHitsB = (*clusters)[max_sec_b_idx].NHits();

				unsigned int testView = geo::kUnknown;
				if ((nSecHitsA > nSecHitsB) && (nSecHitsA >= minSizeCompl))
				{
					mf::LogVerbatim("PMAlgTrackMaker") << "use plane  *** " << sec_view_a << " ***  size: " << nSecHitsA;
					tried_clusters[sec_view_a].push_back(max_sec_a_idx);
					idx = max_sec_a_idx;
					testView = sec_view_b;
				}
				else if (nSecHitsB >= minSizeCompl)
				{
					mf::LogVerbatim("PMAlgTrackMaker") << "use plane  *** " << sec_view_b << " ***  size: " << nSecHitsB;
					tried_clusters[sec_view_b].push_back(max_sec_b_idx);
					idx = max_sec_b_idx;
					testView = sec_view_a;
				}
				else try_build = false;

				if (try_build)
				{
					if (!fGeom->TPC(tpc, cryo).HasPlane(testView)) testView = geo::kUnknown;

					double m0 = 0.0, v0 = 0.0;
					double mseThr = 0.15, validThr = 0.7; // cuts for a good track candidate

					candidate.Clusters.push_back(idx);
					candidate.Track = fProjectionMatchingAlg.buildTrack(v_first, fbp.at(idx));

					if (candidate.Track && // no track if hits from 2 views do not alternate
					    fProjectionMatchingAlg.isContained(*(candidate.Track))) // sticks out of TPC's?
					{
						m0 = candidate.Track->GetMse();
						if (m0 < mseThr) // check validation only if MSE is passing - thanks for Tracy for noticing this
							v0 = validate(*(candidate.Track), testView);
					}
					if (candidate.Track && (m0 < mseThr) && (v0 > validThr)) // good candidate, try to extend it
					{
						mf::LogVerbatim("PMAlgTrackMaker")
							<< "  track candidate, MSE = " << m0 << ", v = " << v0;

						candidate.Mse = m0;
						candidate.Validation = v0;
						candidate.Good = true;

						size_t minSize = 5;      // min size for clusters matching
						double fraction = 0.4;   // min fraction of close hits

						idx = 0;
						while (idx >= 0) // try to collect matching clusters, use **any** plane except validation
						{
							idx = matchCluster(candidate, clusters, fbp, minSize, fraction, geo::kUnknown, testView, tpc, cryo);
							if (idx >= 0)
							{
								// try building extended copy:
								//                src,        hits,    valid.plane, add nodes
								if (extendTrack(candidate, fbp.at(idx),  testView,    true))
								{
									candidate.Clusters.push_back(idx);
								}
								else idx = -1;
							}
						}

						mf::LogVerbatim("PMAlgTrackMaker") << "merge clusters from the validation plane";
						fraction = 0.7; // only well matching the existing track

						idx = 0;
						while ((idx >= 0) && (testView != geo::kUnknown))
						{	//                     match clusters from the plane used previously for the validation
							idx = matchCluster(candidate, clusters, fbp, minSize, fraction, testView, geo::kUnknown, tpc, cryo);
							if (idx >= 0)
							{
								// no validation, no new nodes:
								if (extendTrack(candidate, fbp.at(idx), geo::kUnknown, false))
								{
									candidate.Clusters.push_back(idx);
								}
								else idx = -1;
							}
						}

						candidate.Validation = validate(*(candidate.Track), testView);
					}
					else
					{
						mf::LogVerbatim("PMAlgTrackMaker") << "track REJECTED, MSE = " << m0 << "; v = " << v0;
						candidate.Good = false; // save also bad matches to avoid trying again the same pair of clusters
					}
					fCandidates.push_back(candidate);
				}
				else
				{
					mf::LogVerbatim("PMAlgTrackMaker") << "no matching clusters";
				}
			} // end loop over complementary views

			if (fCandidates.size()) // save best candidate, release other tracks and clusters
			{
				int best_trk = -1;
				double f, max_f = 0., min_mse = 10., max_v = 0.;
				for (size_t t = 0; t < fCandidates.size(); t++)
					if (fCandidates[t].Good)
				{
					f = fProjectionMatchingAlg.twoViewFraction(*(fCandidates[t].Track));

					if ((f > max_f) || ((f == max_f) &&
						((fCandidates[t].Validation > max_v) || (fCandidates[t].Mse < min_mse))))
					{
						max_f = f;
						min_mse = fCandidates[t].Mse;
						max_v = fCandidates[t].Validation;
						best_trk = t;
					}
				}

				if ((best_trk > -1) && fCandidates[best_trk].Good && (max_f > fMinTwoViewFraction))
				{
					fCandidates[best_trk].Track->ShiftEndsToHits();

					result.push_back(fCandidates[best_trk].Track);

					for (auto c : fCandidates[best_trk].Clusters)
						used_clusters.push_back(c);
                }

				for (size_t t = 0; t < fCandidates.size(); t++)
				{
					if (int(t) != best_trk) delete fCandidates[t].Track;
				}
				fCandidates.clear();
			}
		}
		else
		{
			mf::LogWarning("PMAlgTrackMaker") << "small clusters only";
		}
	} // end loop over clusters, any view, from the largest
}
// ------------------------------------------------------

int PMAlgTrackMaker::maxCluster(
	art::Handle< std::vector<recob::Cluster> > clusters,
	const art::FindManyP< recob::Hit >& fbp,
	size_t min_clu_size,
	geo::View_t view, unsigned int tpc, unsigned int cryo)
{
	int idx = -1;
	size_t s_max = 0, s;

	for (size_t i = 0; i < clusters->size(); ++i)
	{
		if (has(used_clusters, i) ||
		    has(initial_clusters, i) ||
		    has(tried_clusters[view], i) ||
		   ((view != geo::kUnknown) && ((*clusters)[i].View() != view)))
		continue;

		std::vector< art::Ptr<recob::Hit> > v = fbp.at(i);

		if ((v.front()->WireID().TPC == tpc) &&
		    (v.front()->WireID().Cryostat == cryo))
		{
			s = v.size();
			if ((s >= min_clu_size) && (s > s_max))
			{
				s_max = s; idx = i;
			}
		}
	}
	return idx;
}
// ------------------------------------------------------

int PMAlgTrackMaker::maxCluster(size_t first_idx,
	art::Handle< std::vector<recob::Cluster> > clusters,
	const art::FindManyP< recob::Hit >& fbp,
	float tmin, float tmax, size_t min_clu_size,
	geo::View_t view, unsigned int tpc, unsigned int cryo)
{
	int idx = -1;
	size_t s_max = 0, s;
	double fraction = 0.0;
	float x;

	for (size_t i = 0; i < clusters->size(); ++i)
	{
		if (has(used_clusters, i) ||
		    has(initial_clusters, i) ||
		    has(tried_clusters[view], i) ||
		    ((*clusters)[i].NHits() <  min_clu_size) ||
		    ((*clusters)[i].View() != view)) continue;

		bool pair_checked = false;
		for (auto const & c : fCandidates)
			if (has(c.Clusters, first_idx) && has(c.Clusters, i))
			{
				pair_checked = true; break;
			}
		if (pair_checked) continue;
		    
		std::vector< art::Ptr<recob::Hit> > v = fbp.at(i);

		if ((v.front()->WireID().TPC == tpc) &&
		    (v.front()->WireID().Cryostat == cryo))
		{
			s = 0;
			for (size_t j = 0; j < v.size(); ++j)
			{
				x = fDetProp->ConvertTicksToX(v[j]->PeakTime(), view, tpc, cryo);
				if ((x >= tmin) && (x <= tmax)) s++;
			}

			if (s > s_max)
			{
				s_max = s; idx = i; fraction = s / (double)v.size();
			}
		}
	}
	if (fraction > 0.4) return idx;
	else return -1;
}
// ------------------------------------------------------

void PMAlgTrackMaker::listUsedClusters(art::Handle< std::vector<recob::Cluster> > clusters) const
{
	mf::LogVerbatim("PMAlgTrackMaker") << std::endl << "----------- matched clusters: -----------";
	for (size_t i = 0; i < clusters->size(); ++i)
		if (has(used_clusters, i))
		{
			mf::LogVerbatim("PMAlgTrackMaker")
				<< "    tpc: " << (*clusters)[i].Plane().TPC
				<< ";\tview: " << (*clusters)[i].View()
				<< ";\tsize: " << (*clusters)[i].NHits();
		}
	mf::LogVerbatim("PMAlgTrackMaker") << "--------- not matched clusters: ---------";
	for (size_t i = 0; i < clusters->size(); ++i)
		if (!has(used_clusters, i))
		{
			mf::LogVerbatim("PMAlgTrackMaker")
				<< "    tpc: " << (*clusters)[i].Plane().TPC
				<< ";\tview: " << (*clusters)[i].View()
				<< ";\tsize: " << (*clusters)[i].NHits();
		}
	mf::LogVerbatim("PMAlgTrackMaker") << "-----------------------------------------";
}
// ------------------------------------------------------
// ------------------------------------------------------

int PMAlgTrackMaker::fromExistingAssocs(const art::Event& evt, std::vector< pma::Track3D* >& result)
{
	return 0;
}
// ------------------------------------------------------
// ------------------------------------------------------

DEFINE_ART_MODULE(PMAlgTrackMaker)

} // namespace trkf

