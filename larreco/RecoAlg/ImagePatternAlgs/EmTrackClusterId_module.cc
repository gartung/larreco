/////////////////////////////////////////////////////////////////////////////////
// Class:       EmTrackClusterId
// Module Type: producer
// File:        EmTrackClusterId_module.cc
// Authors:     dorota.stefan@cern.ch pplonski86@gmail.com robert.sulej@cern.ch
//
// Module applies neural net to 2D image made of deconvoluted wire waveforms in
// order to distinguish EM-like activity from track-like objects. Clusters of
// hits that were recognized as EM-like event parts are produced. Module uses
// clusters made with any algorithm as input; optionally (recommended) also
// single, unclustered hits are tested and if they look like EM parts then
// single-hit clusters are produced and added to the output collection.
//
/////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg/PointIdAlg.h"

#include <memory>

namespace nnet {

class EmTrackClusterId : public art::EDProducer {
public:

	struct Config {
		using Name = fhicl::Name;
		using Comment = fhicl::Comment;

		fhicl::Table<nnet::PointIdAlg::Config> PointIdAlg {
			Name("PointIdAlg")
		};

		fhicl::Atom<art::InputTag> WireLabel {
			Name("WireLabel"),
			Comment("tag of deconvoluted ADC on wires (recob::Wire)")
		};

		fhicl::Atom<art::InputTag> HitModuleLabel {
			Name("HitModuleLabel"),
			Comment("tag of hits made to create clusters, unclustered leftovers to be EM/track tagged")
		};

		fhicl::Atom<art::InputTag> ClusterModuleLabel {
			Name("ClusterModuleLabel"),
			Comment("tag of clusters which are to be EM/track tagged")
		};

		fhicl::Atom<double> Threshold {
			Name("Threshold"),
			Comment("tag cluster as EM-like if net output accumulated over hits > threshold")
		};

		fhicl::Sequence<int> Views {
			Name("Views"),
			Comment("tag clusters in selected views only, or in all views if empty list")
		};
	};
	using Parameters = art::EDProducer::Table<Config>;
	explicit EmTrackClusterId(Parameters const & p);

	EmTrackClusterId(EmTrackClusterId const &) = delete;
	EmTrackClusterId(EmTrackClusterId &&) = delete;
	EmTrackClusterId & operator = (EmTrackClusterId const &) = delete;
	EmTrackClusterId & operator = (EmTrackClusterId &&) = delete;

	void produce(art::Event & e) override;

private:
	bool isViewSelected(int view) const;

	PointIdAlg fPointIdAlg;

	art::InputTag fWireProducerLabel;
	art::InputTag fHitModuleLabel;
	art::InputTag fClusterModuleLabel;

	double fThreshold;

	std::vector< int > fViews;
};
// ------------------------------------------------------

EmTrackClusterId::EmTrackClusterId(EmTrackClusterId::Parameters const& config) :
	fPointIdAlg(config().PointIdAlg()),
	fWireProducerLabel(config().WireLabel()),
	fHitModuleLabel(config().HitModuleLabel()),
	fClusterModuleLabel(config().ClusterModuleLabel()),
	fThreshold(config().Threshold()),
	fViews(config().Views())
{
	produces< std::vector<recob::Cluster> >();
	produces< art::Assns<recob::Cluster, recob::Hit> >();
}
// ------------------------------------------------------

void EmTrackClusterId::produce(art::Event & evt)
{
	auto clusters = std::make_unique< std::vector< recob::Cluster > >();
	auto clu2hit = std::make_unique< art::Assns< recob::Cluster, recob::Hit > >();

	auto wireHandle = evt.getValidHandle< std::vector<recob::Wire> >(fWireProducerLabel);
	auto cluListHandle = evt.getValidHandle< std::vector<recob::Cluster> >(fClusterModuleLabel);

	mf::LogVerbatim("EmTrackClusterId") << "next event: " << evt.run() << " / " << evt.id().event();

	unsigned int cidx = 0;
	const unsigned int emTag = 0x10000;

	// ******************* classify clusters ********************
	art::FindManyP< recob::Hit > hitsFromClusters(cluListHandle, evt, fClusterModuleLabel);
	for (size_t i = 0; i < hitsFromClusters.size(); ++i)
	{
		auto v = hitsFromClusters.at(i);
		if (v.empty()) continue;

		int view = v.front()->View();
		int tpc = v.front()->WireID().TPC;
		int cryo = v.front()->WireID().Cryostat;

		if (!isViewSelected(view)) continue;

		// better if clusters are stored as sorted by tpc/view, otherwise this can be costful call:
		fPointIdAlg.setWireDriftData(*wireHandle, view, tpc, cryo);

		// anyway this is where we spend most time:
		auto vout = fPointIdAlg.predictIdVector(v);
		float pvalue = vout[0] / (vout[0] + vout[1]);

		mf::LogVerbatim("EmTrackClusterId") << "cluster in tpc:" << tpc << " view:" << view
			<< " size:" << v.size() << " p:" << pvalue;

		if (pvalue > fThreshold) continue; // identified as track-like, for the moment we save only em-like parts

		clusters->emplace_back(
			recob::Cluster(0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
				v.size(), 0.0F, 0.0F, cidx + emTag, (geo::View_t)v.front()->View(), v.front()->WireID().planeID()));
		util::CreateAssn(*this, evt, *clusters, v, *clu2hit);
		cidx++;
	}
	// **********************************************************

	// ****************** classify single hits ******************
	art::Handle< std::vector<recob::Hit> > hitListHandle;
	if ((fHitModuleLabel != "") && evt.getByLabel(fHitModuleLabel, hitListHandle))
	{
		std::vector< art::Ptr<recob::Hit> > allhitlist;
		art::fill_ptr_vector(allhitlist, hitListHandle);
		for (auto const & hi : allhitlist)
		{
			bool unclustered = true;
			for (size_t k = 0; k < hitsFromClusters.size(); ++k)
			{
				auto v = hitsFromClusters.at(k);
				for (auto const & hj : v)
				{
					if (hi.key() == hj.key()) { unclustered = false; break; }
				}
				if (!unclustered) break;
			}

			if (unclustered)
			{
				int view = hi->View();
				int tpc = hi->WireID().TPC;
				int cryo = hi->WireID().Cryostat;

				if (!isViewSelected(view)) continue;

				// better if hits are stored as sorted by tpc/view, otherwise this can be costful call:
				fPointIdAlg.setWireDriftData(*wireHandle, view, tpc, cryo);

				// anyway this is where we spend most time:
				auto vout = fPointIdAlg.predictIdVector(hi->WireID().Wire, hi->PeakTime());
				float pvalue = vout[0] / (vout[0] + vout[1]);

				mf::LogVerbatim("EmTrackClusterId") << "hit in tpc:" << tpc << " view:" << view
					<< " wire:" << hi->WireID().Wire << " drift:" << hi->PeakTime() << " p:" << pvalue;

				if (pvalue > fThreshold) continue; // identified as track-like, for the moment we save only em-like parts

				art::PtrVector< recob::Hit > cluster_hits;
				cluster_hits.push_back(hi);
				clusters->emplace_back(
					recob::Cluster(0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F,
					1, 0.0F, 0.0F, cidx + emTag, (geo::View_t)hi->View(), hi->WireID().planeID()));
				util::CreateAssn(*this, evt, *clusters, cluster_hits, *clu2hit);
				cidx++;
			}
		}
	}
	// **********************************************************

	evt.put(std::move(clusters));
	evt.put(std::move(clu2hit));
}
// ------------------------------------------------------

bool EmTrackClusterId::isViewSelected(int view) const
{
	if (fViews.empty()) return true;
	else
	{
		bool selected = false;
		for (auto k : fViews) if (k == view) { selected = true; break; }
		return selected;
	}
}
// ------------------------------------------------------

DEFINE_ART_MODULE(EmTrackClusterId)

}

