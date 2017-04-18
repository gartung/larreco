////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgCosmicTagger
// Author:      ..., R.Sulej (..., robert.sulej@cern.ch) March 2017
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/PMAlgCosmicTagger.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

void pma::PMAlgCosmicTagger::tag(pma::TrkCandidateColl& tracks)
{
    mf::LogInfo("pma::PMAlgCosmicTagger") << "Passed " << tracks.size() << " tracks for tagging cosmics.";

    size_t n = 0;

    if (fTagOutOfDriftTracks) n += outOfDriftWindow(tracks);
    if (fTagTopDownTracks) n += topDownCrossing(tracks);
		if (fTagNonBeamT0Tracks) n += nonBeamT0Tag(tracks);

    mf::LogInfo("pma::PMAlgCosmicTagger") << "...tagged " << n << " cosmic-like tracks.";
}


size_t pma::PMAlgCosmicTagger::outOfDriftWindow(pma::TrkCandidateColl& tracks)
{
    mf::LogInfo("pma::PMAlgCosmicTagger") << "   - tag tracks out of 1 drift window;";
    size_t n = 0;

    auto const* geom = lar::providerFrom<geo::Geometry>();

    for (auto & t : tracks.tracks())
    {
        pma::Track3D & trk = *(t.Track());

        double min, max, p;
        bool node_out_of_drift = false;
        for (size_t nidx = 0; nidx < trk.Nodes().size(); ++nidx)
        {
            auto const & node = *(trk.Nodes()[nidx]);
            auto const & tpcGeo = geom->TPC(node.TPC(), node.Cryo());
            // DetectDriftDirection returns a short int, but switch requires an int
            int driftDir = abs(tpcGeo.DetectDriftDirection());
            p = node.Point3D()[driftDir-1];
            switch (driftDir)
            {
                case 1: min = tpcGeo.MinX(); max = tpcGeo.MaxX(); break;
                case 2: min = tpcGeo.MinY(); max = tpcGeo.MaxY(); break;
                case 3: min = tpcGeo.MinZ(); max = tpcGeo.MaxZ(); break;
                default: throw cet::exception("PMAlgCosmicTagger") << "Drift direction unknown: " << driftDir << std::endl;
            }
            if ((p < min - fOutOfDriftMargin) || (p > max + fOutOfDriftMargin)) { node_out_of_drift = true; break; }
        }

        if (node_out_of_drift) { trk.SetTagFlag(pma::Track3D::kCosmic); ++n; }
    }

    return n;
}

size_t pma::PMAlgCosmicTagger::topDownCrossing(pma::TrkCandidateColl& tracks)
{
    mf::LogInfo("pma::PMAlgCosmicTagger") << "   - tag tracks crossing full Y range (if it is not drift);";
    size_t n = 0;


    return n;
}

// Leigh: Make use of the fact that our cathode and anode crossing tracks have a reconstructed T0.
// Check to see if this time is consistent with the beam
size_t pma::PMAlgCosmicTagger::nonBeamT0Tag(pma::TrkCandidateColl &tracks){

	size_t n = 0;

	auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

	// Search through all of the tracks
	for(auto & t : tracks){

		// Non zero T0 means we reconstructed it
		if(t->GetT0() != 0.0){

			if(fabs(t->GetT0() - detprop->GetTriggerOffset()) < fNonBeamT0Margin){
				++n;
				t.SetTagFlag(pma::Track3D::kCosmic);
			}

		}

	}

	mf::LogInfo("pma::PMAlgCosmicTagger") << " - Tagged " << n << " non-beam T0 tracks.";
	return n;

}


