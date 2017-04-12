#ifndef GAUSHITFINDER_H
#define GAUSHITFINDER_H

////////////////////////////////////////////////////////////////////////
//
// GaussHitFinder class
//
// jaasaadi@syr.edu
//
//  This algorithm is designed to find hits on wires after deconvolution.
// -----------------------------------
// This algorithm is based on the FFTHitFinder written by Brian Page,
// Michigan State University, for the ArgoNeuT experiment.
//
//
// The algorithm walks along the wire and looks for pulses above threshold
// The algorithm then attempts to fit n-gaussians to these pulses where n
// is set by the number of peaks found in the pulse
// If the Chi2/NDF returned is "bad" it attempts to fit n+1 gaussians to
// the pulse. If this is a better fit it then uses the parameters of the
// Gaussian fit to characterize the "hit" object
//
// To use this simply include the following in your producers:
// gaushit:     @local::microboone_gaushitfinder
// gaushit: @local::argoneut_gaushitfinder
////////////////////////////////////////////////////////////////////////


// C/C++ standard library
#include <algorithm> // std::accumulate()
#include <vector>
#include <string>
#include <memory> // std::unique_ptr()
#include <utility> // std::move()

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"


// LArSoft Includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/HitFinder/GausHitFinderAlg.h"


// ROOT Includes
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TDecompSVD.h"
#include "TMath.h"
#include "TF1.h"
#include "TTree.h"
#include "TStopwatch.h"

namespace hit {
class GausHitFinderOrganized : public art::EDProducer {

public:

    explicit GausHitFinderOrganized(fhicl::ParameterSet const& pset);

    void produce(art::Event& evt) override;
    void beginJob() override;
    void endJob() override;
    void reconfigure(fhicl::ParameterSet const& p) override;

private:
    std::string         fCalDataModuleLabel;
    GausHitFinderAlg    fHitFinderAlg;

protected:


}; // class GausHitFinderOrganized


//-------------------------------------------------
//-------------------------------------------------
GausHitFinderOrganized::GausHitFinderOrganized(fhicl::ParameterSet const& pset)
    : fHitFinderAlg(pset)
{
    this->reconfigure(pset);

    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionCreator::declare_products(*this);

} // GausHitFinderOrganized::GausHitFinderOrganized()



//-------------------------------------------------
//-------------------------------------------------
void GausHitFinderOrganized::reconfigure(fhicl::ParameterSet const& p)
{
    // Implementation of optional member function here.
    fCalDataModuleLabel = p.get< std::string  >("CalDataModuleLabel");


}

//-------------------------------------------------
//-------------------------------------------------
void GausHitFinderOrganized::beginJob()
{

    fHitFinderAlg.beginJob();

}

//-------------------------------------------------
//-------------------------------------------------
void GausHitFinderOrganized::endJob()
{

}

//  This algorithm uses the fact that deconvolved signals are very smooth
//  and looks for hits as areas between local minima that have signal above
//  threshold.
//-------------------------------------------------
void GausHitFinderOrganized::produce(art::Event& evt)
{
    //==================================================================================================

    TH1::AddDirectory(kFALSE);

    // Instantiate and Reset a stop watch
    //TStopwatch StopWatch;
    //StopWatch.Reset();

    // ################################
    // ### Calling Geometry service ###
    // ################################
    art::ServiceHandle<geo::Geometry> geom;


    // ###############################################
    // ### Making a ptr vector to put on the event ###
    // ###############################################
    // this contains the hit collection
    // and its associations to wires and raw digits
    recob::HitCollectionCreator hcol(*this, evt);

    // ##########################################
    // ### Reading in the Wire List object(s) ###
    // ##########################################
    art::Handle< std::vector<recob::Wire> > wireVecHandle;
    evt.getByLabel(fCalDataModuleLabel, wireVecHandle);

    // #################################################################
    // ### Reading in the RawDigit associated with these wires, too  ###
    // #################################################################
    art::FindOneP<raw::RawDigit> RawDigits
    (wireVecHandle, evt, fCalDataModuleLabel);



    //##############################
    //### Looping over the wires ###
    //##############################
    for (size_t wireIter = 0; wireIter < wireVecHandle->size(); wireIter++)
    {
        // ####################################
        // ### Getting this particular wire ###
        // ####################################
        art::Ptr<recob::Wire>   wire(wireVecHandle, wireIter);
        art::Ptr<raw::RawDigit> rawdigits = RawDigits.at(wireIter);

        std::vector<recob::Hit> _hit_vector;
        fHitFinderAlg.RunOnWire(wire, _hit_vector, *geom);
        for (auto & hit : _hit_vector) {
            hcol.emplace_back(std::move(hit), wire, rawdigits);
        }


    }//<---End looping over all the wires


    //==================================================================================================
    // End of the event

    // move the hit collection and the associations into the event
    hcol.put_into(evt);

} // End of produce()





DEFINE_ART_MODULE(GausHitFinderOrganized)

} // end of hit namespace
#endif // GAUSHITFINDER_H
