////////////////////////////////////////////////////////////////////////
//
// MagDriftAna class
//
// dmckee@phys.ksu.edu
//
//
////////////////////////////////////////////////////////////////////////
// C++ std library includes
#include <string>
#include <algorithm>
// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// Root Includes
#include "TLine.h"
#include "TH2.h"
// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "nug4/MagneticField/MagneticField.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/RecoBase/Hit.h"

///Detector simulation of raw signals on wires
namespace hit {

  /// Base class for creation of raw signals on wires.
  class MagDriftAna : public art::EDAnalyzer {

  public:

    explicit MagDriftAna(fhicl::ParameterSet const& pset);

  private:

    /// read/write access to event
    void analyze (const art::Event& evt);
    void endJob();

    // intilize the histograms
    //
    // Can't be done in Begin job because I want to use LArProperties
    // which used the database, so I test and run on each
    // event. Wasteful and silly, but at least it *works*.
    void ensureHists();

    std::string            fFFTHitFinderModuleLabel;
    std::string            fTrackFinderModuleLabel;
    std::string            fLArG4ModuleLabel;

    // Flag for initialization done, because we set up histograms the
    // first time through beginRun() so that we can use the
    // database...
    bool initDone;

    // Drift properties
    double fDirCosY;
    double fDirCosZ;

    TH1D * fChargeXpos;  // << position of the MC Truth charge deposition
    TH1D * fChargeYpos;
    TH1D * fChargeZpos;
    TH1D * fHitZpos;     // << Z position of the recorded hit (from the
                         //    z-sensitive wire)

    TH1D * fDriftDeltaZ; // << Difference in MC charge Z and recorded hit Z
    TH1D * fDeltaZoverX; // << Delta Z as a function of drift distance
    TH2D * fDeltaZvsX;

    // Same as above, but only for long drift distances (greater than
    // 4/5 of the detector)
    TH1D * fDriftDeltaZAway; // << Difference in MC charge Z and recorded hit Z
    TH1D * fDeltaZoverXAway; // << Delta Z as a function of drift distance

  }; // class MagdriftAna


  //-------------------------------------------------
  MagDriftAna::MagDriftAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , initDone(false)
    , fDirCosY(0.0)
    , fDirCosZ(0.0)
    , fChargeXpos()
    , fChargeYpos()
    , fChargeZpos()
    , fHitZpos()
    , fDriftDeltaZ()
    , fDeltaZoverX()
    , fDeltaZvsX()
    , fDriftDeltaZAway()
    , fDeltaZoverXAway()
  {
    fFFTHitFinderModuleLabel  = pset.get< std::string >("HitsModuleLabel");
//     fTrackFinderModuleLabel   = pset.get< std::string >("TracksModuleLabel");
    fLArG4ModuleLabel         = pset.get< std::string >("LArGeantModuleLabel");
  }

  //-------------------------------------------------
  void MagDriftAna::ensureHists() {
    if (initDone) return; // Bail if we've already done this.
    initDone = true; // Insure that we bail later on

    // get access to the TFile service
    art::ServiceHandle<art::TFileService const> tfs;
    // Find magetic field related corrections
    const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    art::ServiceHandle<mag::MagneticField const> MagField;
    double Efield         = detprop->Efield();
    double Temperature    = detprop->Temperature();
    double DriftVelocity  = detprop->DriftVelocity(Efield,Temperature)/1000.;

    // MagneticField::FieldAtPoint() returns (0, 0, 0) if there is no
    // field at the requested point, so these direction cosines are
    // 0 if there is no field
    fDirCosY = -DriftVelocity * MagField->FieldAtPoint().z() / Efield;
    fDirCosZ = +DriftVelocity * MagField->FieldAtPoint().y() / Efield;
    MF_LOG_VERBATIM("MagDriftAna")
    << "Drift ratios: "
                << "dY/dX = " << fDirCosY << ", "
    << "dZ/dX = " << fDirCosZ;

    // geometry data.
    art::ServiceHandle<geo::Geometry const> geom;
    // assumes all TPCs are the same
    double width  = 2 * geom->TPC(0).HalfWidth();
    double halfHeight = geom->TPC(0).HalfHeight();
    double length = geom->TPC(0).Length();

    double zScale = std::max(fDirCosZ/2.0,4e-4);

    // Assumes microboone dimensions. Ideally we'd fix this later...
    fChargeXpos  = tfs->make<TH1D>("hChargeXpos",
                                   "MC X charge depositions; X (cm); Events",
                                   101, 0.0, width);
    fChargeYpos  = tfs->make<TH1D>("hChargeYpos",
                                   "MC Y charge depositions; Y (cm); Events",
                                   101, -halfHeight, halfHeight);
    fChargeZpos  = tfs->make<TH1D>("hChargeZpos",
                                   "MC Z charge depositions; Z (cm); Events",
                                   101, 0.0, length);
    fHitZpos     = tfs->make<TH1D>("hHitZpos",
                                   "Z charge collection; Z (cm); Events",
                                   101, 0.0, length);

    fDriftDeltaZ = tfs->make<TH1D>("hDriftDeltaZ",
                                   "Z drift of charge; delta Z (cm); Events",
                                   101, -5*zScale*width, 5*zScale*width);
    fDeltaZoverX = tfs->make<TH1D>("hDeltaZoverX",
                                   "Z drift of charge; delta Z/X; Events",
                                   51, -10*zScale, 10*zScale);
    fDeltaZvsX = tfs->make<TH2D>("hDeltaZvsX",
                                 "delta Z vs X; X (cm); delta Z (cm), Events",
                                 51, 0.0, width,
                                 51, -20*zScale, 20*zScale);

    // Some stats only when the Xdrift is large (more than 3/4)
    fDriftDeltaZAway =
    tfs->make<TH1D>("hDriftDeltaZAway",
                    "Z drift of charge (long drift); delta Z (cm); Events",
                    101, -5*zScale*width, 5*zScale*width);
    fDeltaZoverXAway =
    tfs->make<TH1D>("hDeltaZoverXAway",
                    "Z drift of charge (long drift); delta Z/X; Events",
                    51, -10*zScale, 10*zScale);

    return;

  }

  //-------------------------------------------------
  void MagDriftAna::endJob()
  {
    // Add a line on the deltaZ/X graph to denote the calculated value
    // of the drift ration
    TLine *l = new TLine(fDirCosZ, 0,
                         fDirCosZ, 1.05*fDeltaZoverX->GetMaximum());
    l->SetLineColor(kRed);
    l->SetLineStyle(kDotted);
    fDeltaZoverX->GetListOfFunctions()->Add(l);

    // I know this looks like a memory leak, but each historgram needs
    // it's own copy of the line to prevent double freeing by the
    // framework...
    ////////////
    l = new TLine(fDirCosZ, 0,
                  fDirCosZ, 1.05*fDeltaZoverX->GetMaximum());
    ///////////
    //
    l->SetLineColor(kRed);
    l->SetLineStyle(kDotted);
    fDeltaZoverXAway->GetListOfFunctions()->Add(l);
  }

  //-------------------------------------------------
  void MagDriftAna::analyze(const art::Event& evt)
  {

    if (evt.isRealData()) {
      throw cet::exception("MagDriftAna: ") << "Not for use on Data yet...\n";
    }

    ensureHists();

    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(fFFTHitFinderModuleLabel,hitHandle);

//     art::Handle< std::vector<recob::Track> > trackHandle;
//     evt.getByLabel(fTrackFinderModuleLabel,trackHandle);

    art::ServiceHandle<geo::Geometry const> geom;

    // We're going to want to compare the reconstructed Z with the
    // simulted Z. For that purpose we use the simultion backtracking.
    //

    art::ServiceHandle<cheat::BackTrackerService const> bt_serv;

    //    art::PtrVector<recob::Hit> hits;
    std::vector< art::Ptr<recob::Hit> > hits;
    art::fill_ptr_vector(hits, hitHandle);
//     std::vector< art::Ptr<recob::Track> > tracks;
//     art::fill_ptr_vector(tracks, trackHandle);

    geo::WireID hitWireID;

    //++++++++++
    // Loop over the hits (data) and fill histos
    //++++++++++
    for ( auto itr : hits) {

      hitWireID = itr->WireID();
      // By assumption the drift occurs only in the z-direction, so
      // we can get all the info we need from the z-measug plane.
      if (hitWireID.Plane != (geom->Nplanes()-1)) continue;

      // Charge collected at the wire
      //
      // Exactly once for each recob::Hit
      double w0pos[3] = {0.};
      geom->TPC(hitWireID.TPC).Plane(hitWireID.Plane).Wire(0).GetCenter(w0pos);
      double HitZpos = w0pos[2] + hitWireID.Wire * geom->TPC(hitWireID.TPC).WirePitch();
      double Charge = itr->Integral();
      fHitZpos->Fill(HitZpos,Charge);

      // Charge deposition in the detector
      std::vector<double> xyz = bt_serv->HitToXYZ(itr);
      fChargeXpos->Fill(xyz[0],Charge);
      fChargeYpos->Fill(xyz[1],Charge);
      double ChargeZpos = xyz[2];
      fChargeZpos->Fill(ChargeZpos,Charge);

//      // Delta-Y
//      fDriftDeltaY->Fill(HitYpos-ChargeYpos,Charge);
//      // Delta Y correlation with X
//      fDeltaYoverX->Fill((HitYpos-ChargeYpos)/xyz[0],Charge);
//      fDeltaYvsX->Fill(xyz[0],HitYpos-ChargeZpos,Charge);

      // Delta-Z
      //
      // Compares the collected z position from the wire to the
      // simulated z position
      double DeltaZ = HitZpos-ChargeZpos;
      fDriftDeltaZ->Fill(DeltaZ,Charge);
      // Delta Z correlation with X
      fDeltaZoverX->Fill(DeltaZ/xyz[0],Charge);
      fDeltaZvsX->Fill(xyz[0],DeltaZ,Charge);
      // The X related histograms know the dimensions of the
      // detector, so we use them to set the "away" limit
      if (xyz[0] >  (fChargeYpos->GetXaxis()->GetXmax() * 0.80) ) {
        fDriftDeltaZAway->Fill(DeltaZ,Charge);
        fDeltaZoverXAway->Fill(DeltaZ/xyz[0],Charge);
      }

    } // loop on Hits

    return;
  }//end analyze method

  DEFINE_ART_MODULE(MagDriftAna)

} // end of hit namespace
