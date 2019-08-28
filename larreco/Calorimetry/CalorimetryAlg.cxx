////////////////////////////////////////////////////////////////////////
//  \file CalorimetryAlg.cxx
//
//  \brief Functions to calculate dE/dx. Based on code in Calorimetry.cxx
//
// andrzej.szelc@yale.edu
//
////////////////////////////////////////////////////////////////////////

// LArSoft includes
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larevt/CalibrationDBI/Interface/ElectronLifetimeProvider.h"
#include "larevt/CalibrationDBI/Interface/ElectronLifetimeService.h"

namespace calo {

  //--------------------------------------------------------------------
  CalorimetryAlg::CalorimetryAlg(const Config& config)
    : detprop{art::ServiceHandle<detinfo::DetectorPropertiesService const>()
                ->provider()}
    , fCalAmpConstants{config.CalAmpConstants()}
    , fCalAreaConstants{config.CalAreaConstants()}
    , fUseModBox{config.CaloUseModBox()}
    , fLifeTimeForm{config.CaloLifeTimeForm()}
    , fDoLifeTimeCorrection{config.CaloDoLifeTimeCorrection()}
  {
    if (fLifeTimeForm != 0 || fLifeTimeForm != 1) {
      throw cet::exception("CalorimetryAlg")
        << "Unknow CaloLifeTimeForm " << fLifeTimeForm << '\n'
        << "Must select either '0' for exponential or '1' for exponential + constant.\n";
    }
  }

  //------------------------------------------------------------------------------------//
  // Functions to calculate the dEdX based on the AMPLITUDE of the pulse
  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AMP(recob::Hit const& hit,
                           double const pitch,
                           double const T0) const
  {
    return dEdx_AMP(
      hit.PeakAmplitude() / pitch, hit.PeakTime(), hit.WireID().Plane, T0);
  }

  ///\todo The plane argument should really be for a view instead
  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AMP(double const dQ,
                           double const time,
                           double const pitch,
                           unsigned int const plane,
                           double const T0) const
  {
    double const dQdx = dQ / pitch; // in ADC/cm
    return dEdx_AMP(dQdx, time, plane, T0);
  }

  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AMP(double const dQdx,
                           double const time,
                           unsigned int const plane,
                           double const T0) const
  {
    double const fADCtoEl = fCalAmpConstants[plane];
    double const dQdx_e = dQdx / fADCtoEl; // Conversion from ADC/cm to e/cm
    return dEdx_from_dQdx_e(dQdx_e, time, T0);
  }

  //------------------------------------------------------------------------------------//
  // Functions to calculate the dEdX based on the AREA of the pulse
  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AREA(recob::Hit const& hit,
                            double const pitch,
                            double const T0) const
  {
    return dEdx_AREA(
      hit.Integral() / pitch, hit.PeakTime(), hit.WireID().Plane, T0);
  }

  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AREA(double const dQ,
                            double const time,
                            double const pitch,
                            unsigned int const plane,
                            double const T0) const
  {
    double const dQdx = dQ / pitch; // in ADC/cm
    return dEdx_AREA(dQdx, time, plane, T0);
  }

  // ----------------------------------------------------------------------------------//
  double
  CalorimetryAlg::dEdx_AREA(double const dQdx,
                            double const time,
                            unsigned int const plane,
                            double const T0) const
  {
    double const fADCtoEl = fCalAreaConstants[plane];
    double const dQdx_e = dQdx / fADCtoEl; // Conversion from ADC/cm to e/cm
    return dEdx_from_dQdx_e(dQdx_e, time, T0);
  }

  // ----------------- apply Lifetime and recombination correction.
  // -----------------//
  double
  CalorimetryAlg::dEdx_from_dQdx_e(double dQdx_e,
                                   double const time,
                                   double const T0) const
  {
    if (fDoLifeTimeCorrection) {
      dQdx_e *= LifetimeCorrection(time, T0); // (dQdx_e in e/cm)
    }

    if (fUseModBox) {
      return detprop->ModBoxCorrection(dQdx_e);
    }

    return detprop->BirksCorrection(dQdx_e);
  }

  //------------------------------------------------------------------------------------//
  // for the time being copying from Calorimetry.cxx - should be decided where
  // to keep it.
  // ----------------------------------------------------------------------------------//
  double
  calo::CalorimetryAlg::LifetimeCorrection(double const time, double const T0) const
  {
    float const t = time - detprop->TriggerOffset();
    double const timetick = detprop->SamplingRate() * 1.e-3; // time sample in microsec
    double const adjusted_time = t * timetick - T0 * 1e-3; //  (in microsec)

    assert(fLifeTimeForm < 2);
    if (fLifeTimeForm == 0) {
      // Exponential form
      double const tau = detprop->ElectronLifetime();
      return exp(adjusted_time / tau);
    }

    // Exponential+constant form
    auto const& elifetime_provider =
      art::ServiceHandle<lariov::ElectronLifetimeService const>()
      ->GetProvider();
    return elifetime_provider.Lifetime(adjusted_time);
  }

} // namespace
