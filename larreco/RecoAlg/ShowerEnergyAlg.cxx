////////////////////////////////////////////////////////////////////////
// Class: ShowerEnergyAlg
// File:  ShowerEnergyAlg.cxx
// Author: Mike Wallbank (m.wallbank@sheffield.ac.uk), November 2015
//
// Shower energy finding class
////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/ShowerEnergyAlg.h"

shower::ShowerEnergyAlg::ShowerEnergyAlg(fhicl::ParameterSet const& pset)
  : detprop(lar::providerFrom<detinfo::DetectorPropertiesService>())
{
  fUGradient  = pset.get<double>("UGradient");
  fUIntercept = pset.get<double>("UIntercept");
  fVGradient  = pset.get<double>("VGradient");
  fVIntercept = pset.get<double>("VIntercept");
  fZGradient  = pset.get<double>("ZGradient");
  fZIntercept = pset.get<double>("ZIntercept");
}

double shower::ShowerEnergyAlg::ShowerEnergy(std::vector<art::Ptr<recob::Hit> > const& hits, int plane, double t0) {

  double totalCharge = 0, totalEnergy = 0;
  double time = 0;

  for (art::PtrVector<recob::Hit>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit) {
    if (int((*hit)->WireID().Plane) != plane)
      continue;
    time = (*hit)->PeakTime() * detprop->SamplingRate() * 1e-3 /*us*/ - t0 * 1e-3 /*us*/;
    totalCharge += (*hit)->Integral() * TMath::Exp(time/detprop->ElectronLifetime());
  }

  switch (plane) {
  case 0:
    totalEnergy = (totalCharge * fUGradient) + fUIntercept;
    break;
  case 1:
    totalEnergy = (totalCharge * fVGradient) + fVIntercept;
    break;
  case 2:
    totalEnergy = (totalCharge * fZGradient) + fZIntercept;
    break;
  }

  return totalEnergy;

}
