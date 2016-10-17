////////////////////////////////////////////////////////////////////////
// Class: ShowerEnergyAlg
// File:  ShowerEnergyAlg.h
// Author: Mike Wallbank (m.wallbank@sheffield.ac.uk), November 2015
//
// Shower energy finding class
////////////////////////////////////////////////////////////////////////

#ifndef ShowerEnergyAlg_hxx
#define ShowerEnergyAlg_hxx

// Framework
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

// larsoft
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

// ROOT
#include "TMath.h"

namespace shower {
  class ShowerEnergyAlg;
}

class shower::ShowerEnergyAlg {
 public:

  ShowerEnergyAlg(fhicl::ParameterSet const& pset);

  /// Finds the total energy deposited by the shower in this view
  /// Performs T0 correction is T0 is supplied (in ns)
  double ShowerEnergy(std::vector<art::Ptr<recob::Hit> > const& hits, int plane, double t0 = 0);

 private:

  double fUGradient, fUIntercept;
  double fVGradient, fVIntercept;
  double fZGradient, fZIntercept;

  detinfo::DetectorProperties const* detprop = nullptr;

};

#endif
