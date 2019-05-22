#ifndef SBNShowerAlg_hxx
#define SBNShowerAlg_hxx

//Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"

//C++ Includes

#include <iostream>
#include <vector>

//Root Includes

#include "TVector3.h"
#include "TMath.h"
#include "TPrincipal.h"
#include "TVector.h"


namespace shower {
  class SBNShowerAlg;
}

class shower::SBNShowerAlg {
 public:
  SBNShowerAlg(const fhicl::ParameterSet& pset);
  void OrderShowerHits(int test=0);

  void OrderShowerSpacePoints(std::vector<art::Ptr<recob::SpacePoint> >& showerspcs,
			      TVector3& vertex, TVector3& direction);


  TVector3 ShowerCentre(std::vector<art::Ptr<recob::SpacePoint> >& showersps,
			art::FindManyP<recob::Hit>& fmh, float& totalCharge);

  
  TVector3 ShowerCentre(std::vector<art::Ptr<recob::SpacePoint> >& showerspcs,
			art::FindManyP<recob::Hit>& fmh);   

  TVector3 SpacePointPosition(const art::Ptr<recob::SpacePoint>& sp);

  double SpacePointCharge(art::Ptr<recob::SpacePoint> sp, art::FindManyP<recob::Hit>& fmh);



 private:

  bool fUseCollectionOnly;

};

#endif
