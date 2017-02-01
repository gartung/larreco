////////////////////////////////////////////////////////////////////
// Implementation of the EMShower algorithm
//
// Forms EM showers from clusters and associated tracks.
// Also provides methods for finding the vertex and further
// properties of the shower.
//
// Mike Wallbank (m.wallbank@sheffield.ac.uk), September 2015
////////////////////////////////////////////////////////////////////

#ifndef EMShowerAlg_hxx
#define EMShowerAlg_hxx

// framework
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// larsoft
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/ArtDataHelper/TrackUtils.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"
#include "larreco/RecoAlg/ProjectionMatchingAlg.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "larreco/RecoAlg/ShowerEnergyAlg.h"

// C++
#include <iostream>
#include <map>
#include <iterator>

// ROOT
#include "TVector2.h"
#include "TMath.h"

//temp
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TString.h"
#include "TF1.h"
#include "larsim/MCCheater/BackTracker.h"
#include "TH1I.h"
#include "TText.h"
#include "TFile.h"
#include "TPrincipal.h"
#include "TProfile.h"

namespace shower {
  class EMShowerAlg;
  class HitPosition;
}

class shower::EMShowerAlg {
public:

  EMShowerAlg(const fhicl::ParameterSet& pset);

  /// Map associated tracks and clusters together given their associated hits
  void AssociateClustersAndTracks(const std::vector<art::Ptr<recob::Cluster> >& clusters,
				  const art::FindManyP<recob::Hit>& fmh,
				  const std::unique_ptr<art::FindManyP<recob::Track> >& fmt,
				  const std::vector<art::Ptr<recob::Vertex> >& vertices,
				  std::map<int,std::vector<int> >& clusterToTracks,
				  std::map<int,std::vector<int> >& trackToClusters);

  /// Map associated tracks and clusters together given their associated hits, whilst ignoring certain clusters
  void AssociateClustersAndTracks(const std::vector<art::Ptr<recob::Cluster> >& clusters,
				  const art::FindManyP<recob::Hit>& fmh,
				  const std::unique_ptr<art::FindManyP<recob::Track> >& fmt,
				  const std::vector<art::Ptr<recob::Vertex> >& vertices,
				  const std::vector<int>& clustersToIgnore,
				  std::map<int,std::vector<int> >& clusterToTracks,
				  std::map<int,std::vector<int> >& trackToClusters);

  /// Takes the initial showers found and tries to resolve issues where one bad view ruins the event
  std::vector<int> CheckShowerPlanes(const std::vector<std::vector<int> >& initialShowers,
  				     const std::vector<art::Ptr<recob::Cluster> >& clusters,
  				     const art::FindManyP<recob::Hit>& fmh);

  /// Constructs a recob::Track from sets of hits in two views. Intended to be used to construct the initial first part of a shower.
  /// All PMA methods taken from the pma tracking algorithm (R. Sulej and D. Stefan).
  std::unique_ptr<recob::Track> ConstructTrack(const std::vector<art::Ptr<recob::Hit> >& track1,
					       const std::vector<art::Ptr<recob::Hit> >& track2);

  /// Finds the initial track-like part of the shower and the hits in all views associated with it
  void FindInitialTrack(const std::map<int,std::vector<art::Ptr<recob::Hit> > >& hits,//art::PtrVector<recob::Hit> const& hits,
			std::unique_ptr<recob::Track>& initialTrack,
			std::map<int,std::vector<art::Ptr<recob::Hit> > >& initialTrackHits, int plane);

  /// Makes showers given a map between tracks and all clusters associated with them
  std::vector<std::vector<int> > FindShowers(const std::map<int,std::vector<int> >& trackToClusters);

  /// Makes a recob::Shower object given the hits in the shower and the initial track-like part
  recob::Shower MakeShower(const art::PtrVector<recob::Hit>& hits,
			   const std::unique_ptr<recob::Track>& initialTrack,
			   const std::map<int,std::vector<art::Ptr<recob::Hit> > >& initialTrackHits,
			   const std::unique_ptr<TVector3>& defaultDirection);

  /// <Tingjun to document>
  recob::Shower MakeShower(const art::PtrVector<recob::Hit>& hits,
			   const art::Ptr<recob::Vertex>& vertex,
			   int & iok);

  /// Makes space points from the shower hits in each plane
  std::vector<recob::SpacePoint> MakeSpacePoints(std::map<int,std::vector<art::Ptr<recob::Hit> > > hits, std::vector<std::vector<art::Ptr<recob::Hit> > >& hitAssns);

  /// Takes the hits associated with a shower and orders them so they follow the direction of the shower
  //std::map<int,std::vector<art::Ptr<recob::Hit> > > OrderShowerHits(art::PtrVector<recob::Hit> const& shower, int plane);
  std::map<int,std::vector<art::Ptr<recob::Hit> > > OrderShowerHits(const art::PtrVector<recob::Hit>& shower, const std::vector<art::Ptr<recob::Vertex> >& vertices, int plane);

  /// Returns the 3D shower centre
  /// If space point-hit assocations are valid, the centre is charge weighted
  TVector3 ShowerCentre(const std::vector<recob::SpacePoint>& showerSpacePoints, const std::vector<std::vector<art::Ptr<recob::Hit> > >& hitAssns);

  /// <Tingjun to document>
  void FindInitialTrackHits(const std::vector<art::Ptr<recob::Hit> >& showerHits,
			    const art::Ptr<recob::Vertex>& vertex,
			    std::vector<art::Ptr<recob::Hit> >& trackHits);
  
  /// <Tingjun to document>
  Int_t WeightedFit(const Int_t n, const Double_t *x, const Double_t *y, const Double_t *w,  Double_t *parm);

  /// <Tingjun to document>
  bool isCleanShower(const std::vector<art::Ptr<recob::Hit> >& hits);

  int fDebug;

private:

  /// Checks the hits across the views in a given shower to determine if there is one in the incorrect TPC
  void CheckIsolatedHits(std::map<int,std::vector<art::Ptr<recob::Hit> > >& showerHitsMap);

  /// Takes the shower hits in all views and ensure the ordering is consistent
  /// Returns bool, indicating whether or not everything makes sense!
  bool CheckShowerHits(const std::map<int,std::vector<art::Ptr<recob::Hit> > >& showerHitsMap);

  /// Takes the correctly oriented hit map and attempts to clean things up
  /// Fixes issues such as where tracks near the vertex are accidentally reconstructed as part of the shower
  /// and is evident in the positioning of the vertex in each plane
  std::map<int,std::vector<art::Ptr<recob::Hit> > > CleanShowerHits(const std::map<int,std::vector<art::Ptr<recob::Hit> > >& showerHits);

  /// Constructs a 3D point (in [cm]) to represent the hits given in two views
  TVector3 Construct3DPoint(const art::Ptr<recob::Hit>& hit1, const art::Ptr<recob::Hit>& hit2);

  /// Constructs a 3D point (in [cm]) given two sets of wire IDs and x coordinates
  TVector3 Construct3DPoint(const std::pair<geo::WireID,double>& hit1, const std::pair<geo::WireID,double>& hit2);

  /// Finds dE/dx for the track given a set of hits
  double FinddEdx(const std::vector<art::Ptr<recob::Hit> >& trackHits, const std::unique_ptr<recob::Track>& track);

  /// Orders hits along the best fit line through the charge-weighted centre of the hits.
  /// Orders along the line perpendicular to the least squares line if perpendicular is set to true.
  std::vector<art::Ptr<recob::Hit> > FindOrderOfHits(const std::vector<art::Ptr<recob::Hit> >& hits, bool perpendicular = false);

  /// Takes a map of the shower hits on each plane (ordered from what has been decided to be the start)
  /// Returns a map of the initial track-like part of the shower on each plane
  std::map<int,std::vector<art::Ptr<recob::Hit> > > FindShowerStart(const std::map<int,std::vector<art::Ptr<recob::Hit> > >& orderedShowerMap);

  /// Takes all the shower hits, ready ordered, and returns information to help with the orientation of the shower in each view
  /// Returns map of most likely permutations of reorientation
  /// Starts at 0,0,0 (i.e. don't need to reorient any plane) and ends with 1,1,1 (i.e. every plane needs reorienting)
  /// Every permutation inbetween represent increasing less likely changes to satisfy the correct orientation criteria
  std::map<int,std::map<int,bool> > GetPlanePermutations(const std::map<int,std::vector<art::Ptr<recob::Hit> > >& showerHitsMap);

  /// Find the global wire position
  double GlobalWire(const geo::WireID& wireID);

  /// Return the coordinates of this hit in global wire/tick space
  TVector2 HitCoordinates(const art::Ptr<recob::Hit>& hit);

  /// Return the coordinates of this hit in units of cm
  TVector2 HitPosition(const art::Ptr<recob::Hit>& hit);

  /// Return the coordinates of this hit in units of cm
  TVector2 HitPosition(const TVector2& pos, geo::PlaneID planeID);

  /// Return the coordiantes of this hit in units geo::WireID and x coordinate
  std::pair<geo::WireID,double> HitWireX(const art::Ptr<recob::Hit>& hit);

  /// Takes initial track hits from multiple views and forms a track object which best represents the start of the shower
  std::unique_ptr<recob::Track> MakeInitialTrack(const std::map<int,std::vector<art::Ptr<recob::Hit> > >& initialHitsMap,
						 const std::map<int,std::vector<art::Ptr<recob::Hit> > >& showerHitsMap);

  /// Takes the hits associated with a shower and orders then so they follow the direction of the shower
  void OrderShowerHits(const std::vector<art::Ptr<recob::Hit> >& shower,
		       std::vector<art::Ptr<recob::Hit> >& orderedShower,
		       const art::Ptr<recob::Vertex>& vertex);

  /// Projects a 3D point (units [cm]) onto a 2D plane
  /// Returns 2D point (units [cm])
  TVector2 Project3DPointOntoPlane(const TVector3& point, int plane, int cryostat = 0);

  /// Projects a 3D point (units [cm]) onto a 2D wire
  /// Returns a pair containing the wire ID and the x coordinate of the point
  std::pair<geo::WireID, double> Project3DPointOntoWire(const TVector3& point, int plane, int cryostat = 0);

  /// Determines the 'relative wire width', i.e. how spread a shower is across wires of each plane relative to the others
  /// If a shower travels along the wire directions in a particular view, it will have a smaller wire width in that view
  /// Returns a map relating these widths to each plane
  std::map<double,int> RelativeWireWidth(const std::map<int,std::vector<art::Ptr<recob::Hit> > >& showerHitsMap);

  /// Returns the charge-weighted shower centre
  TVector2 ShowerCentre(const std::vector<art::Ptr<recob::Hit> >& showerHits);

  /// Returns a rough charge-weighted shower 'direction' given the hits in the shower
  TVector2 ShowerDirection(const std::vector<art::Ptr<recob::Hit> >& showerHits);

  /// Returns the RMS of the hits from the central shower 'axis' along the length of the shower
  double ShowerHitRMS(const std::vector<art::Ptr<recob::Hit> >& showerHits);

  /// Returns the gradient of the RMS vs shower segment graph
  double ShowerHitRMSGradient(const std::vector<art::Ptr<recob::Hit> >& showerHits, TVector2 trueStart = TVector2(0,0));

  /// Returns the least square/d.o.f. for a linear fit through the hits
  double TrackLeastSquare(const std::vector<art::Ptr<recob::Hit> >& hits);

  /// Returns the plane which is determined to be the least likely to be correct
  int WorstPlane(const std::map<int,std::vector<art::Ptr<recob::Hit> > >& showerHitsMap);


  // Parameters
  double fMinTrackLength;
  int fMinTrackShowerHits;
  int fNumShowerSegments;
  double fdEdxTrackLength;
  double fSpacePointSize;
  std::string fBestPlaneMetric;

  // Parameters to fit wire vs time
  unsigned int         fNfitpass;
  std::vector<unsigned int>     fNfithits;
  std::vector<double>  fToler;

  // Services used by this class
  art::ServiceHandle<geo::Geometry> fGeom;
  detinfo::DetectorProperties const* fDetProp;
  art::ServiceHandle<art::TFileService> tfs;

  // Algs used by this class
  shower::ShowerEnergyAlg fShowerEnergyAlg;
  calo::CalorimetryAlg fCalorimetryAlg;
  pma::ProjectionMatchingAlg fProjectionMatchingAlg;

  std::string fDetector;


  // tmp
  void MakePicture();
  bool fMakeGradientPlot, fMakeRMSGradientPlot;

};

class shower::HitPosition {
 public:

  TVector2 WireTick;
  TVector2 Cm;

  // default constructor
  HitPosition();
  // contructors
  HitPosition(TVector2 wiretick, TVector2 cm): HitPosition()
  {
    WireTick = wiretick;
    Cm = cm;
  }
  HitPosition(TVector2 wiretick, geo::PlaneID planeID): HitPosition() {
    WireTick = wiretick;
    Cm = ConvertWireTickToCm(wiretick, planeID);
  }

  TVector2 ConvertWireTickToCm(TVector2 wiretick, geo::PlaneID planeID) {
    return TVector2(wiretick.X() * fGeom->WirePitch(planeID),
		    fDetProp->ConvertTicksToX(wiretick.Y(), planeID));
  }

 private:

  geo::GeometryCore const* fGeom = nullptr;
  detinfo::DetectorProperties const* fDetProp = nullptr;

};

#endif
