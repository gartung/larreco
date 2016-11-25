////////////////////////////////////////////////////////////////////////
// Class: TrackShowerSepAlg
// File:  TrackShowerSepAlg.h
// Author: Mike Wallbank (m.wallbank@sheffield.ac.uk), November 2015
//
// Track/shower separation class.
// Provides methods for removing hits associated with track-like
// objects.
// To be run after track reconstruction, before shower reconstruction.
////////////////////////////////////////////////////////////////////////

#ifndef TrackShowerSepAlg_hxx
#define TrackShowerSepAlg_hxx

// Framework
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// larsoft
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"

// ROOT
#include "TVector2.h"
#include "TTree.h"
#include "TMath.h"
#include "TPrincipal.h"
#include "TGraph2D.h"
#include "TPolyLine3D.h"
#include "TCanvas.h"

namespace shower {
  class TrackShowerSepAlg;
  class ReconTrack;
}

class shower::ReconTrack {
 public:

  ReconTrack(int id) {
    fID = id;
    fTrack = false;
    fShower = false;
    fShowerTrack = false;
    fShowerCone = false;

    fNumHits = 0;
    fNumRectangleHits = 0;

    fMinForwardSpacePoint = 1e3;
    fMinBackwardSpacePoint = 1e3;
  }

  // Setters
  void AddPlane(int plane) { fPlanes.push_back(plane); }
  void SetVertex(TVector3 vertex) { fVertex = vertex; }
  void SetVertex2D(std::map<int,TVector2> vertices) { fVertex2D = vertices; }
  void SetVertex2D(int plane, TVector2 vertex) { fVertex2D[plane] = vertex; }
  void SetEnd(TVector3 end) { fEnd = end; }
  void SetEnd2D(std::map<int,TVector2> ends) { fEnd2D = ends; }
  void SetEnd2D(int plane, TVector2 end) { fEnd2D[plane] = end; }
  void SetLength(double length) { fLength = length; }
  void SetVertexDir(TVector3 vertexDir) { fVertexDir = vertexDir; }
  void SetDirection3D(TVector3 direction) { fDirection3D = direction; }
  void SetDirection2D(std::map<int,TVector2> directions) { fDirection2D = directions; }
  void SetDirection2D(int plane, TVector2 direction) { fDirection2D[plane] = direction; }
  void SetHits(std::vector<art::Ptr<recob::Hit> > hits) {
    fNumHits = hits.size();
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt)
      fHits[(*hitIt)->WireID().Plane].push_back(*hitIt);
  }
  void SetSpacePoints(std::vector<art::Ptr<recob::SpacePoint> > spacePoints) { fSpacePoints = spacePoints; }
  void SetLeastSquareNDOF(double ls) { fLeastSquareNDOF = ls; }

  void AddForwardTrack(int track) {
    if (std::find(fForwardConeTracks.begin(), fForwardConeTracks.end(), track) == fForwardConeTracks.end())
      fForwardConeTracks.push_back(track);
  }
  void AddBackwardTrack(int track) {
    if (std::find(fBackwardConeTracks.begin(), fBackwardConeTracks.end(), track) == fBackwardConeTracks.end())
      fBackwardConeTracks.push_back(track);
  }
  void AddShowerTrack(int track) { fShowerTracks.push_back(track); }

  void AddForwardSpacePoint(int spacePoint, double distance) {
    fForwardSpacePoints.push_back(spacePoint);
    if (distance < fMinForwardSpacePoint)
      fMinForwardSpacePoint = distance;
  }
  void AddBackwardSpacePoint(int spacePoint, double distance) {
    fBackwardSpacePoints.push_back(spacePoint);
    if (distance < fMinBackwardSpacePoint)
      fMinBackwardSpacePoint = distance;
  }
  void AddCylinderSpacePoint(int spacePoint) { fCylinderSpacePoints.push_back(spacePoint); }
  void AddSphereSpacePoint(int spacePoint) { fSphereSpacePoints.push_back(spacePoint); }
  void AddIsolationSpacePoint(int spacePoint, double distance) { fIsolationSpacePoints[spacePoint] = distance; }
  void AddRectangleHit(const art::Ptr<recob::Hit>& hit) { fRectangleHits[hit->WireID().Plane].push_back(hit.key()); ++fNumRectangleHits; }

  // Getters
  int ID() const { return fID; }
  TVector3 Vertex() const { return fVertex; }
  TVector3 End() const { return fEnd; }
  double Length() const { return fLength; }
  TVector3 VertexDirection() const { return fVertexDir; }
  TVector3 Direction3D() const { return fDirection3D; }
  TVector3 Centre3D() const {
    TVector3 centre = TVector3(0,0,0);
    for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spIt = fSpacePoints.begin(); spIt != fSpacePoints.end(); ++spIt) {
      const double* xyz = (*spIt)->XYZ();
      centre += TVector3(xyz[0], xyz[1], xyz[2]);
    }
    centre *= 1./(double)fSpacePoints.size();
    return centre;
  }
  TVector2 Vertex2D(int plane) { return fVertex2D[plane]; }
  TVector2 End2D(int plane) { return fEnd2D[plane]; }
  TVector2 Direction2D(int plane) { return fDirection2D[plane]; }
  const std::vector<art::Ptr<recob::Hit> >& Hits(int plane) { return fHits[plane]; }
  const std::map<int,std::vector<art::Ptr<recob::Hit> > >& Hits() const { return fHits; }
  const std::vector<art::Ptr<recob::SpacePoint> >& SpacePoints() const { return fSpacePoints; }

  void FlipTrack() {
    // Vertex and end
    TVector3 tmp3D = fEnd;
    fEnd = fVertex;
    fVertex = tmp3D;
    fDirection3D *= -1;
    for (std::vector<int>::iterator plane = fPlanes.begin(); plane != fPlanes.end(); ++plane) {
      TVector2 tmp2D = fEnd2D[*plane];
      fEnd2D[*plane] = fVertex2D[*plane];
      fVertex2D[*plane] = tmp2D;
      fDirection2D[*plane] *= -1;
    }
    // Cone tracks and space points
    std::vector<int> tmpVec = fForwardSpacePoints;
    fForwardSpacePoints = fBackwardSpacePoints;
    fBackwardSpacePoints = tmpVec;
    tmpVec = fForwardConeTracks;
    fForwardConeTracks = fBackwardConeTracks;
    fBackwardConeTracks = tmpVec;
    double tmpDouble;
    tmpDouble = fMinForwardSpacePoint;
    fMinForwardSpacePoint = fMinBackwardSpacePoint;
    fMinBackwardSpacePoint = tmpDouble;
  }

  void MakeShower() {
    if (fTrack)
      this->MakeShowerTrack();
    else
      this->MakeShowerCone();
  }
  void MakeShowerTrack() {
    fShower = true;
    fShowerTrack = true;
    fShowerCone = false;
    fTrack = false;
  }
  void MakeShowerCone() {
    fShower = true;
    fShowerCone = true;
    fShowerTrack = false;
    fTrack = false;
  }
  void MakeTrack() {
    fTrack = true;
    fShower = false;
    fShowerTrack = false;
    fShowerCone = false;
  }

  bool IsShower() const { return fShower; }
  bool IsShowerTrack() const { return fShowerTrack; }
  bool IsShowerCone() const { return fShowerCone; }
  bool IsTrack() const { return fTrack; }
  bool IsUndetermined() const { return !fTrack and !fShower; }

  int TrackConeSize() const { return (int)fForwardConeTracks.size() - (int)fBackwardConeTracks.size(); }
  bool ShowerTrackCandidate() const { return TrackConeSize() > 5; }
  const std::vector<int>& ShowerTracks() const { return fShowerTracks; }
  const std::vector<int>& ForwardConeTracks() const { return fForwardConeTracks; }

  int ConeSize() const { return (int)fForwardSpacePoints.size() - (int)fBackwardSpacePoints.size(); }
  int ForwardSpacePoints() const { return fForwardSpacePoints.size(); }
  double MinForwardSpacePointDistance() const { return fMinForwardSpacePoint; }
  double MinBackwardSpacePointDistance() const { return fMinBackwardSpacePoint; }
  int NumCylinderSpacePoints() const { return fCylinderSpacePoints.size(); }
  double CylinderSpacePointRatio() const { return (double)fCylinderSpacePoints.size()/(double)fSpacePoints.size(); }
  int NumSphereSpacePoints() const { return fSphereSpacePoints.size(); }
  //double SphereSpacePointDensity() const { return (double)fSphereSpacePoints.size()/((double)fSpacePoints.size()); }
  double SphereSpacePointDensity(double scale) const { return (double)fSphereSpacePoints.size()/(4*TMath::Pi()*TMath::Power((scale*fLength/2.),3)/3.); }
  int NumRectangleHits() const { return fNumRectangleHits; }
  int NumRectangleHits(int plane) { return fRectangleHits[plane].size(); }
  double RectangleHitRatio() const { return (double)TMath::Power(fNumRectangleHits,1)/(double)fNumHits; }
  double RectangleHitRatio(int plane) { return (double)TMath::Power(fRectangleHits[plane].size(),1)/fHits[plane].size(); }
  double IsolationSpacePointDistance() const {
    std::vector<double> distances;
    std::transform(fIsolationSpacePoints.begin(), fIsolationSpacePoints.end(), std::back_inserter(distances), [](const std::pair<int,double>& p){return p.second;});
    return TMath::Mean(distances.begin(), distances.end());
  }
  double LeastSquareNDOF() const { return fLeastSquareNDOF; }

 private:

  int fID;
  std::vector<int> fPlanes;
  TVector3 fVertex;
  std::map<int,TVector2> fVertex2D;
  TVector3 fEnd;
  std::map<int,TVector2> fEnd2D;
  double fLength;
  TVector3 fVertexDir;
  TVector3 fDirection3D;
  std::map<int,TVector2> fDirection2D;
  std::map<int,std::vector<art::Ptr<recob::Hit> > > fHits;
  int fNumHits;
  std::vector<art::Ptr<recob::SpacePoint> > fSpacePoints;

  double fLeastSquareNDOF;

  std::vector<int> fForwardConeTracks;
  std::vector<int> fBackwardConeTracks;
  std::vector<int> fShowerTracks;

  std::vector<int> fForwardSpacePoints;
  std::vector<int> fBackwardSpacePoints;
  std::vector<int> fCylinderSpacePoints;
  std::vector<int> fSphereSpacePoints;
  std::map<int,double> fIsolationSpacePoints;
  std::map<int,std::vector<int> > fRectangleHits;
  int fNumRectangleHits;
  double fMinForwardSpacePoint;
  double fMinBackwardSpacePoint;

  bool fShower;
  bool fShowerTrack;
  bool fShowerCone;
  bool fTrack;

};

class shower::TrackShowerSepAlg {
 public:

  TrackShowerSepAlg(fhicl::ParameterSet const& pset);

  /// Read in configurable parameters from provided parameter set
  void reconfigure(fhicl::ParameterSet const& pset);

  /// Run the tracks shower separation over the given hits
  /// Uses tracks and space points (and their associations) to assist
  void RunTrackShowerSep(int event,
			 const std::vector<art::Ptr<recob::Hit> >& hits,
			 const std::vector<art::Ptr<recob::Track> >& tracks,
			 const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
			 const art::FindManyP<recob::Hit>& fmht,
			 const art::FindManyP<recob::Track>& fmth,
			 const art::FindManyP<recob::SpacePoint>& fmspt,
			 const art::FindManyP<recob::Track>& fmtsp);

  /// Returns track-like tracks
  std::vector<art::Ptr<recob::Track> > TrackTracks();

  /// Returns shower-like hits
  std::vector<art::Ptr<recob::Hit> > ShowerHits();

  /// Return shower shower points
  std::vector<TVector3> ShowerStarts();

 private:

  ///
  std::vector<int> InitialTrackLikeSegment(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks);

  ///
  TVector2 Gradient(const std::vector<art::Ptr<recob::Hit> >& hits);

  ///
  TVector2 Gradient(const std::vector<art::Ptr<recob::Hit> >& hits, const std::unique_ptr<TVector2>& end);

  ///
  TVector3 Gradient(const std::vector<TVector3>& points, const std::unique_ptr<TVector3>& dir);

  ///
  TVector3 Gradient(const art::Ptr<recob::Track>& track);

  ///
  TVector3 Gradient(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints);

  /// Projects the 3D point given along the line defined by the specified direction
  /// Coordinate system is assumed to be centred at the origin unless a difference point is specified
  TVector3 ProjPoint(const TVector3& point, const TVector3& direction, const TVector3& origin = TVector3(0,0,0));

  /// Return 3D point of this space point
  TVector3 SpacePointPos(const art::Ptr<recob::SpacePoint>& spacePoint);

  ///
  double SpacePointsRMS(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints);

  // Copied from EMShower -- would be good to have them in a common place
  TVector2 HitCoordinates(art::Ptr<recob::Hit> const& hit);
  TVector2 HitPosition(art::Ptr<recob::Hit> const& hit);
  TVector2 HitPosition(TVector2 const& pos, geo::PlaneID planeID);
  double GlobalWire(const geo::WireID& wireID);
  TVector2 Project3DPointOntoPlane(TVector3 const& point, int plane, int cryostat = 0);






  // All the information about all reconstructed tracks
  std::map<int,std::unique_ptr<ReconTrack> > fReconTracks;

  // Output data products
  std::vector<art::Ptr<recob::Hit> >   fShowerHits;
  std::vector<art::Ptr<recob::Track> > fTrackTracks;
  std::vector<TVector3>                fShowerStarts;

  // Parameters
  int fDebug;
  std::string fDetector;

  // Properties
  double fConeAngle;
  double fCylinderRadius;
  double fRectangleWidth;

  // Cuts
  double fTrackVertexCut;
  double fCylinderCut;
  double fShowerConeCut;
  double fRectangleCut;
  double fLeastSquareCut;

  // art services
  art::ServiceHandle<geo::Geometry> fGeom;
  detinfo::DetectorProperties const* fDetProp;

};

#endif













 /*  // --------------------------- OLD (late 2015) ------------------------------- */

 /* public: */

 /*  /// Takes specific previously reconstructed quantites and removes hits which are considered track-like */
 /*  /// Returns a vector of hits which are determined to be shower-like */
 /*  std::vector<art::Ptr<recob::Hit> > RemoveTrackHits(const std::vector<art::Ptr<recob::Hit> >& initialHits, */
 /* 						     const std::vector<art::Ptr<recob::Track> >& tracks, */
 /* 						     const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints, */
 /* 						     const std::vector<art::Ptr<recob::Vertex> >& vertices, */
 /* 						     const art::FindManyP<recob::Track>& fmth, */
 /* 						     const art::FindManyP<recob::Track>& fmtsp, */
 /* 						     const art::FindManyP<recob::Hit>& fmh, */
 /* 						     int event, */
 /* 						     int run); */

 /*  /// Uses information from Pandora reconstruction to separate track-like and shower-like hits */
 /*  /// Returns a vector of hits which are determined to be shower-like */
 /*  std::vector<art::Ptr<recob::Hit> > RemoveTrackHits(const std::vector<art::Ptr<recob::Hit> >& hits, */
 /* 						     const std::vector<art::Ptr<recob::PFParticle> > pfParticles, */
 /* 						     const art::FindManyP<recob::Cluster>& fmc, */
 /* 						     const art::FindManyP<recob::Hit>& fmh); */


 /* private: */

 /*  /// Fill the output container with all the hits not associated with track-like objects */
 /*  std::vector<art::Ptr<recob::Hit> > FillHitsToCluster(const std::vector<art::Ptr<recob::Hit> >& initialHits, */
 /* 						       const art::FindManyP<recob::Track>& fmt); */

 /*  /// Find the true track most likely associated with this hit */
 /*  int FindTrackID(const art::Ptr<recob::Hit>& hit); */

 /*  /// Find the true track most likely associated with this set of hits */
 /*  int FindTrueTrack(const std::vector<art::Ptr<recob::Hit> >& trackHits); */

 /*  /// Finds the space points surrounding the track but not part of it */
 /*  std::vector<art::Ptr<recob::SpacePoint> > GetSurroundingSpacePoints(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints, */
 /* 								      const art::FindManyP<recob::Track>& fmt, */
 /* 								      unsigned int trackID); */

 /*  /// Look for space points near the track and within a narrow cone */
 /*  std::vector<art::Ptr<recob::SpacePoint> > GetSpacePointsInCone(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints, */
 /* 								 const TVector3& trackEnd, */
 /* 								 const TVector3& trackDirection); */

 /*  /// Determines whether or not a track is actually the start of a shower */
 /*  bool IdentifyShowerLikeTrack(const TVector3& end, */
 /* 			       const TVector3& direction, */
 /* 			       const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints); */

 /*  /// Attempt to identify tracks in the event by using the centre of the event */
 /*  void IdentifyTracksFromEventCentre(const std::vector<art::Ptr<recob::Track> >& tracks, */
 /* 				     const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints, */
 /* 				     const art::FindManyP<recob::Track>& fmtsp); */

 /*  /// Identifies tracks which start just after previously identified tracks end */
 /*  void IdentifyTracksNearTracks(const std::vector<art::Ptr<recob::Track> >& tracks); */

 /*  /// Identifies hadron-like tracks which originate from near the interaction vertex */
 /*  void IdentifyTracksNearVertex(const art::Ptr<recob::Vertex>& vertex, */
 /* 				const std::vector<art::Ptr<recob::Track> >& tracks, */
 /* 				const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints, */
 /* 				const art::FindManyP<recob::Track>& fmtsp); */

 /*  /// Finds the spread of a set of space points about their central axis */
 /*  double SpacePointSpread(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints); */

 /*  std::vector<int> fTrackLikeIDs, fShowerLikeIDs; */

 /*  // Configurable parameters */
 /*  double fAngleCut, fDistanceCut, fVertexProximityCut, fTrackProximityCut, fAvTrackHitDistance; */

 /*  art::ServiceHandle<cheat::BackTracker> backtracker; */
 /*  art::ServiceHandle<art::TFileService> tfs; */

 /*  TTree* ftree; */
 /*  double Distance, Angle, Length, AvDistance; */
 /*  int Event, Run, TrackID, pdg, NSpacePoints; */
