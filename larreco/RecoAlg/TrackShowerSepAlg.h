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
  struct TrackShowerSepParameters;
}

class shower::ReconTrack {
 public:

  ReconTrack(int id) {
    fID = id;

    fTrack = false;
    fShower = false;
    fShowerTrack = false;
    fShowerCone = false;

    fTrackLike = 0;
    fShowerLike = false;

    fNumHits = 0;
    fNumRectangleHits = 0;
  }

  // Setters -------------------------------------------------------------------------------------------------

  // Track properties
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

  // Determined properties
  void AddForwardConeSpacePoint(int spacePoint) { fForwardConeSpacePoints.push_back(spacePoint); }
  void AddBackwardConeSpacePoint(int spacePoint) { fBackwardConeSpacePoints.push_back(spacePoint); }
  void AddForwardShowerConeSpacePoint(int spacePoint) { fForwardShowerConeSpacePoints.push_back(spacePoint); }
  void AddBackwardShowerConeSpacePoint(int spacePoint) { fBackwardShowerConeSpacePoints.push_back(spacePoint); }
  void AddForwardConeTrack(int track) {
    if (std::find(fForwardConeTracks.begin(), fForwardConeTracks.end(), track) == fForwardConeTracks.end())
      fForwardConeTracks.push_back(track);
  }
  void AddBackwardConeTrack(int track) {
    if (std::find(fBackwardConeTracks.begin(), fBackwardConeTracks.end(), track) == fBackwardConeTracks.end())
      fBackwardConeTracks.push_back(track);
  }
  void AddForwardShowerConeTrack(int track) {
    if (std::find(fForwardShowerConeTracks.begin(), fForwardShowerConeTracks.end(), track) == fForwardShowerConeTracks.end())
      fForwardShowerConeTracks.push_back(track);
  }
  void AddBackwardShowerConeTrack(int track) {
    if (std::find(fBackwardShowerConeTracks.begin(), fBackwardShowerConeTracks.end(), track) == fBackwardShowerConeTracks.end())
      fBackwardShowerConeTracks.push_back(track);
  }
  void AddForwardConeSpacePointEnd(int spacePoint) { fForwardConeSpacePointsEnd.push_back(spacePoint); }
  void AddBackwardConeSpacePointEnd(int spacePoint) { fBackwardConeSpacePointsEnd.push_back(spacePoint); }
  void AddForwardShowerConeSpacePointEnd(int spacePoint) { fForwardShowerConeSpacePointsEnd.push_back(spacePoint); }
  void AddBackwardShowerConeSpacePointEnd(int spacePoint) { fBackwardShowerConeSpacePointsEnd.push_back(spacePoint); }
  void AddForwardConeTrackEnd(int track) {
    if (std::find(fForwardConeTracksEnd.begin(), fForwardConeTracksEnd.end(), track) == fForwardConeTracksEnd.end())
      fForwardConeTracksEnd.push_back(track);
  }
  void AddBackwardConeTrackEnd(int track) {
    if (std::find(fBackwardConeTracksEnd.begin(), fBackwardConeTracksEnd.end(), track) == fBackwardConeTracksEnd.end())
      fBackwardConeTracksEnd.push_back(track);
  }
  void AddForwardShowerConeTrackEnd(int track) {
    if (std::find(fForwardShowerConeTracksEnd.begin(), fForwardShowerConeTracksEnd.end(), track) == fForwardShowerConeTracksEnd.end())
      fForwardShowerConeTracksEnd.push_back(track);
  }
  void AddBackwardShowerConeTrackEnd(int track) {
    if (std::find(fBackwardShowerConeTracksEnd.begin(), fBackwardShowerConeTracksEnd.end(), track) == fBackwardShowerConeTracksEnd.end())
      fBackwardShowerConeTracksEnd.push_back(track);
  }

  void AddTrackConversion(int track) { fTrackConversions.push_back(track); }
  void AddShowerTrack(int track) { fShowerTracks.push_back(track); }
  void AddShowerCone(int track) { fShowerCones.push_back(track); }

  void AddCylinderSpacePoint(int spacePoint) { fCylinderSpacePoints.push_back(spacePoint); }
  void AddSphereSpacePoint(int spacePoint) { fSphereSpacePoints.push_back(spacePoint); }
  void AddIsolationSpacePoint(int spacePoint, double distance) { fIsolationSpacePoints[spacePoint] = distance; }
  void AddRectangleHit(const art::Ptr<recob::Hit>& hit) { fRectangleHits[hit->WireID().Plane].push_back(hit.key()); ++fNumRectangleHits; }

  // Getters -------------------------------------------------------------------------------------------------

  // Track properties
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

  // Determined space point properties
  int ConeSpacePointSize() const { return (int)fForwardConeSpacePoints.size() - (int)fBackwardConeSpacePoints.size(); }
  int NumForwardConeSpacePoints() const { return fForwardConeSpacePoints.size(); }
  int NumBackwardConeSpacePoints() const { return fBackwardConeSpacePoints.size(); }
  const std::vector<int>& ForwardConeSpacePoints() const { return fForwardConeSpacePoints; }
  const std::vector<int>& BackwardConeSpacePoints() const { return fBackwardConeSpacePoints; }
  int ShowerConeSpacePointSize() const { return (int)fForwardShowerConeSpacePoints.size() - (int)fBackwardShowerConeSpacePoints.size(); }
  int NumForwardShowerConeSpacePoints() const { return fForwardShowerConeSpacePoints.size(); }
  int NumBackwardShowerConeSpacePoints() const { return fBackwardShowerConeSpacePoints.size(); }
  const std::vector<int>& ForwardShowerConeSpacePoints() const { return fForwardShowerConeSpacePoints; }
  const std::vector<int>& BackwardShowerConeSpacePoints() const { return fBackwardShowerConeSpacePoints; }
  int ConeSpacePointSizeEnd() const { return (int)fForwardConeSpacePointsEnd.size() - (int)fBackwardConeSpacePointsEnd.size(); }
  int ShowerConeSpacePointSizeEnd() const { return (int)fForwardShowerConeSpacePointsEnd.size() - (int)fBackwardShowerConeSpacePointsEnd.size(); }

  // Determined track properties
  int ConeTrackSize() const { return (int)fForwardConeTracks.size() - (int)fBackwardConeTracks.size(); }
  int NumForwardConeTracks() const { return fForwardConeTracks.size(); }
  int NumBackwardConeTracks() const { return fBackwardConeTracks.size(); }
  const std::vector<int>& ForwardConeTracks() const { return fForwardConeTracks; }
  const std::vector<int>& BackwardConeTracks() const { return fBackwardConeTracks; }
  int ShowerConeTrackSize() const { return (int)fForwardShowerConeTracks.size() - (int)fBackwardShowerConeTracks.size(); }
  int NumForwardShowerConeTracks() const { return fForwardShowerConeTracks.size(); }
  int NumBackwardShowerConeTracks() const { return fBackwardShowerConeTracks.size(); }
  const std::vector<int>& ForwardShowerConeTracks() const { return fForwardShowerConeTracks; }
  const std::vector<int>& BackwardShowerConeTracks() const { return fBackwardShowerConeTracks; }
  int ConeTrackSizeEnd() const { return (int)fForwardConeTracksEnd.size() - (int)fBackwardConeTracksEnd.size(); }
  int ShowerConeTrackSizeEnd() const { return (int)fForwardShowerConeTracksEnd.size() - (int)fBackwardShowerConeTracksEnd.size(); }

  //
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

  //
  const std::vector<int>& ShowerTracks() const { return fShowerTracks; }
  const std::vector<int>& ShowerCones() const { return fShowerCones; }
  const std::vector<int>& TrackConversions() const { return fTrackConversions; }

  // Decision
  bool IsShower() const { return fShower; }
  bool IsShowerLike() const { return fShowerLike; }
  bool IsShowerTrack() const { return fShowerTrack; }
  bool IsShowerCone() const { return fShowerCone; }
  bool IsTrack() const { return fTrack; }
  bool IsTrackLike() const { return fTrackLike > 3; }
  bool IsUndetermined() const { return !fTrack and !fShower; }
  int Trackiness() const { return fTrackLike; }

  // Doers -------------------------------------------------------------------------------------------------
  void MakeShower() {
    fShowerLike = true;
    if (fTrack)
      this->MakeShowerTrack();
    else
      this->MakeShowerCone();
  }
  void MakeShowerTrack() {
    fShower = true;
    fShowerLike = true;
    fShowerTrack = true;
    fShowerCone = false;
    fTrack = false;
  }
  void MakeShowerCone() {
    fShower = true;
    fShowerLike = true;
    fShowerCone = true;
    fShowerTrack = false;
    fTrack = false;
  }
  void MakeShowerLike() { fShowerLike = true; }
  void MakeTrack() {
    fTrack = true;
    fShower = false;
    fShowerTrack = false;
    fShowerCone = false;
  }
  void MakeUndetermined() {
    fTrack = false;
    fShower = false;
    fShowerTrack = false;
    fShowerCone = false;
  }
  void IncreaseTrackiness() { ++fTrackLike; }

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
    std::vector<int> tmpVec;
    tmpVec = fForwardConeSpacePointsEnd;
    fForwardConeSpacePointsEnd = fForwardConeSpacePoints;
    fForwardConeSpacePoints = tmpVec;
    tmpVec = fBackwardConeSpacePointsEnd;
    fBackwardConeSpacePointsEnd = fBackwardConeSpacePoints;
    fBackwardConeSpacePoints = tmpVec;
    tmpVec = fForwardShowerConeSpacePointsEnd;
    fForwardShowerConeSpacePointsEnd = fForwardShowerConeSpacePoints;
    fForwardShowerConeSpacePoints = tmpVec;
    tmpVec = fBackwardShowerConeSpacePointsEnd;
    fBackwardShowerConeSpacePointsEnd = fBackwardShowerConeSpacePoints;
    fBackwardShowerConeSpacePoints = tmpVec;
    tmpVec = fForwardConeTracksEnd;
    fForwardConeTracksEnd = fForwardConeTracks;
    fForwardConeTracks = tmpVec;
    tmpVec = fBackwardConeTracksEnd;
    fBackwardConeTracksEnd = fBackwardConeTracks;
    fBackwardConeTracks = tmpVec;
    tmpVec = fForwardShowerConeTracksEnd;
    fForwardShowerConeTracksEnd = fForwardShowerConeTracks;
    fForwardShowerConeTracks = tmpVec;
    tmpVec = fBackwardShowerConeTracksEnd;
    fBackwardShowerConeTracksEnd = fBackwardShowerConeTracks;
    fBackwardShowerConeTracks = tmpVec;
  }

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

  std::vector<int> fTrackConversions;
  std::vector<int> fShowerTracks;
  std::vector<int> fShowerCones;

  std::vector<int> fForwardConeSpacePoints, fForwardConeSpacePointsEnd;
  std::vector<int> fBackwardConeSpacePoints, fBackwardConeSpacePointsEnd;
  std::vector<int> fForwardShowerConeSpacePoints, fForwardShowerConeSpacePointsEnd;
  std::vector<int> fBackwardShowerConeSpacePoints, fBackwardShowerConeSpacePointsEnd;
  std::vector<int> fForwardConeTracks, fForwardConeTracksEnd;
  std::vector<int> fBackwardConeTracks, fBackwardConeTracksEnd;
  std::vector<int> fForwardShowerConeTracks, fForwardShowerConeTracksEnd;
  std::vector<int> fBackwardShowerConeTracks, fBackwardShowerConeTracksEnd;

  std::vector<int> fCylinderSpacePoints;
  std::vector<int> fSphereSpacePoints;
  std::map<int,double> fIsolationSpacePoints;
  int fNumRectangleHits;
  std::map<int,std::vector<int> > fRectangleHits;

  bool fTrack;
  bool fShower;
  bool fShowerTrack;
  bool fShowerCone;

  int fTrackLike;
  bool fShowerLike;

};

struct shower::TrackShowerSepParameters {
  double AverageRectangleHits;
  double AverageCylinderSpacePoints;
  double AverageConeSize;
  TrackShowerSepParameters() {
    AverageRectangleHits = 0;
    AverageCylinderSpacePoints = 0;
    AverageConeSize = 0;
  }    
  TrackShowerSepParameters(double avRectangleHits, double avCylinderSpacePoints, double avConeSize) {
    AverageRectangleHits = avRectangleHits;
    AverageCylinderSpacePoints = avCylinderSpacePoints;
    AverageConeSize = avConeSize;
  }
  void SetAvRectangleHits(double avRectangleHits) { AverageRectangleHits = avRectangleHits; }
  void SetAvCylinderSpacePoints(double avCylinderSpacePoints) { AverageCylinderSpacePoints = avCylinderSpacePoints; }
  void SetAvConeSize(double avConeSize) { AverageConeSize = avConeSize; }
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
  void TrackProperties2D(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks,
			 TrackShowerSepParameters& parameters,
			 const std::vector<art::Ptr<recob::Hit> >& hits,
			 const art::FindManyP<recob::Track> fmth);

  ///
  void TrackProperties3D(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks,
			 TrackShowerSepParameters& parameters,
			 const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
			 const art::FindManyP<recob::Track>& fmtsp);

  ///
  void IdentifyTracks(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks, const TrackShowerSepParameters& parameters);

  ///
  void ConeProperties(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks,
		      TrackShowerSepParameters& parameters,
		      const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
		      const art::FindManyP<recob::Track>& fmtsp);

  ///
  void IdentifyShowers(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks, const TrackShowerSepParameters& parameters);

  ///
  void IdentifyShowerTracks(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks);

  ///
  void DowngradeTrack(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks, int track);

  ///
  void IdentifyShowerCones(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks);

  ///
  void ResolveUndeterminedTracks(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks);

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






  // Output data products
  // Note anything added here needs to be cleared at the start of each event
  std::vector<art::Ptr<recob::Hit> >   fShowerHits;
  std::vector<art::Ptr<recob::Track> > fTrackTracks;
  std::vector<TVector3>                fShowerStarts;

  // Parameters
  bool fRunTrackShowerSep;
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
