////////////////////////////////////////////////////////////////////////
// Class: TrackShowerSepAlg
// File:  TrackShowerSepAlg.cxx
// Author: Mike Wallbank (m.wallbank@sheffield.ac.uk), November 2015
//
// Track/shower separation class.
// Provides methods for removing hits associated with track-like
// objects.
// To be run after track reconstruction, before shower reconstruction.
////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/TrackShowerSepAlg.h"

shower::TrackShowerSepAlg::TrackShowerSepAlg(fhicl::ParameterSet const& pset) : fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>()) {
  this->reconfigure(pset);
}

void shower::TrackShowerSepAlg::reconfigure(fhicl::ParameterSet const& pset) {
  fConeAngle      = pset.get<double>("ConeAngle");
  fCylinderRadius = pset.get<double>("CylinderRadius");
  fRectangleWidth = pset.get<double>("RectangleWidth");
  fTrackVertexCut = pset.get<double>("TrackVertexCut");
  fCylinderCut    = pset.get<double>("CylinderCut");
  fShowerConeCut  = pset.get<double>("ShowerConeCut");
  fRectangleCut   = pset.get<double>("RectangleCut");
  fLeastSquareCut = pset.get<double>("LeastSquareCut");

  fDebug = pset.get<int>("Debug",0);
  fDetector = pset.get<std::string>("Detector");
}

void shower::TrackShowerSepAlg::RunTrackShowerSep(int event,
						  const std::vector<art::Ptr<recob::Hit> >& hits,
						  const std::vector<art::Ptr<recob::Track> >& tracks,
						  const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
						  const art::FindManyP<recob::Hit>& fmht,
						  const art::FindManyP<recob::Track>& fmth,
						  const art::FindManyP<recob::SpacePoint>& fmspt,
						  const art::FindManyP<recob::Track>& fmtsp) {

  // Ok, here we are again
  // Playing the game in which no one wins, but everyone loses
  // Trying to separate tracks from showers
  //    Ode to track/shower separation
  //    M Wallbank, Oct 2016

  // Showers are comprised of two topologically different parts:
  //  -- the 'shower track' (the initial part of the shower)
  //  -- the 'shower cone' (the part of the shower after the initial track when showering happens)

  //                                            -
  //                                          --
  //                 -             --       --
  //               --         ---         ---
  //             --      ----       --
  //  -----------   -----       -  --        ---
  //                    ---             ---
  //                       --               ---
  //  {=========}{==============================}
  // shower track          shower cone

  fReconTracks.clear();
  fShowerHits.clear();
  fTrackTracks.clear();

  // Save info about the tracks
  for (std::vector<art::Ptr<recob::Track> >::const_iterator trackIt = tracks.begin(); trackIt != tracks.end(); ++trackIt) {

    std::unique_ptr<ReconTrack> track = std::make_unique<ReconTrack>(trackIt->key());

    // 3D
    track->SetVertex((*trackIt)->Vertex());
    track->SetEnd((*trackIt)->End());
    track->SetVertexDir((*trackIt)->VertexDirection());
    track->SetLength((*trackIt)->Length());
    track->SetDirection3D(Gradient(*trackIt));
    track->SetHits(fmht.at(trackIt->key()));
    track->SetSpacePoints(fmspt.at(trackIt->key()));

    // 2D
    const std::map<int,std::vector<art::Ptr<recob::Hit> > >& trackHits = track->Hits();
    for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator hitIt = trackHits.begin(); hitIt != trackHits.end(); ++hitIt) {
      track->AddPlane(hitIt->first);
      track->SetVertex2D(hitIt->first, Project3DPointOntoPlane(track->Vertex(), hitIt->first));
      track->SetEnd2D(hitIt->first, Project3DPointOntoPlane(track->End(), hitIt->first));
      track->SetDirection2D(hitIt->first, Gradient(hitIt->second, std::make_unique<TVector2>(track->End2D(hitIt->first)-track->Vertex2D(hitIt->first))));
    }

    // const std::vector<art::Ptr<recob::SpacePoint> > spsss = fmspt.at(trackIt->key());
    // std::cout << "Track " << trackIt->key() << " has " << spsss.size() << " space points and " << (*trackIt)->NumberTrajectoryPoints() << " traj points" << std::endl;
    // if (trackIt->key() == 5)
    //   for (unsigned int i = 0; i < (*trackIt)->NumberTrajectoryPoints(); ++i)
    // 	std::cout << "Traj point " << i << " has position (" << (*trackIt)->LocationAtPoint(i).X() << ", " << (*trackIt)->LocationAtPoint(i).Y() << ", " << (*trackIt)->LocationAtPoint(i).Z() << ")" << std::endl;

    fReconTracks[trackIt->key()] = std::move(track);

  }

  // std::vector<int> showerLikeTracks, trackLikeTracks;
  // std::vector<int> showerTracks = InitialTrackLikeSegment(reconTracks);

  // Consider the hit rectangle situation
  double avRectangleHits = 0;
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = fReconTracks.begin(); trackIt != fReconTracks.end(); ++trackIt) {
    for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
      const std::vector<art::Ptr<recob::Track> >& hitTracks = fmth.at(hitIt->key());
      if (find_if(hitTracks.begin(), hitTracks.end(), [&trackIt](const art::Ptr<recob::Track>& t){ return (int)t.key() == trackIt->first; }) != hitTracks.end())
	continue;
      TVector2 point = trackIt->second->Vertex2D((*hitIt)->WireID().Plane);
      TVector2 end = trackIt->second->End2D((*hitIt)->WireID().Plane);
      TVector2 direction = trackIt->second->Direction2D((*hitIt)->WireID().Plane);
      TVector2 pos = HitPosition(*hitIt);
      TVector2 proj = (pos-point).Proj(direction) + point;
      if ((pos-proj).Mod() < fRectangleWidth)//  and
	  // (pos-point)*direction > 0 and
	  // (pos-point)*direction < (end-point)*direction)
	trackIt->second->AddRectangleHit(*hitIt);
    }
    avRectangleHits += trackIt->second->RectangleHitRatio();
  }
  avRectangleHits /= (double)fReconTracks.size();

  if (fDebug > 1) {
    std::cout << std::endl << "Rectangle hit ratios:" << std::endl;
    for (std::map<int,std::unique_ptr<ReconTrack> >::const_iterator trackIt = fReconTracks.begin(); trackIt != fReconTracks.end(); ++trackIt) {
      std::cout << "  Track " << trackIt->first << " has " << trackIt->second->NumRectangleHits() << ", ratio " << trackIt->second->RectangleHitRatio() << " (event average ratio " << avRectangleHits << ")" << std::endl;
      // const std::map<int,std::vector<art::Ptr<recob::Hit> > >& trackHits = trackIt->second->Hits();
      // for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator hitIt = trackHits.begin(); hitIt != trackHits.end(); ++hitIt)
      //   std::cout << "    Plane " << hitIt->first << " has " << trackIt->second->NumRectangleHits(hitIt->first) << " (ratio " << trackIt->second->RectangleHitRatio(hitIt->first) << ")" << std::endl;
    }
  }

  // Consider the space point cylinder situation
  double avCylinderSpacePoints = 0;
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = fReconTracks.begin(); trackIt != fReconTracks.end(); ++trackIt) {
    // Get the 3D properties of the track
    TVector3 point = trackIt->second->Vertex();
    TVector3 direction = trackIt->second->Direction3D();
    // if (trackIt->second->Vertex().X() > 250 and trackIt->second->Vertex().X() < 252 and
    // 	trackIt->second->Vertex().Y() > -440 and trackIt->second->Vertex().Y() < -430 and
    // 	trackIt->second->Vertex().Z() > 1080 and trackIt->second->Vertex().Z() < 1090)
    //   std::cout << "Track " << trackIt->first << " begins at the most upstream vertex" << std::endl;
    // if (trackIt->second->Vertex().X() > 254 and trackIt->second->Vertex().X() < 258 and
    // 	trackIt->second->Vertex().Y() > -440 and trackIt->second->Vertex().Y() < -420 and
    // 	trackIt->second->Vertex().Z() > 1115 and trackIt->second->Vertex().Z() < 1130)
    //   std::cout << "Track " << trackIt->first << " begins at the supposed vertex" << std::endl;
    // if (trackIt->second->End().X() > 254 and trackIt->second->End().X() < 258 and
    // 	trackIt->second->End().Y() > -440 and trackIt->second->End().Y() < -420 and
    // 	trackIt->second->End().Z() > 1115 and trackIt->second->End().Z() < 1130)
    //   std::cout << "Track " << trackIt->first << " ends at the supposed vertex" << std::endl;
    // std::cout << "Track " << trackIt->first << " has vertex (" << trackIt->second->Vertex().X() << ", " << trackIt->second->Vertex().Y() << ", " << trackIt->second->Vertex().Z() << ") and end (" << trackIt->second->End().X() << ", " << trackIt->second->End().Y() << ", " << trackIt->second->End().Z() << "), with vertex direction (" << trackIt->second->VertexDirection().X() << ", " << trackIt->second->VertexDirection().Y() << ", " << trackIt->second->VertexDirection().Z() << ")" << std::endl;
    // Count space points in the volume around the track
    for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt) {
      const std::vector<art::Ptr<recob::Track> > spTracks = fmtsp.at(spacePointIt->key());
      if (find_if(spTracks.begin(), spTracks.end(), [&trackIt](const art::Ptr<recob::Track>& t){ return (int)t.key() == trackIt->first; }) != spTracks.end())
	continue;
      // Get the properties of this space point
      TVector3 pos = SpacePointPos(*spacePointIt);
      TVector3 proj = ProjPoint(pos, direction, point);
      if ((pos-proj).Mag() < fCylinderRadius)
	trackIt->second->AddCylinderSpacePoint(spacePointIt->key());
      // if ((pos-proj).Mag() < fCylinderRadius and
      // 	  (pos-point)*direction > 0 and
      // 	  (pos-point)*direction < trackIt->second->Length())
      // 	std::cout << "Space point " << spacePointIt->key() << " (" << pos.X() << ", " << pos.Y() << ", " << pos.Z() << ") in cylinder around track " << trackIt->first << " (assocatied with track " << spTracks.at(0).key() << "); point is (" << point.X() << ", " << point.Y() << ", " << point.Z() << "), proj end (" << ((trackIt->second->Length()*trackIt->second->VertexDirection())+point).X() << ", " << ((trackIt->second->Length()*trackIt->second->VertexDirection())+point).Y() << ", " << ((trackIt->second->Length()*trackIt->second->VertexDirection())+point).Z() << ")" << std::endl;
    }
    avCylinderSpacePoints += trackIt->second->CylinderSpacePointRatio();
  }
  avCylinderSpacePoints /= (double)fReconTracks.size();

  if (fDebug > 1) {
    std::cout << std::endl << "Cylinder space point ratios:" << std::endl;
    for (std::map<int,std::unique_ptr<ReconTrack> >::const_iterator trackIt = fReconTracks.begin(); trackIt != fReconTracks.end(); ++trackIt)
      std::cout << "  Track " << trackIt->first << " has cylinder space point ratio " << trackIt->second->CylinderSpacePointRatio() << " (event average " << avCylinderSpacePoints << ")" << std::endl;
  }

  double avCylinderRectangle = 0;
  for (std::map<int,std::unique_ptr<ReconTrack> >::const_iterator trackIt = fReconTracks.begin(); trackIt != fReconTracks.end(); ++trackIt)
    avCylinderRectangle += trackIt->second->CylinderSpacePointRatio()*trackIt->second->RectangleHitRatio();
  avCylinderRectangle /= (double)fReconTracks.size();

  if (fDebug > 1) {
    std::cout << std::endl << "Combinded cylinder/rectangle:" << std::endl;
    for (std::map<int,std::unique_ptr<ReconTrack> >::const_iterator trackIt = fReconTracks.begin(); trackIt != fReconTracks.end(); ++trackIt)
      std::cout << "  Track " << trackIt->first << ": " << (trackIt->second->CylinderSpacePointRatio()) * (trackIt->second->RectangleHitRatio()) << "\tEvent average is " << avCylinderRectangle << std::endl;
  }

  if (fDebug > 1)
    std::cout << std::endl << "Looking at track straightness:" << std::endl;
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = fReconTracks.begin(); trackIt != fReconTracks.end(); ++trackIt) {
    const std::vector<art::Ptr<recob::SpacePoint> >& trackPoints = trackIt->second->SpacePoints();
    double least_sq = 0;
    for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = trackPoints.begin(); spacePointIt != trackPoints.end(); ++spacePointIt) {
      TVector3 pos = SpacePointPos(*spacePointIt);
      TVector3 proj = ProjPoint(pos, trackIt->second->Direction3D(), trackIt->second->Centre3D());
      least_sq += TMath::Power((proj-pos).Mag(),2);
    }
    least_sq /= (double)((int)trackPoints.size()+2);
    trackIt->second->SetLeastSquareNDOF(least_sq);
    if (fDebug > 1)
      std::cout << "  Track " << trackIt->first << " has least squared / dof " << least_sq << std::endl;
  }

  // Identify tracks
  if (fDebug > 0)
    std::cout << std::endl << "Identifying tracks:" << std::endl;
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = fReconTracks.begin(); trackIt != fReconTracks.end(); ++trackIt)
    if (trackIt->second->CylinderSpacePointRatio() / avCylinderSpacePoints < fCylinderCut and
	trackIt->second->RectangleHitRatio() / avRectangleHits < fRectangleCut and
	trackIt->second->LeastSquareNDOF() < fLeastSquareCut) {
      if (fDebug > 0)
	std::cout << "  Making track " << trackIt->first << " a track (Type I) (" << trackIt->second->CylinderSpacePointRatio()/avCylinderSpacePoints << " < " << fCylinderCut << ", " << trackIt->second->RectangleHitRatio()/avRectangleHits << " < " << fRectangleCut << ", " << trackIt->second->LeastSquareNDOF() << " < " << fLeastSquareCut << ")" << std::endl;
      trackIt->second->MakeTrack();
    }

  // for (std::map<int,std::unique_ptr<ReconTrack> >::const_iterator trackIt = fReconTracks.begin(); trackIt != fReconTracks.end(); ++trackIt) {
  //   if (!trackIt->second->IsTrack())
  //     continue;
  //   std::cout << "Track " << trackIt->first << " (tagged as track) is this close to other tracks:" << std::endl;
  //   for (std::map<int,std::unique_ptr<ReconTrack> >::const_iterator otherTrackIt = fReconTracks.begin(); otherTrackIt != fReconTracks.end(); ++otherTrackIt) {
  //     if (trackIt->first == otherTrackIt->first)
  // 	continue;
  //     std::cout << "  Other track " << otherTrackIt->first << " from track " << trackIt->first << " vertex: " << (otherTrackIt->second->Vertex()-trackIt->second->Vertex()).Mag() << " (vertex) and " << (otherTrackIt->second->End()-trackIt->second->Vertex()).Mag() << " (end); end: " << (otherTrackIt->second->Vertex()-trackIt->second->End()).Mag() << " (vertex) and " << (otherTrackIt->second->End()-trackIt->second->End()).Mag() << " (end)" << std::endl;
  //   }
  // }

  // Identify further tracks
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = fReconTracks.begin(); trackIt != fReconTracks.end(); ++trackIt) {
    if (trackIt->second->IsTrack())
      continue;
    for (std::map<int,std::unique_ptr<ReconTrack> >::const_iterator otherTrackIt = fReconTracks.begin(); otherTrackIt != fReconTracks.end(); ++otherTrackIt) {
      if (trackIt->first == otherTrackIt->first or !otherTrackIt->second->IsTrack())
  	continue;
      if (fDebug > 3) {
	std::cout << "Distances from tracks to tracks:" << std::endl;
	double minDist = TMath::Min((trackIt->second->Vertex()-otherTrackIt->second->Vertex()).Mag(),
				    TMath::Min((trackIt->second->Vertex()-otherTrackIt->second->End()).Mag(),
					       TMath::Min((trackIt->second->End()-otherTrackIt->second->Vertex()).Mag(),
							  (trackIt->second->End()-otherTrackIt->second->End()).Mag())));

	std::cout << "  Min distance between track " << trackIt->first << " and " << otherTrackIt->first << " is " << minDist << std::endl;
      }
      if ((trackIt->second->Vertex() - otherTrackIt->second->Vertex()).Mag() < fTrackVertexCut or
  	  (trackIt->second->Vertex() - otherTrackIt->second->End()).Mag() < fTrackVertexCut or
  	  (trackIt->second->End() - otherTrackIt->second->Vertex()).Mag() < fTrackVertexCut or
  	  (trackIt->second->End() - otherTrackIt->second->End()).Mag() < fTrackVertexCut) {
	if (fDebug > 0)
	  std::cout << "  Making track " << trackIt->first << " a track (Type II) (close to " << otherTrackIt->first << ")" << std::endl;
  	trackIt->second->MakeTrack();
      }
    }
  }

  // Consider removing false tracks by looking at their closest approach to any other track

  // Consider the space point cone situation
  std::vector<art::Ptr<recob::SpacePoint> > showerSpacePoints;
  for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt) {
    bool showerSpacePoint = true;
    const std::vector<art::Ptr<recob::Track> > spacePointTracks = fmtsp.at(spacePointIt->key());
    for (std::vector<art::Ptr<recob::Track> >::const_iterator trackIt = spacePointTracks.begin(); trackIt != spacePointTracks.end(); ++trackIt)
      if (fReconTracks[trackIt->key()]->IsTrack())
	showerSpacePoint = false;
    if (showerSpacePoint)
      showerSpacePoints.push_back(*spacePointIt);
  }

  // Identify tracks which slipped through and shower tracks
  // For the moment, until the track tagging gets better at least, don't try to identify tracks from this
  double avConeSize = 0;
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = fReconTracks.begin(); trackIt != fReconTracks.end(); ++trackIt) {
    for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = showerSpacePoints.begin(); spacePointIt != showerSpacePoints.end(); ++spacePointIt) {
      bool associatedSpacePoint = false;
      const std::vector<art::Ptr<recob::Track> > spTracks = fmtsp.at(spacePointIt->key());
      for (std::vector<art::Ptr<recob::Track> >::const_iterator spTrackIt = spTracks.begin(); spTrackIt != spTracks.end(); ++spTrackIt)
	if ((int)spTrackIt->key() == trackIt->first)
	  associatedSpacePoint = true;
      if (associatedSpacePoint)
	continue;
      if ((SpacePointPos(*spacePointIt) - trackIt->second->Vertex()).Angle(trackIt->second->Direction3D()) < fConeAngle * TMath::Pi() / 180) {
	//std::cout << "Space point at (" << SpacePointPos(*spacePointIt).X() << ", " << SpacePointPos(*spacePointIt).Y() << ", " << SpacePointPos(*spacePointIt).Z() << ") is " << (SpacePointPos(*spacePointIt)-trackIt->second->End())*trackIt->second->Direction3D() << " cm away from track " << trackIt->first << " end" << std::endl;
	trackIt->second->AddForwardSpacePoint(spacePointIt->key(), (SpacePointPos(*spacePointIt)-trackIt->second->End())*trackIt->second->Direction3D());
	trackIt->second->AddForwardTrack(spTracks.at(0).key());
      }
      if ((SpacePointPos(*spacePointIt) - trackIt->second->Vertex()).Angle(-1*trackIt->second->Direction3D()) < fConeAngle * TMath::Pi() / 180) {
	//std::cout << "Space point at (" << SpacePointPos(*spacePointIt).X() << ", " << SpacePointPos(*spacePointIt).Y() << ", " << SpacePointPos(*spacePointIt).Z() << ") is " << (SpacePointPos(*spacePointIt)-trackIt->second->Vertex())*(-1*trackIt->second->Direction3D()) << " cm away from track " << trackIt->first << " vertex" << std::endl;
	trackIt->second->AddBackwardSpacePoint(spacePointIt->key(), (SpacePointPos(*spacePointIt)-trackIt->second->Vertex())*(-1*trackIt->second->Direction3D()));
	trackIt->second->AddBackwardTrack(spTracks.at(0).key());
      }
    }
    avConeSize += trackIt->second->ConeSize();
  }
  avConeSize /= (double)fReconTracks.size();
  if (fDebug > 0)
    std::cout << std::endl << "Identifying showers:" << std::endl;
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = fReconTracks.begin(); trackIt != fReconTracks.end(); ++trackIt) {
    // if (trackIt->second->ForwardSpacePoints() == 0)
    //   trackIt->second->MakeTrack();
    // double distanceFromAverage = (trackIt->second->ConeSize() - avConeSize) / TMath::Abs(avConeSize);
    // if (fDebug > 1)
    //   std::cout << "    Track " << trackIt->first << " has cone size " << trackIt->second->ConeSize() << " (forward " << trackIt->second->ForwardSpacePoints() << ") and tracks " << trackIt->second->TrackConeSize() << " (average " << avConeSize << ", distance from average " << distanceFromAverage << ")" << std::endl;
    // if (distanceFromAverage > fShowerConeCut) {
    //   trackIt->second->MakeShower();
    //   if (fDebug > 0)
    // 	std::cout << "  Making track " << trackIt->first << " a shower (Type I)" << std::endl;
    // }
    if (fDebug > 1)
      std::cout << "    Track " << trackIt->first << " has space point cone size " << trackIt->second->ConeSize() << " and track cone size " << trackIt->second->TrackConeSize() << std::endl;
    if (TMath::Abs(trackIt->second->ConeSize()) > 30 and TMath::Abs(trackIt->second->TrackConeSize()) > 3) {
      if (trackIt->second->ConeSize() < 0)
	trackIt->second->FlipTrack();
      // std::cout << "    Track " << trackIt->first << " has min forward space point distance " << trackIt->second->MinForwardSpacePointDistance() << std::endl;
      // if (TMath::Abs(trackIt->second->MinForwardSpacePointDistance() < 10)) {
	trackIt->second->MakeShower();
	if (fDebug > 0)
	  std::cout << "  Making track " << trackIt->first << " a shower (Type I)" << std::endl;
	//}
    }
  }

  // Look for shower cones
  std::cout << std::endl << "  Shower cones:" << std::endl;
  for (std::map<int,std::unique_ptr<ReconTrack> >::const_iterator trackIt = fReconTracks.begin(); trackIt != fReconTracks.end(); ++trackIt) {
    if (trackIt->second->IsShower()) {
      if (fDebug > 1)
	std::cout << "    Track " << trackIt->first << std::endl;
      for (std::map<int,std::unique_ptr<ReconTrack> >::iterator otherTrackIt = fReconTracks.begin(); otherTrackIt != fReconTracks.end(); ++otherTrackIt) {
	if (trackIt->first == otherTrackIt->first or !otherTrackIt->second->IsUndetermined())
	  continue;
	if (fDebug > 2)
	  std::cout << "      Other track " << otherTrackIt->first << " has angle from vertex " << (otherTrackIt->second->Vertex()-trackIt->second->Vertex()).Angle(trackIt->second->Direction3D()) * 180/TMath::Pi() << " and angle from end " << (otherTrackIt->second->End()-trackIt->second->Vertex()).Angle(trackIt->second->Direction3D()) * 180/TMath::Pi() << std::endl;
        if ((otherTrackIt->second->Vertex()-trackIt->second->Vertex()).Angle(trackIt->second->Direction3D()) < fConeAngle * TMath::Pi() / 180 or
	    (otherTrackIt->second->End()-trackIt->second->Vertex()).Angle(trackIt->second->Direction3D()) < fConeAngle * TMath::Pi() / 180) {
	  std::cout << "      " << otherTrackIt->first << std::endl;
	  otherTrackIt->second->MakeShower();
	  if (fDebug > 0)
	    std::cout << "      Making track " << otherTrackIt->first << " a shower (Type II)" << std::endl;
	}
      }
    }
  }

  // Look at remaining undetermined tracks
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = fReconTracks.begin(); trackIt != fReconTracks.end(); ++trackIt) {
    if (trackIt->second->IsUndetermined()) {
      // std::cout << "Remaining undetermined track " << trackIt->first << std::endl;
      // std::cout << "Length is " << trackIt->second->Length() << " and rms of space points is " << SpacePointsRMS(trackIt->second->SpacePoints()) << std::endl;
      trackIt->second->MakeShower();
    }
  }

  std::cout << std::endl << "Event " << event << " track shower separation:" << std::endl;
  std::cout << "Shower initial tracks are:" << std::endl;
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = fReconTracks.begin(); trackIt != fReconTracks.end(); ++trackIt)
    if (trackIt->second->IsShowerTrack())
      std::cout << "  " << trackIt->first << std::endl;

  std::cout << "Track tracks are:" << std::endl;
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = fReconTracks.begin(); trackIt != fReconTracks.end(); ++trackIt)
    if (trackIt->second->IsTrack())
      std::cout << "  " << trackIt->first << "\t\t\t\t\tStart (" << trackIt->second->Vertex().X() << ", " << trackIt->second->Vertex().Y() << ", " << trackIt->second->Vertex().Z() << "), end (" << trackIt->second->End().X() << ", " << trackIt->second->End().Y() << ", " << trackIt->second->End().Z() << "), length " << trackIt->second->Length() << std::endl;

  // Select all hits which aren't associated with a determined track
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    bool showerHit = true;
    const std::vector<art::Ptr<recob::Track> > hitTracks = fmth.at(hitIt->key());
    for (std::vector<art::Ptr<recob::Track> >::const_iterator hitTrackIt = hitTracks.begin(); hitTrackIt != hitTracks.end(); ++hitTrackIt)
      if (fReconTracks[hitTrackIt->key()]->IsTrack())
	showerHit = false;
    if (showerHit)
      fShowerHits.push_back(*hitIt);
  }

  // Select all tracks which were determined to be tracks
  for (std::vector<art::Ptr<recob::Track> >::const_iterator trackIt = tracks.begin(); trackIt != tracks.end(); ++trackIt)
    if (fReconTracks[trackIt->key()]->IsTrack())
      fTrackTracks.push_back(*trackIt);

  // Save shower starts
  for (std::map<int,std::unique_ptr<ReconTrack> >::const_iterator trackIt = fReconTracks.begin(); trackIt != fReconTracks.end(); ++trackIt)
    if (trackIt->second->IsShowerTrack())
      fShowerStarts.push_back(trackIt->second->Vertex());

}

std::vector<int> shower::TrackShowerSepAlg::InitialTrackLikeSegment(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks) {

  // Consider the cones for each track
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt) {

    for (std::map<int,std::unique_ptr<ReconTrack> >::iterator otherTrackIt = reconTracks.begin(); otherTrackIt != reconTracks.end(); ++otherTrackIt) {

      if (trackIt->first == otherTrackIt->first)
	continue;
      if ((otherTrackIt->second->Vertex() - trackIt->second->Vertex()).Angle(trackIt->second->VertexDirection()) < fConeAngle * TMath::Pi() / 180 or
	  (otherTrackIt->second->End() - trackIt->second->Vertex()).Angle(trackIt->second->VertexDirection()) < fConeAngle * TMath::Pi() / 180) {
	trackIt->second->AddForwardTrack(otherTrackIt->first);
	otherTrackIt->second->AddShowerTrack(trackIt->first);
      }
      if ((otherTrackIt->second->Vertex() - trackIt->second->Vertex()).Angle(-1*trackIt->second->VertexDirection()) < fConeAngle * TMath::Pi() / 180 or
	  (otherTrackIt->second->End() - trackIt->second->Vertex()).Angle(-1*trackIt->second->VertexDirection()) < fConeAngle * TMath::Pi() / 180)
	trackIt->second->AddBackwardTrack(otherTrackIt->first);

    }

  }

  // Determine if any of these tracks are actually shower tracks
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt) {

    if (!trackIt->second->ShowerTrackCandidate())
      continue;

    std::cout << "Track " << trackIt->first << " is a candidate, with shower cone tracks:" << std::endl;
    const std::vector<int>& showerConeTracks = trackIt->second->ForwardConeTracks();
    for (std::vector<int>::const_iterator showerConeTrackIt = showerConeTracks.begin(); showerConeTrackIt != showerConeTracks.end(); ++showerConeTrackIt)
      std::cout << "  " << *showerConeTrackIt << std::endl;

    bool isBestCandidate = true;
    const std::vector<int>& showerTracks = trackIt->second->ShowerTracks();
    for (std::vector<int>::const_iterator showerTrackIt = showerTracks.begin(); showerTrackIt != showerTracks.end(); ++showerTrackIt) {
      if (!reconTracks[*showerTrackIt]->ShowerTrackCandidate())
	continue;
      if (std::find(showerConeTracks.begin(), showerConeTracks.end(), *showerTrackIt) == showerConeTracks.end())
	continue;
      if (reconTracks[*showerTrackIt]->IsShowerTrack())
	continue;
      if (trackIt->second->TrackConeSize() < reconTracks[*showerTrackIt]->TrackConeSize())
	isBestCandidate = false;
    }

    if (isBestCandidate)
      trackIt->second->MakeShowerTrack();

  }

  // Determine which tracks are shower cones
  std::vector<int> showerTracks;
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt) {
    if (trackIt->second->IsShowerTrack()) {
      showerTracks.push_back(trackIt->first);
      const std::vector<int>& coneTracks = trackIt->second->ForwardConeTracks();
      for (std::vector<int>::const_iterator coneTrackIt = coneTracks.begin(); coneTrackIt != coneTracks.end(); ++coneTrackIt)
	reconTracks[*coneTrackIt]->MakeShowerCone();
    }
  }

  return showerTracks;

}

TVector2 shower::TrackShowerSepAlg::Gradient(const std::vector<art::Ptr<recob::Hit> >& hits, const std::unique_ptr<TVector2>& end) {

  TVector2 pos;
  int nhits = 0;
  double sumx=0., sumy=0., sumx2=0., sumxy=0.;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hit = hits.begin(); hit != hits.end(); ++hit) {
    ++nhits;
    pos = HitPosition(*hit);
    sumx += pos.X();
    sumy += pos.Y();
    sumx2 += pos.X() * pos.X();
    sumxy += pos.X() * pos.Y();
  }
  double gradient = (nhits * sumxy - sumx * sumy) / (nhits * sumx2 - sumx * sumx);
  TVector2 direction = TVector2(1,gradient).Unit();

  if (end and (*end)*direction < 0)
    direction *= -1;

  return direction;

}

TVector2 shower::TrackShowerSepAlg::Gradient(const std::vector<art::Ptr<recob::Hit> >& hits) {

  std::unique_ptr<TVector2> end;

  return Gradient(hits, end);

}

TVector3 shower::TrackShowerSepAlg::Gradient(const std::vector<TVector3>& points, const std::unique_ptr<TVector3>& dir) {

  int nhits = 0;
  double sumx=0., sumy=0., sumz=0., sumx2=0., sumy2=0., sumxy=0., sumxz=0., sumyz=0.;
  for (std::vector<TVector3>::const_iterator pointIt = points.begin(); pointIt != points.end(); ++pointIt) {
    ++nhits;
    sumx += pointIt->X();
    sumy += pointIt->Y();
    sumz += pointIt->Z();
    sumx2 += pointIt->X() * pointIt->X();
    sumy2 += pointIt->Y() * pointIt->Y();
    sumxy += pointIt->X() * pointIt->Y();
    sumxz += pointIt->X() * pointIt->Z();
    sumyz += pointIt->Y() * pointIt->Z();
  }

  double dydx = (nhits * sumxy - sumx * sumy) / (nhits * sumx2 - sumx * sumx);
  double yint = (sumy * sumx2 - sumx * sumxy) / (nhits * sumx2 - sumx * sumx);
  double dzdx = (nhits * sumxz - sumx * sumz) / (nhits * sumx2 - sumx * sumx);
  double zint = (sumz * sumx2 - sumx * sumxz) / (nhits * sumx2 - sumx * sumx);
  TVector2 directionXY = TVector2(1,dydx).Unit(), directionXZ = TVector2(1,dzdx).Unit();
  TVector3 direction = TVector3(1,dydx,dzdx).Unit();
  TVector3 intercept = TVector3(0,yint,zint);

  // Make sure the best fit direction is pointing correctly
  if (dir and TMath::Abs(direction.Angle(*dir)) > TMath::Pi() / 2.)
    direction *= -1;

  return direction;

}

TVector3 shower::TrackShowerSepAlg::Gradient(const art::Ptr<recob::Track>& track) {

  std::vector<TVector3> points;
  std::unique_ptr<TVector3> dir;

  for (unsigned int traj = 0; traj < track->NumberTrajectoryPoints(); ++traj)
    points.push_back(track->LocationAtPoint(traj));
  dir = std::make_unique<TVector3>(track->VertexDirection());

  return Gradient(points, dir);

}

TVector3 shower::TrackShowerSepAlg::Gradient(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints) {

  std::vector<TVector3> points;
  std::unique_ptr<TVector3> dir;

  for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt)
    points.push_back(SpacePointPos(*spacePointIt));

  return Gradient(points, dir);

}

TVector3 shower::TrackShowerSepAlg::ProjPoint(const TVector3& point, const TVector3& direction, const TVector3& origin) {
  return (point-origin).Dot(direction) * direction + origin;
}

TVector3 shower::TrackShowerSepAlg::SpacePointPos(const art::Ptr<recob::SpacePoint>& spacePoint) {
  const double* xyz = spacePoint->XYZ();
  return TVector3(xyz[0], xyz[1], xyz[2]);
}

double shower::TrackShowerSepAlg::SpacePointsRMS(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints) {

  TVector3 point = SpacePointPos(spacePoints.at(0));
  TVector3 direction = Gradient(spacePoints);

  std::vector<double> distances;
  for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt) {
    TVector3 pos = SpacePointPos(*spacePointIt);
    TVector3 proj = ProjPoint(pos, direction, point);
    distances.push_back((pos-proj).Mag());
  }

  double rms = TMath::RMS(distances.begin(), distances.end());

  return rms;

}

std::vector<art::Ptr<recob::Hit> > shower::TrackShowerSepAlg::ShowerHits() {

  if (fReconTracks.size() == 0)
    std::cout << "Warning! (TrackShowerSepAlg): no shower hits found. (Did you RunTrackShowerSep first?)" << std::endl;

  return fShowerHits;

}

std::vector<TVector3> shower::TrackShowerSepAlg::ShowerStarts() {

  if (fReconTracks.size() == 0)
    std::cout << "Warning! (TrackShowerSepAlg): no shower starts found. (Did you RunTrackShowerSep first?)" << std::endl;

  return fShowerStarts;  

}

std::vector<art::Ptr<recob::Track> > shower::TrackShowerSepAlg::TrackTracks() {

  if (fReconTracks.size() == 0)
    std::cout << "Warning! (TrackShowerSepAlg): no track tracks found. (Did you RunTrackShowerSep first?)" << std::endl;

  return fTrackTracks;

}



///////////////////////////////////////////////////////////////////////////////////
// Copied from EMShower -- would be good to have these in a common place for use

TVector2 shower::TrackShowerSepAlg::HitCoordinates(art::Ptr<recob::Hit> const& hit) {

  return TVector2(GlobalWire(hit->WireID()), hit->PeakTime());

}

TVector2 shower::TrackShowerSepAlg::HitPosition(art::Ptr<recob::Hit> const& hit) {

  geo::PlaneID planeID = hit->WireID().planeID();

  return HitPosition(HitCoordinates(hit), planeID);

}

TVector2 shower::TrackShowerSepAlg::HitPosition(TVector2 const& pos, geo::PlaneID planeID) {

  return TVector2(pos.X() * fGeom->WirePitch(planeID),
		  fDetProp->ConvertTicksToX(pos.Y(), planeID));

}

double shower::TrackShowerSepAlg::GlobalWire(const geo::WireID& wireID) {

  double globalWire = -999;

  // Induction
  if (fGeom->SignalType(wireID) == geo::kInduction) {
    double wireCentre[3];
    fGeom->WireIDToWireGeo(wireID).GetCenter(wireCentre);
    if (wireID.TPC % 2 == 0) globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 0, wireID.Cryostat);
    else globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 1, wireID.Cryostat);
  }

  // Collection
  else {
    // FOR COLLECTION WIRES, HARD CODE THE GEOMETRY FOR GIVEN DETECTORS
    // THIS _SHOULD_ BE TEMPORARY. GLOBAL WIRE SUPPORT IS BEING ADDED TO THE LARSOFT GEOMETRY AND SHOULD BE AVAILABLE SOON
    if (fDetector == "dune35t") {
      unsigned int nwires = fGeom->Nwires(wireID.Plane, 0, wireID.Cryostat);
      if (wireID.TPC == 0 or wireID.TPC == 1) globalWire = wireID.Wire;
      else if (wireID.TPC == 2 or wireID.TPC == 3 or wireID.TPC == 4 or wireID.TPC == 5) globalWire = nwires + wireID.Wire;
      else if (wireID.TPC == 6 or wireID.TPC == 7) globalWire = (2*nwires) + wireID.Wire;
      else mf::LogError("BlurredClusterAlg") << "Error when trying to find a global induction plane coordinate for TPC " << wireID.TPC << " (geometry" << fDetector << ")";
    }
    else if (fDetector == "dune10kt") {
      unsigned int nwires = fGeom->Nwires(wireID.Plane, 0, wireID.Cryostat);
      // Detector geometry has four TPCs, two on top of each other, repeated along z...
      int block = wireID.TPC / 4;
      globalWire = (nwires*block) + wireID.Wire;
    }
    else {
      double wireCentre[3];
      fGeom->WireIDToWireGeo(wireID).GetCenter(wireCentre);
      if (wireID.TPC % 2 == 0) globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 0, wireID.Cryostat);
      else globalWire = fGeom->WireCoordinate(wireCentre[1], wireCentre[2], wireID.Plane, 1, wireID.Cryostat);
    }
  }

  return globalWire;

}

TVector2 shower::TrackShowerSepAlg::Project3DPointOntoPlane(TVector3 const& point, int plane, int cryostat) {

  TVector2 wireTickPos = TVector2(-999., -999.);

  double pointPosition[3] = {point.X(), point.Y(), point.Z()};

  geo::TPCID tpcID = fGeom->FindTPCAtPosition(pointPosition);
  int tpc = 0;
  if (tpcID.isValid)
    tpc = tpcID.TPC;
  else
    return wireTickPos;

  // Construct wire ID for this point projected onto the plane
  geo::PlaneID planeID = geo::PlaneID(cryostat, tpc, plane);
  geo::WireID wireID = fGeom->NearestWireID(point, planeID);

  wireTickPos = TVector2(GlobalWire(wireID),
			 fDetProp->ConvertXToTicks(point.X(), planeID));

  // wireTickPos = TVector2(fGeom->WireCoordinate(point.Y(), point.Z(), planeID.Plane, tpc % 2, planeID.Cryostat),
  // 			 fDetProp->ConvertXToTicks(point.X(), planeID.Plane, tpc % 2, planeID.Cryostat));

  //return wireTickPos;
  return HitPosition(wireTickPos, planeID);

}








// // --------------------------- OLD (late 2015) -------------------------------


// shower::TrackShowerSepAlg::TrackShowerSepAlg(fhicl::ParameterSet const& pset) {
//   this->reconfigure(pset);
//   // tree
//   // ftree = tfs->make<TTree>("tree","tree");
//   // ftree->Branch("Run",&Run);
//   // ftree->Branch("Event",&Event);
//   // ftree->Branch("Distance",&Distance);
//   // ftree->Branch("Angle",&Angle);
//   // ftree->Branch("Length",&Length);
//   // ftree->Branch("TrackID",&TrackID);
//   // ftree->Branch("PDG",&pdg);
//   // ftree->Branch("NSpacePoints",&NSpacePoints);
//   // ftree->Branch("AvDistance",&AvDistance);
// }

// void shower::TrackShowerSepAlg::reconfigure(fhicl::ParameterSet const& pset) {
//   fConeAngle = pset.get<double>("ConeAngle");

//   // fAngleCut           = pset.get<double>("AngleCut");
//   // fDistanceCut        = pset.get<double>("DistanceCut");
//   // fVertexProximityCut = pset.get<double>("VertexProximityCut");
//   // fTrackProximityCut  = pset.get<double>("TrackProximityCut");
//   // fAvTrackHitDistance = pset.get<double>("AvTrackHitDistance");
// }

// void shower::TrackShowerSepAlg::IdentifyTracksFromEventCentre(const std::vector<art::Ptr<recob::Track> >& tracks,
// 								     const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
// 								     const art::FindManyP<recob::Track>& fmtsp) {

//   // Find the charge weighted centre of the entire event!
//   // ---- except space points don't have charge (could look at associated hits, but cba!)
//   TVector3 centre = TVector3(0,0,0);
//   for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt)
//     centre += TVector3((*spacePointIt)->XYZ()[0], (*spacePointIt)->XYZ()[1], (*spacePointIt)->XYZ()[2]);
//   centre *= 1/(double)spacePoints.size();

//   TVector3 trackVertex, trackEnd, trackDirection;

//   // Look at all tracks
//   for (std::vector<art::Ptr<recob::Track> >::const_iterator trackIt = tracks.begin(); trackIt != tracks.end(); ++trackIt) {

//     trackVertex = ( ((*trackIt)->Vertex()-centre).Mag() > ((*trackIt)->End()-centre).Mag() ) ? (*trackIt)->Vertex() : (*trackIt)->End();
//     trackEnd = ( ((*trackIt)->Vertex()-centre).Mag() > ((*trackIt)->End()-centre).Mag() ) ? (*trackIt)->End() : (*trackIt)->Vertex();
//     trackDirection = ( ((*trackIt)->Vertex()-centre).Mag() > ((*trackIt)->End()-centre).Mag() ) ? (*trackIt)->VertexDirection() : (-1)*(*trackIt)->EndDirection();

//     // Get the space points not associated with the current track
//     std::vector<art::Ptr<recob::SpacePoint> > surroundingSpacePoints = GetSurroundingSpacePoints(spacePoints, fmtsp, (*trackIt)->ID());

//     // // TRUE -- use truth to find which particle this track is associated with --------
//     // std::vector<art::Ptr<recob::Hit> > trackHits = fmht.at(trackIt->key());
//     // int trueTrackID = FindTrueTrack(trackHits);
//     // const simb::MCParticle* trueParticle = backtracker->TrackIDToParticle(trueTrackID);
//     // pdg = trueParticle->PdgCode();
//     // Length = (*trackIt)->Length();
//     // TrackID = (*trackIt)->ID();
//     // // -------------------------------------------------------------------------------

//     bool showerLike = IdentifyShowerLikeTrack(trackEnd, trackDirection, surroundingSpacePoints);

//     if (showerLike) {
//       fShowerLikeIDs.push_back((*trackIt)->ID());
//       continue;
//     }

//     else
//       fTrackLikeIDs.push_back((*trackIt)->ID());

//   }

//   return;

// }

// void shower::TrackShowerSepAlg::IdentifyTracksNearTracks(const std::vector<art::Ptr<recob::Track> >& tracks) {

//   std::vector<art::Ptr<recob::Track> > identifiedTracks;

//   for (std::vector<art::Ptr<recob::Track> >::const_iterator trackIt = tracks.begin(); trackIt != tracks.end(); ++trackIt)
//     if (std::find(fTrackLikeIDs.begin(), fTrackLikeIDs.end(), (*trackIt)->ID()) != fTrackLikeIDs.end())
//       identifiedTracks.push_back(*trackIt);

//   // Look through tracks
//   bool allTracksRemoved = false;
//   while (!allTracksRemoved) {

//     int tracksRemoved = 0;

//     // Look through all the tracks
//     for (std::vector<art::Ptr<recob::Track> >::const_iterator trackIt = tracks.begin(); trackIt != tracks.end(); ++trackIt) {

//       // Don't consider tracks already tagged as tracks!
//       if (std::find(fTrackLikeIDs.begin(), fTrackLikeIDs.end(), (*trackIt)->ID()) != fTrackLikeIDs.end())
// 	continue;

//       bool trackShouldBeRemoved = false;

//       // Find tracks which are close to previously identified tracks
//       for (std::vector<art::Ptr<recob::Track> >::iterator identifiedTrackIt = identifiedTracks.begin(); identifiedTrackIt != identifiedTracks.end(); ++identifiedTrackIt) {
// 	if (std::find(fShowerLikeIDs.begin(), fShowerLikeIDs.end(), (*trackIt)->ID()) != fShowerLikeIDs.end())
// 	  continue;
// 	if ( ( ((*trackIt)->Vertex() - (*identifiedTrackIt)->Vertex()).Mag() < fTrackProximityCut ) or
// 	     ( ((*trackIt)->Vertex() - (*identifiedTrackIt)->End()   ).Mag() < fTrackProximityCut ) or
// 	     ( ((*trackIt)->End()    - (*identifiedTrackIt)->Vertex()).Mag() < fTrackProximityCut ) or
// 	     ( ((*trackIt)->End()    - (*identifiedTrackIt)->End()   ).Mag() < fTrackProximityCut ) )
// 	  trackShouldBeRemoved = true;
//       }

//       // Tag this track as a 'track-like' object
//       if (trackShouldBeRemoved) {
// 	fTrackLikeIDs.push_back((*trackIt)->ID());
// 	identifiedTracks.push_back(*trackIt);
// 	++tracksRemoved;
//       }

//     }

//     // If there were no tracks removed then we'll call it a day
//     if (tracksRemoved == 0)
//       allTracksRemoved = true;

//   }

//   return;

// }

// void shower::TrackShowerSepAlg::IdentifyTracksNearVertex(const art::Ptr<recob::Vertex>& vertex,
// 								const std::vector<art::Ptr<recob::Track> >& tracks,
// 								const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
// 								const art::FindManyP<recob::Track>& fmtsp) {

//   double xyz[3];
//   vertex->XYZ(xyz);
//   TVector3 vertexPos = TVector3(xyz[0], xyz[1], xyz[2]);

//   // Look at each track
//   for (std::vector<art::Ptr<recob::Track> >::const_iterator trackIt = tracks.begin(); trackIt != tracks.end(); ++trackIt) {

//     TVector3 end, direction;

//     // See if either end is close to the vertex
//     if ( ((*trackIt)->Vertex() - vertexPos).Mag() < fVertexProximityCut or
// 	 ((*trackIt)->End() - vertexPos).Mag() < fVertexProximityCut) {
//       end = ((*trackIt)->Vertex() - vertexPos).Mag() < ((*trackIt)->End() - vertexPos).Mag() ? (*trackIt)->End() : (*trackIt)->VertexDirection();
//       direction = ((*trackIt)->Vertex() - vertexPos).Mag() < ((*trackIt)->End() - vertexPos).Mag() ? (*trackIt)->VertexDirection() : (-1)*(*trackIt)->VertexDirection();
//     }

//     else
//       continue;

//     // Get the space points not associated with the current track
//     std::vector<art::Ptr<recob::SpacePoint> > surroundingSpacePoints = GetSurroundingSpacePoints(spacePoints, fmtsp, (*trackIt)->ID());

//     // Make sure this track start isn't a shower start
//     bool showerLike = IdentifyShowerLikeTrack(end, direction, surroundingSpacePoints);

//     if (showerLike) {
//       fShowerLikeIDs.push_back((*trackIt)->ID());
//       continue;
//     }

//     // This track originates from near the interaction vertex and is not shower-like
//     fTrackLikeIDs.push_back((*trackIt)->ID());

//   }

//   return;

// }

// bool shower::TrackShowerSepAlg::IdentifyShowerLikeTrack(const TVector3& end,
// 							       const TVector3& direction,
// 							       const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints) {

//   std::vector<art::Ptr<recob::SpacePoint> > spacePointsInCone = GetSpacePointsInCone(spacePoints, end, direction);

//   if (spacePointsInCone.size() < 2)
//     return false;

//   // Get the average spread of these space points
//   double spread = SpacePointSpread(spacePointsInCone);

//   if (spread < fAvTrackHitDistance)
//     return false;

//   return true;

// }

// std::vector<art::Ptr<recob::Hit> > shower::TrackShowerSepAlg::FillHitsToCluster(const std::vector<art::Ptr<recob::Hit> >& initialHits,
// 										       const art::FindManyP<recob::Track>& fmt) {

//   // Container to fill with shower-like hits
//   std::vector<art::Ptr<recob::Hit> > hitsToCluster;

//   for (std::vector<art::Ptr<recob::Hit> >::const_iterator initialHit = initialHits.begin(); initialHit != initialHits.end(); ++initialHit) {
//     std::vector<art::Ptr<recob::Track> > showerTracks = fmt.at(initialHit->key());
//     if ( (showerTracks.size() and (std::find(fTrackLikeIDs.begin(), fTrackLikeIDs.end(), showerTracks.at(0)->ID()) == fTrackLikeIDs.end()))//hit on track-like track
// 	 || !showerTracks.size() )//hit not on any track
//       hitsToCluster.push_back(*initialHit);
//   }

//   return hitsToCluster;

// }

// int shower::TrackShowerSepAlg::FindTrackID(const art::Ptr<recob::Hit>& hit) {
//   double particleEnergy = 0;
//   int likelyTrackID = 0;
//   std::vector<sim::TrackIDE> trackIDs = backtracker->HitToTrackID(hit);
//   for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
//     if (trackIDs.at(idIt).energy > particleEnergy) {
//       particleEnergy = trackIDs.at(idIt).energy;
//       likelyTrackID = TMath::Abs(trackIDs.at(idIt).trackID);
//     }
//   }
//   return likelyTrackID;
// }

// int shower::TrackShowerSepAlg::FindTrueTrack(const std::vector<art::Ptr<recob::Hit> >& trackHits) {
//   std::map<int,double> trackMap;
//   for (std::vector<art::Ptr<recob::Hit> >::const_iterator trackHitIt = trackHits.begin(); trackHitIt != trackHits.end(); ++trackHitIt) {
//     art::Ptr<recob::Hit> hit = *trackHitIt;
//     int trackID = FindTrackID(hit);
//     trackMap[trackID] += hit->Integral();
//   }
//   //return std::max_element(trackMap.begin(), trackMap.end(), [](const std::pair<int,double>& p1, const std::pair<int,double>& p2) {return p1.second < p2.second;} )->first;
//   double highestCharge = 0;
//   int clusterTrack = 0;
//   for (std::map<int,double>::iterator trackIt = trackMap.begin(); trackIt != trackMap.end(); ++trackIt)
//     if (trackIt->second > highestCharge) {
//       highestCharge = trackIt->second;
//       clusterTrack  = trackIt->first;
//     }
//   return clusterTrack;
// }

// std::vector<art::Ptr<recob::SpacePoint> > shower::TrackShowerSepAlg::GetSurroundingSpacePoints(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
// 												      const art::FindManyP<recob::Track>& fmt,
// 												      unsigned int trackID) {

//   // The space points to return
//   std::vector<art::Ptr<recob::SpacePoint> > surroundingSpacePoints;

//   for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt) {

//     bool spacePointIsInCurrentTrack = false;
  
//     std::vector<art::Ptr<recob::Track> > spacePointTracks = fmt.at(spacePointIt->key());
//     for (std::vector<art::Ptr<recob::Track> >::iterator spacePointTrackIt = spacePointTracks.begin(); spacePointTrackIt != spacePointTracks.end(); ++spacePointTrackIt)
//       if (spacePointTrackIt->key() == trackID) spacePointIsInCurrentTrack = true;
    
//     if (!spacePointIsInCurrentTrack)
//       surroundingSpacePoints.push_back(*spacePointIt);

//   }

//   return surroundingSpacePoints;

// }

// std::vector<art::Ptr<recob::SpacePoint> > shower::TrackShowerSepAlg::GetSpacePointsInCone(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
// 												 const TVector3& trackEnd,
// 												 const TVector3& trackDirection) {

//   // The space points in cone to return
//   std::vector<art::Ptr<recob::SpacePoint> > spacePointsInCone;

//   TVector3 spacePointPos, spacePointProj;
//   double displacementFromAxis, distanceFromTrackEnd, projDistanceFromTrackEnd, angleFromTrackEnd;

//   for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt) {

//     // Get the properties of this space point
//     const double* xyz = (*spacePointIt)->XYZ();
//     spacePointPos = TVector3(xyz[0], xyz[1], xyz[2]);
//     spacePointProj = ((spacePointPos-trackEnd).Dot(trackDirection))*trackDirection + trackEnd;
//     displacementFromAxis = (spacePointProj - spacePointPos).Mag();
//     distanceFromTrackEnd = (spacePointPos - trackEnd).Mag();
//     projDistanceFromTrackEnd = (spacePointProj - trackEnd).Mag();
//     angleFromTrackEnd = TMath::ASin(displacementFromAxis/distanceFromTrackEnd) * 180 / TMath::Pi();

//     if ( (projDistanceFromTrackEnd < fDistanceCut) and (angleFromTrackEnd < fAngleCut) and ((spacePointProj-trackEnd).Dot(trackDirection) > 0) )
//       spacePointsInCone.push_back(*spacePointIt);

//   }

//   return spacePointsInCone;

// }

// std::vector<art::Ptr<recob::Hit> > shower::TrackShowerSepAlg::RemoveTrackHits(const std::vector<art::Ptr<recob::Hit> >& initialHits,
// 										     const std::vector<art::Ptr<recob::Track> >& tracks,
// 										     const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
// 										     const std::vector<art::Ptr<recob::Vertex> >& vertices,
// 										     const art::FindManyP<recob::Track>& fmth,
// 										     const art::FindManyP<recob::Track>& fmtsp,
// 										     const art::FindManyP<recob::Hit>& fmh,
// 										     int event,
// 										     int run) {

//   Event = event;
//   Run = run;

//   if (spacePoints.size() == 0)
//     return initialHits;

//   // Container for shower-like hits
//   std::vector<art::Ptr<recob::Hit> > hitsToCluster;

//   fTrackLikeIDs.clear();
//   fShowerLikeIDs.clear();

//   // Find the vertex furthest upstream (if it exists)
//   art::Ptr<recob::Vertex> vertex;
//   for (std::vector<art::Ptr<recob::Vertex> >::const_iterator vertexIt = vertices.begin(); vertexIt != vertices.end(); ++vertexIt) {
//     double xyzNew[3], xyzOld[3];
//     (*vertexIt)->XYZ(xyzNew);
//     if (vertex.isNull())
//       vertex = *vertexIt;
//     else {
//       vertex->XYZ(xyzOld);
//       if (xyzNew[2] < xyzOld[2])
// 	vertex = *vertexIt;
//     }
//   }

//   // If we have a vertex then, for the time being, use it to remove tracks which are too close!
//   if (!vertex.isNull())
//     this->IdentifyTracksNearVertex(vertex, tracks, spacePoints, fmtsp);

//   // If no vertex, things are harder
//   else
//     this->IdentifyTracksFromEventCentre(tracks, spacePoints, fmtsp);

//   // Once we've identified some tracks, can look for others at the ends
//   this->IdentifyTracksNearTracks(tracks);

//   hitsToCluster = FillHitsToCluster(initialHits, fmth);

//   return hitsToCluster;

// }

// std::vector<art::Ptr<recob::Hit> > shower::TrackShowerSepAlg::RemoveTrackHits(const std::vector<art::Ptr<recob::Hit> >& hits,
// 										     const std::vector<art::Ptr<recob::PFParticle> > pfParticles,
// 										     const art::FindManyP<recob::Cluster>& fmc,
// 										     const art::FindManyP<recob::Hit>& fmh) {

//   std::vector<art::Ptr<recob::Hit> > showerHits;

//   // Use information from Pandora to identify shower-like hits
//   for (std::vector<art::Ptr<recob::PFParticle> >::const_iterator pfParticleIt = pfParticles.begin(); pfParticleIt != pfParticles.end(); ++pfParticleIt) {

//     // See if this is a shower particle
//     if ((*pfParticleIt)->PdgCode() == 11) {
//       std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(pfParticleIt->key());
//       for (std::vector<art::Ptr<recob::Cluster> >::iterator clusterIt = clusters.begin(); clusterIt != clusters.end(); ++clusterIt) {
//         std::vector<art::Ptr<recob::Hit> > hits = fmh.at(clusterIt->key());
//         for (std::vector<art::Ptr<recob::Hit> >::iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt)
//           showerHits.push_back(*hitIt);
//       }
//     }

//   }

//   return showerHits;

// }

// double shower::TrackShowerSepAlg::SpacePointSpread(const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints) {

//   /// Finds the spread of a set of space points about their central axis

//   // Find the centre of the space points
//   TVector3 centre = TVector3(0,0,0);
//   for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt)
//     centre += TVector3((*spacePointIt)->XYZ()[0], (*spacePointIt)->XYZ()[1], (*spacePointIt)->XYZ()[2]);
//   centre *= 1/(double)spacePoints.size();

//   // Find the central axis of the space points
//   TPrincipal* pca = new TPrincipal(3,"");
//   for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt)
//     pca->AddRow((*spacePointIt)->XYZ());
//   pca->MakePrincipals();
//   const TMatrixD* eigenvectors = pca->GetEigenVectors();
//   TVector3 spacePointDir = TVector3((*eigenvectors)[0][0], (*eigenvectors)[1][0], (*eigenvectors)[2][0]);
//   delete pca;

//   // See if the space points form something which may resemble a track (i.e. straight line) or a shower
//   double avDistanceFromCentralAxis = 0;
//   for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt) {
//     const double* xyz = (*spacePointIt)->XYZ();
//     TVector3 spacePointPos = TVector3(xyz[0], xyz[1], xyz[2]);
//     TVector3 projectionOntoCentralAxis = ((spacePointPos-centre).Dot(spacePointDir))*spacePointDir + centre;
//     avDistanceFromCentralAxis += (spacePointPos - projectionOntoCentralAxis).Mag();
//   }
//   avDistanceFromCentralAxis /= spacePoints.size();

//   return avDistanceFromCentralAxis;

// }
