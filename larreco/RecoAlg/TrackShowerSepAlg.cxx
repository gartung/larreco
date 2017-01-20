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
  fRunTrackShowerSep = false;
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

  fRunTrackShowerSep = true;
  fShowerHits.clear();
  fTrackTracks.clear();
  fShowerStarts.clear();

  // Save info about the tracks
  std::map<int,std::unique_ptr<ReconTrack> > reconTracks;

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

    reconTracks[trackIt->key()] = std::move(track);

  }

  // Make structure to hold parameters
  TrackShowerSepParameters parameters;

  // Find properties of the tracks
  TrackProperties2D(reconTracks, parameters, hits, fmth);
  TrackProperties3D(reconTracks, parameters, spacePoints, fmtsp);

  // Look for tracks
  IdentifyTracks(reconTracks, parameters);

  // Consider removing false tracks by looking at their closest approach to any other track

  // Look at the track cone properties
  ConeProperties(reconTracks, parameters, spacePoints, fmtsp);

  // Look for showers
  IdentifyShowers(reconTracks, parameters);

  // Look for shower tracks
  IdentifyShowerTracks(reconTracks);

  // Look for shower cones
  IdentifyShowerCones(reconTracks);

  // Look at any remaining tracks
  ResolveUndeterminedTracks(reconTracks);

  if (fDebug > 0) {
    std::cout << std::endl << "Event " << event << " track shower separation:" << std::endl;
    std::cout << "Shower initial tracks are:" << std::endl;
    for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt)
      if (trackIt->second->IsShowerTrack())
	std::cout << "  " << trackIt->first << std::endl;
    std::cout << "Track tracks are:" << std::endl;
    for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt)
      if (trackIt->second->IsTrack())
	std::cout << "  " << trackIt->first << "\t\t\t\t\t"
		  << "Start (" << trackIt->second->Vertex().X() << ", " << trackIt->second->Vertex().Y() << ", " << trackIt->second->Vertex().Z() << "), "
		  << "end (" << trackIt->second->End().X() << ", " << trackIt->second->End().Y() << ", " << trackIt->second->End().Z() << "), "
		  << "length " << trackIt->second->Length() << std::endl;
  }

  // Select all hits which aren't associated with a determined track
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    bool showerHit = true;
    const std::vector<art::Ptr<recob::Track> > hitTracks = fmth.at(hitIt->key());
    for (std::vector<art::Ptr<recob::Track> >::const_iterator hitTrackIt = hitTracks.begin(); hitTrackIt != hitTracks.end(); ++hitTrackIt)
      if (reconTracks[hitTrackIt->key()]->IsTrack())
	showerHit = false;
    if (showerHit)
      fShowerHits.push_back(*hitIt);
  }

  // Select all tracks which were determined to be tracks
  for (std::vector<art::Ptr<recob::Track> >::const_iterator trackIt = tracks.begin(); trackIt != tracks.end(); ++trackIt)
    if (reconTracks[trackIt->key()]->IsTrack())
      fTrackTracks.push_back(*trackIt);

  // Save shower starts
  for (std::map<int,std::unique_ptr<ReconTrack> >::const_iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt)
    if (trackIt->second->IsShowerTrack())
      fShowerStarts.push_back(trackIt->second->Vertex());

}

void shower::TrackShowerSepAlg::TrackProperties2D(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks,
						  TrackShowerSepParameters& parameters,
						  const std::vector<art::Ptr<recob::Hit> >& hits,
						  const art::FindManyP<recob::Track> fmth) {

  // Look at hits in the rectangle around each tracks
  double avRectangleHits = 0;
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt) {
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
  avRectangleHits /= (double)reconTracks.size();

  if (fDebug > 2) {
    std::cout << std::endl << "Rectangle hit ratios:" << std::endl;
    for (std::map<int,std::unique_ptr<ReconTrack> >::const_iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt) {
      std::cout << "  Track " << trackIt->first << " has " << trackIt->second->NumRectangleHits() << ", ratio " << trackIt->second->RectangleHitRatio()
		<< " (event average ratio " << avRectangleHits << ") -- "
		<< trackIt->second->RectangleHitRatio()/(double)avRectangleHits << " (cut " << fRectangleCut << ")" << std::endl;
      // const std::map<int,std::vector<art::Ptr<recob::Hit> > >& trackHits = trackIt->second->Hits();
      // for (std::map<int,std::vector<art::Ptr<recob::Hit> > >::const_iterator hitIt = trackHits.begin(); hitIt != trackHits.end(); ++hitIt)
      //   std::cout << "    Plane " << hitIt->first << " has " << trackIt->second->NumRectangleHits(hitIt->first)
      // 		  << " (ratio " << trackIt->second->RectangleHitRatio(hitIt->first) << ")" << std::endl;
    }
  }

  parameters.SetAvRectangleHits(avRectangleHits);

  return;

}

void shower::TrackShowerSepAlg::TrackProperties3D(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks,
						  TrackShowerSepParameters& parameters,
						  const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
						  const art::FindManyP<recob::Track>& fmtsp) {

  // Consider the space point cylinder situation
  double avCylinderSpacePoints = 0;
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt) {

    // Get the 3D properties of the track
    TVector3 point = trackIt->second->Vertex();
    TVector3 direction = trackIt->second->Direction3D();

    // Count space points in the volume around the track
    for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt) {

      // Look at all space points which are not associated with the track
      const std::vector<art::Ptr<recob::Track> > spTracks = fmtsp.at(spacePointIt->key());
      if (find_if(spTracks.begin(), spTracks.end(), [&trackIt](const art::Ptr<recob::Track>& t){ return (int)t.key() == trackIt->first; }) != spTracks.end())
	continue;

      // Get the properties of this space point
      TVector3 pos = SpacePointPos(*spacePointIt);
      TVector3 proj = ProjPoint(pos, direction, point);
      if ((pos-proj).Mag() < TMath::Min(fCylinderRadius,trackIt->second->Length()))
	trackIt->second->AddCylinderSpacePoint(spacePointIt->key());
      // if ((pos-proj).Mag() < fCylinderRadius and
      // 	  (pos-point)*direction > 0 and
      // 	  (pos-point)*direction < trackIt->second->Length())
      // 	std::cout << "Space point " << spacePointIt->key() << " (" << pos.X() << ", " << pos.Y() << ", " << pos.Z()
      // 		  << ") in cylinder around track " << trackIt->first << " (assocatied with track " << spTracks.at(0).key()
      // 		  << "); point is (" << point.X() << ", " << point.Y() << ", " << point.Z() << "), proj end ("
      // 		  << ((trackIt->second->Length()*trackIt->second->VertexDirection())+point).X() << ", "
      // 		  << ((trackIt->second->Length()*trackIt->second->VertexDirection())+point).Y() << ", "
      // 		  << ((trackIt->second->Length()*trackIt->second->VertexDirection())+point).Z() << ")" << std::endl;
    }
    avCylinderSpacePoints += trackIt->second->CylinderSpacePointRatio();
  }
  avCylinderSpacePoints /= (double)reconTracks.size();

  if (fDebug > 2) {
    std::cout << std::endl << "Cylinder space point ratios:" << std::endl;
    for (std::map<int,std::unique_ptr<ReconTrack> >::const_iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt)
      std::cout << "  Track " << trackIt->first << " has cylinder space point ratio " << trackIt->second->CylinderSpacePointRatio()
		<< " (event average " << avCylinderSpacePoints << ") -- "
		<< trackIt->second->CylinderSpacePointRatio()/(double)avCylinderSpacePoints << " (cut " << fCylinderCut << ")" << std::endl;
  }

  double avCylinderRectangle = 0;
  for (std::map<int,std::unique_ptr<ReconTrack> >::const_iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt)
    avCylinderRectangle += trackIt->second->CylinderSpacePointRatio()*trackIt->second->RectangleHitRatio();
  avCylinderRectangle /= (double)reconTracks.size();

  if (fDebug > 2) {
    std::cout << std::endl << "Combinded cylinder/rectangle:" << std::endl;
    for (std::map<int,std::unique_ptr<ReconTrack> >::const_iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt)
      std::cout << "  Track " << trackIt->first << ": " << (trackIt->second->CylinderSpacePointRatio()) * (trackIt->second->RectangleHitRatio())
		<< "\tEvent average is " << avCylinderRectangle << std::endl;
  }

  if (fDebug > 2)
    std::cout << std::endl << "Looking at track straightness:" << std::endl;
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt) {
    const std::vector<art::Ptr<recob::SpacePoint> >& trackPoints = trackIt->second->SpacePoints();
    double least_sq = 0;
    for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = trackPoints.begin(); spacePointIt != trackPoints.end(); ++spacePointIt) {
      TVector3 pos = SpacePointPos(*spacePointIt);
      TVector3 proj = ProjPoint(pos, trackIt->second->Direction3D(), trackIt->second->Centre3D());
      least_sq += TMath::Power((proj-pos).Mag(),2);
    }
    least_sq /= (double)((int)trackPoints.size()-2);
    trackIt->second->SetLeastSquareNDOF(least_sq);
    if (fDebug > 2)
      std::cout << "  Track " << trackIt->first << " has least squared / dof " << least_sq << std::endl;
  }

  parameters.SetAvCylinderSpacePoints(avCylinderSpacePoints);

  return;

}

void shower::TrackShowerSepAlg::IdentifyTracks(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks, const TrackShowerSepParameters& parameters) {

  // These methods have been developed over time and work pretty well. I've identified some variables that work well at discriminating tracks from shower cones
  // Really, they need to be implemented in the form of a likelihood
  // Also, rather than boolean TrackLike, VeryTrackLike etc, there should be a variable on a continuous scale that determines the 'track-iness' of tracks
  // (Notes mainly to myself if I ever revisit this)
  // MW, 1/18/2017

  if (fDebug > 1)
    std::cout << std::endl << "Identifying tracks:" << std::endl;

  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt) {
    double trackCylinder = trackIt->second->CylinderSpacePointRatio() / parameters.AverageCylinderSpacePoints;
    double trackRectangle = trackIt->second->RectangleHitRatio() / parameters.AverageRectangleHits;
    double trackLeastSq = trackIt->second->LeastSquareNDOF();
    if (trackCylinder < fCylinderCut) {
      trackIt->second->IncreaseTrackiness();
      if (trackCylinder < fCylinderCut / 5.)
	trackIt->second->IncreaseTrackiness();
      if (trackCylinder < fCylinderCut / 10.)
	trackIt->second->IncreaseTrackiness();
    }
    if (trackRectangle < fRectangleCut) {
      trackIt->second->IncreaseTrackiness();
      if (trackRectangle < fRectangleCut / 5.)
	trackIt->second->IncreaseTrackiness();
      if (trackRectangle < fRectangleCut / 10.)
	trackIt->second->IncreaseTrackiness();
    }
    if (trackLeastSq < fLeastSquareCut) {
      trackIt->second->IncreaseTrackiness();
      if (trackLeastSq < fLeastSquareCut / 10.)
	trackIt->second->IncreaseTrackiness();
    }

    if (trackIt->second->Trackiness() > 4) {
      if (fDebug > 1)
	std::cout << "  Making track " << trackIt->first << " a track (Type I)"
		  << " -- Cylinder " << trackCylinder << " (cut " << fCylinderCut << ")"
		  << " Rectangle " << trackRectangle << " (cut " << fRectangleCut << ")"
		  << " Least square " << trackLeastSq << " (cut " << fLeastSquareCut << ")" << std::endl;
      trackIt->second->MakeTrack();
    }

  }

  // for (std::map<int,std::unique_ptr<ReconTrack> >::const_iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt) {
  //   if (!trackIt->second->IsTrack())
  //     continue;
  //   std::cout << "Track " << trackIt->first << " (tagged as track) is this close to other tracks:" << std::endl;
  //   for (std::map<int,std::unique_ptr<ReconTrack> >::const_iterator otherTrackIt = reconTracks.begin(); otherTrackIt != reconTracks.end(); ++otherTrackIt) {
  //     if (trackIt->first == otherTrackIt->first)
  // 	continue;
  //     std::cout << "  Other track " << otherTrackIt->first << " from track " << trackIt->first << " vertex: " << (otherTrackIt->second->Vertex()-trackIt->second->Vertex()).Mag() << " (vertex) and " << (otherTrackIt->second->End()-trackIt->second->Vertex()).Mag() << " (end); end: " << (otherTrackIt->second->Vertex()-trackIt->second->End()).Mag() << " (vertex) and " << (otherTrackIt->second->End()-trackIt->second->End()).Mag() << " (end)" << std::endl;
  //   }
  // }

  // // Identify further tracks
  // for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt) {
  //   if (trackIt->second->IsTrack())
  //     continue;
  //   for (std::map<int,std::unique_ptr<ReconTrack> >::const_iterator otherTrackIt = reconTracks.begin(); otherTrackIt != reconTracks.end(); ++otherTrackIt) {
  //     if (trackIt->first == otherTrackIt->first or !otherTrackIt->second->IsTrack())
  // 	continue;
  //     if (fDebug > 3) {
  // 	std::cout << "Distances from tracks to tracks:" << std::endl;
  // 	double minDist = TMath::Min((trackIt->second->Vertex()-otherTrackIt->second->Vertex()).Mag(),
  // 				    TMath::Min((trackIt->second->Vertex()-otherTrackIt->second->End()).Mag(),
  // 					       TMath::Min((trackIt->second->End()-otherTrackIt->second->Vertex()).Mag(),
  // 							  (trackIt->second->End()-otherTrackIt->second->End()).Mag())));

  // 	std::cout << "  Min distance between track " << trackIt->first << " and " << otherTrackIt->first << " is " << minDist << std::endl;
  //     }
  //     if ((trackIt->second->Vertex() - otherTrackIt->second->Vertex()).Mag() < fTrackVertexCut or
  // 	  (trackIt->second->Vertex() - otherTrackIt->second->End()).Mag() < fTrackVertexCut or
  // 	  (trackIt->second->End() - otherTrackIt->second->Vertex()).Mag() < fTrackVertexCut or
  // 	  (trackIt->second->End() - otherTrackIt->second->End()).Mag() < fTrackVertexCut) {
  // 	if (fDebug > 1)
  // 	  std::cout << "  Making track " << trackIt->first << " a track (Type II) (close to " << otherTrackIt->first << ")" << std::endl;
  // 	trackIt->second->MakeTrack();
  // 	otherTrackIt->second->AddTrackConversion(trackIt->first);
  //     }
  //   }
  // }

  return;

}

void shower::TrackShowerSepAlg::ConeProperties(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks,
					       TrackShowerSepParameters& parameters,
					       const std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
					       const art::FindManyP<recob::Track>& fmtsp) {

  // Identify tracks which slipped through and shower tracks
  // For the moment, until the track tagging gets better at least, don't try to identify tracks from this

  double avConeSize = 0;
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt) {

    // Do everything via space points
    for (std::vector<art::Ptr<recob::SpacePoint> >::const_iterator spacePointIt = spacePoints.begin(); spacePointIt != spacePoints.end(); ++spacePointIt) {

      // Ignore space points associated with tracks
      bool associatedSpacePoint = false, trackLikeSpacePoint = false;
      const std::vector<art::Ptr<recob::Track> > spTracks = fmtsp.at(spacePointIt->key());
      for (std::vector<art::Ptr<recob::Track> >::const_iterator spTrackIt = spTracks.begin(); spTrackIt != spTracks.end(); ++spTrackIt) {
	if ((int)spTrackIt->key() == trackIt->first)
	  associatedSpacePoint = true;
	if (reconTracks[(int)spTrackIt->key()]->Trackiness() > 6)
	  trackLikeSpacePoint = true;
      }
      if (associatedSpacePoint)
	continue;

      // Add space points and tracks in both the forward and backward cones
      if ((SpacePointPos(*spacePointIt) - trackIt->second->Vertex()).Angle(trackIt->second->Direction3D()) < fConeAngle * TMath::Pi() / 180) {
	trackIt->second->AddForwardConeSpacePoint(spacePointIt->key());
	trackIt->second->AddForwardConeTrack(spTracks.at(0).key());
	if (!trackLikeSpacePoint) {
	  trackIt->second->AddForwardShowerConeSpacePoint(spacePointIt->key());
	  trackIt->second->AddForwardShowerConeTrack(spTracks.at(0).key());
	}
      }
      if ((SpacePointPos(*spacePointIt) - trackIt->second->Vertex()).Angle(-1*trackIt->second->Direction3D()) < fConeAngle * TMath::Pi() / 180) {
	trackIt->second->AddBackwardConeSpacePoint(spacePointIt->key());
	trackIt->second->AddBackwardConeTrack(spTracks.at(0).key());
	if (!trackLikeSpacePoint) {
	  trackIt->second->AddBackwardShowerConeSpacePoint(spacePointIt->key());
	  trackIt->second->AddBackwardShowerConeTrack(spTracks.at(0).key());
	}
      }
      if ((SpacePointPos(*spacePointIt) - trackIt->second->End()).Angle(-1*trackIt->second->Direction3D()) < fConeAngle * TMath::Pi() / 180) {
	trackIt->second->AddForwardConeSpacePointEnd(spacePointIt->key());
	trackIt->second->AddForwardConeTrackEnd(spTracks.at(0).key());
	if (!trackLikeSpacePoint) {
	  trackIt->second->AddForwardShowerConeSpacePointEnd(spacePointIt->key());
	  trackIt->second->AddForwardShowerConeTrackEnd(spTracks.at(0).key());
	}
      }
      if ((SpacePointPos(*spacePointIt) - trackIt->second->End()).Angle(trackIt->second->Direction3D()) < fConeAngle * TMath::Pi() / 180) {
	trackIt->second->AddBackwardConeSpacePointEnd(spacePointIt->key());
	trackIt->second->AddBackwardConeTrackEnd(spTracks.at(0).key());
	if (!trackLikeSpacePoint) {
	  trackIt->second->AddBackwardShowerConeSpacePointEnd(spacePointIt->key());
	  trackIt->second->AddBackwardShowerConeTrackEnd(spTracks.at(0).key());
	}
      }
    }
    avConeSize += trackIt->second->ShowerConeSpacePointSize();
  }
  avConeSize /= (double)reconTracks.size();

  parameters.SetAvConeSize(avConeSize);

  return;

}

void shower::TrackShowerSepAlg::IdentifyShowers(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks, const TrackShowerSepParameters& parameters) {

  if (fDebug > 1)
    std::cout << std::endl << "Identifying showers:" << std::endl;

  // Find tracks which have large shower cones
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt) {

    // if (trackIt->second->NumForwardShowerConeSpacePoints() == 0)
    //   trackIt->second->MakeTrack();
    // double distanceFromAverage = (trackIt->second->ConeSize() - avConeSize) / TMath::Abs(avConeSize);
    // if (fDebug > 2)
    //   std::cout << "    Track " << trackIt->first << " has cone size " << trackIt->second->ConeSize() << " (forward " << trackIt->second->ForwardSpacePoints() << ") and tracks " << trackIt->second->TrackConeSize() << " (average " << avConeSize << ", distance from average " << distanceFromAverage << ")" << std::endl;
    // if (distanceFromAverage > fShowerConeCut) {
    //   trackIt->second->MakeShower();
    //   if (fDebug > 1)
    // 	std::cout << "  Making track " << trackIt->first << " a shower (Type I)" << std::endl;
    // }

    if (fDebug > 2)
      std::cout << "  Track " << trackIt->first << " has space point shower cone size " << trackIt->second->ShowerConeSpacePointSize()
		<< " (average " << parameters.AverageConeSize << ") and track shower cone size " << trackIt->second->ShowerConeTrackSize()
		<< " (trackiness " << trackIt->second->Trackiness() << ")" << std::endl;

    // Large shower cones...
    if (TMath::Abs(trackIt->second->ShowerConeSpacePointSize()) > 30 and //parameters.AverageConeSize and
	TMath::Abs(trackIt->second->ShowerConeTrackSize()) > 2 and
	trackIt->second->Trackiness() < 8) {
    //if (trackIt->second->ShowerConeSpacePointSize() > 30 and trackIt->second->ShowerConeTrackSize() > 3 and !trackIt->second->IsVeryTrackLike()) {
      if (trackIt->second->IsTrackLike() and 
      	  (trackIt->second->ShowerConeSpacePointSize() < 0 or trackIt->second->ShowerConeSpacePointSizeEnd() > trackIt->second->ShowerConeSpacePointSize())) {
      	if (fDebug > 1)
      	  std::cout << "      Flipping track: shower cone size "
      		    << trackIt->second->ShowerConeSpacePointSize() << "; from end " << trackIt->second->ShowerConeSpacePointSizeEnd() << std::endl;
      	trackIt->second->FlipTrack();
      }
      if (trackIt->second->ShowerConeSpacePointSize() > parameters.AverageConeSize * 5)
	trackIt->second->MakeTrack();
      trackIt->second->MakeShower();
      if (fDebug > 1)
	std::cout << "    Making track " << trackIt->first << " a shower (Type I)" << std::endl;
      const std::vector<int>& forwardTracks = trackIt->second->ForwardConeTracks();
      for (std::vector<int>::const_iterator forwardTrackIt = forwardTracks.begin(); forwardTrackIt != forwardTracks.end(); ++forwardTrackIt)
	reconTracks[*forwardTrackIt]->MakeShowerLike();
    }
  }

  return;

}

void shower::TrackShowerSepAlg::IdentifyShowerTracks(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks) {

  if (fDebug > 1)
    std::cout << std::endl << "  Looking for competing shower tracks" << std::endl;

  // Look at the shower tracks and make sure there is only one per shower

  // First, identify 'competing shower tracks' which have the same tracks in their cones
  bool presentShowerTrack = false, presentShower = false;
  std::vector<std::pair<int,int> > competingShowerTracks;
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt) {
    if (trackIt->second->IsShower())
      presentShower = true;
    if (!trackIt->second->IsShowerTrack())
      continue;
    presentShowerTrack = true;
    const std::vector<int>& forwardTracks = trackIt->second->ForwardConeTracks();
    for (std::map<int,std::unique_ptr<ReconTrack> >::iterator otherTrackIt = reconTracks.begin(); otherTrackIt != reconTracks.end(); ++otherTrackIt) {
      if (trackIt->first == otherTrackIt->first or !otherTrackIt->second->IsShowerTrack())
	continue;
      const std::vector<int>& otherForwardTracks = otherTrackIt->second->ForwardConeTracks();
      for (std::vector<int>::const_iterator forwardTrackIt = forwardTracks.begin(); forwardTrackIt != forwardTracks.end(); ++forwardTrackIt)
	if (std::find(otherForwardTracks.begin(), otherForwardTracks.end(), *forwardTrackIt) != otherForwardTracks.end()) {
	  competingShowerTracks.push_back(std::make_pair(trackIt->first, otherTrackIt->first));
	  if (fDebug > 1)
	    std::cout << "    Adding competing shower track pair: " << trackIt->first << " and " << otherTrackIt->first << std::endl;
	}
    }
  }

  // If no shower track, make track with the largest shower cone the shower track
  if (presentShower and !presentShowerTrack) {
    int largestShowerTrack = -1, showerTrackSize = 0;
    for (std::map<int,std::unique_ptr<ReconTrack> >::const_iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt) {
      if (trackIt->second->ShowerConeSpacePointSize() > showerTrackSize) {
	largestShowerTrack = trackIt->first;
	showerTrackSize = trackIt->second->ShowerConeSpacePointSize();
      }
    }
    reconTracks[largestShowerTrack]->MakeShowerTrack();
    if (fDebug > 1)
      std::cout << std::endl << "  Making " << largestShowerTrack << " a shower track (there were previously none)"<< std::endl;
  }

  // Now look at each of these pairs of competing shower tracks and decide which should be kept
  for (std::vector<std::pair<int,int> >::const_iterator trackPairIt = competingShowerTracks.begin(); trackPairIt != competingShowerTracks.end(); ++trackPairIt) {
    if (!reconTracks.at(trackPairIt->first)->IsShowerTrack() or !reconTracks.at(trackPairIt->second)->IsShowerTrack())
      continue;
    int nonShowerTrack;
    const std::unique_ptr<ReconTrack>& firstTrack = reconTracks.at(trackPairIt->first);
    const std::unique_ptr<ReconTrack>& secondTrack = reconTracks.at(trackPairIt->second);
    if ((firstTrack->Length() < 5) xor (secondTrack->Length() < 5))
      nonShowerTrack = firstTrack->Length() < secondTrack->Length() ? trackPairIt->first : trackPairIt->second;
    else if (TMath::Abs(firstTrack->Vertex().Z() - secondTrack->Vertex().Z()) > 50)
      nonShowerTrack = firstTrack->Vertex().Z() - secondTrack->Vertex().Z() > 0 ? trackPairIt->first : trackPairIt->second;
    else if ((firstTrack->ShowerConeSpacePointSize() - secondTrack->ShowerConeSpacePointSize()) / (double)firstTrack->ShowerConeSpacePointSize() < 0.05 and
	     (firstTrack->Length() / (double)secondTrack->Length() > 10 or secondTrack->Length() / (double)firstTrack->Length() > 10))
      nonShowerTrack = firstTrack->Length () < secondTrack->Length() ? trackPairIt->first : trackPairIt->second;
    else
      nonShowerTrack = firstTrack->ShowerConeSpacePointSize() > secondTrack->ShowerConeSpacePointSize() ? trackPairIt->second : trackPairIt->first;
    if (fDebug > 1)
      std::cout << "    Downgrading " << nonShowerTrack << " to non-shower track (comparison between " << trackPairIt->first << " and "
		<< trackPairIt->second << ")" << std::endl;
    DowngradeTrack(reconTracks, nonShowerTrack);
  }

  // Make sure there aren't any tracks left which were made track by shower tracks and are in the cone
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt) {
    if (!trackIt->second->IsShowerTrack())
      continue;
    const std::vector<int>& trackConversions = trackIt->second->TrackConversions();
    const std::vector<int>& forwardTracks = trackIt->second->ForwardConeTracks();
    for (std::vector<int>::const_iterator trackConversionIt = trackConversions.begin(); trackConversionIt != trackConversions.end(); ++trackConversionIt)
      if (std::find(forwardTracks.begin(), forwardTracks.end(), *trackConversionIt) != forwardTracks.end())
	DowngradeTrack(reconTracks, *trackConversionIt);
  }

  return;

}

void shower::TrackShowerSepAlg::DowngradeTrack(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks, int track) {

  // Downgrade the identified track
  reconTracks[track]->MakeUndetermined();

  // Make all tracks which used this track to become tracks undetermined
  const std::vector<int>& trackConversions = reconTracks.at(track)->TrackConversions();
  for (std::vector<int>::const_iterator trackConversionIt = trackConversions.begin(); trackConversionIt != trackConversions.end(); ++trackConversionIt) {
    if (fDebug > 1)
    std::cout << "      Downgrading " << *trackConversionIt << " (was made a track by " << track << ", which is no longer" << std::endl;
    DowngradeTrack(reconTracks, *trackConversionIt);
  }

  return;

}

void shower::TrackShowerSepAlg::IdentifyShowerCones(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks) {

  if (fDebug > 2)
    std::cout << std::endl << "  Shower cones:" << std::endl;

  // Look through all tracks which are determined to be part of showers
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt) {
    if (!trackIt->second->IsShower())
      continue;

    if (fDebug > 2)
      std::cout << "    Track " << trackIt->first << std::endl;

    // Look through all other tracks
    for (std::map<int,std::unique_ptr<ReconTrack> >::iterator otherTrackIt = reconTracks.begin(); otherTrackIt != reconTracks.end(); ++otherTrackIt) {

      // Only consider tracks which haven't been determined
      if (trackIt->first == otherTrackIt->first or
	  !otherTrackIt->second->IsUndetermined() or
	  otherTrackIt->second->Length() > 100 or
	  otherTrackIt->second->Trackiness() > 6 or
	  (otherTrackIt->second->Vertex() - trackIt->second->Vertex()).Mag() < 5 or
	  (otherTrackIt->second->IsTrackLike() and !otherTrackIt->second->IsShowerLike()))
	continue;

      if (fDebug > 3)
	std::cout << "      Other track " << otherTrackIt->first << " has angle from vertex "
		  << (otherTrackIt->second->Vertex()-trackIt->second->Vertex()).Angle(trackIt->second->Direction3D()) * 180/TMath::Pi() << " and angle from end "
		  << (otherTrackIt->second->End()-trackIt->second->Vertex()).Angle(trackIt->second->Direction3D()) * 180/TMath::Pi() << std::endl;

      // Make all tracks in the shower cone shower-like tracks
      if ((otherTrackIt->second->Vertex()-trackIt->second->Vertex()).Angle(trackIt->second->Direction3D()) < fConeAngle * TMath::Pi() / 180 or
	  (otherTrackIt->second->End()-trackIt->second->Vertex()).Angle(trackIt->second->Direction3D()) < fConeAngle * TMath::Pi() / 180) {

	if (fDebug > 2)
	  std::cout << "      " << otherTrackIt->first << std::endl;
	otherTrackIt->second->MakeShower();
	otherTrackIt->second->AddShowerTrack(trackIt->first);
	trackIt->second->AddShowerCone(otherTrackIt->first);
	if (fDebug > 1)
	  std::cout << "      Making track " << otherTrackIt->first << " a shower (Type II)" << std::endl;

      }
    }
  }

  return;

}

void shower::TrackShowerSepAlg::ResolveUndeterminedTracks(std::map<int,std::unique_ptr<ReconTrack> >& reconTracks) {

  // Look at remaining undetermined tracks
  for (std::map<int,std::unique_ptr<ReconTrack> >::iterator trackIt = reconTracks.begin(); trackIt != reconTracks.end(); ++trackIt) {

    if (!trackIt->second->IsUndetermined())
      continue;

    if (trackIt->second->Trackiness() >= 3)
      trackIt->second->MakeTrack();

    else
      trackIt->second->MakeShower();

  }

  return;

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

  if (!fRunTrackShowerSep)
    std::cout << "Warning! (TrackShowerSepAlg): no shower hits found. (Did you RunTrackShowerSep first?)" << std::endl;

  return fShowerHits;

}

std::vector<TVector3> shower::TrackShowerSepAlg::ShowerStarts() {

  if (!fRunTrackShowerSep)
    std::cout << "Warning! (TrackShowerSepAlg): no shower starts found. (Did you RunTrackShowerSep first?)" << std::endl;

  return fShowerStarts;  

}

std::vector<art::Ptr<recob::Track> > shower::TrackShowerSepAlg::TrackTracks() {

  if (!fRunTrackShowerSep)
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
