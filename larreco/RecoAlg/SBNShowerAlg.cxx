#include "larreco/RecoAlg/SBNShowerAlg.h"

shower::SBNShowerAlg::SBNShowerAlg(const fhicl::ParameterSet& pset):  
  fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())
  //,tfs(art::ServiceHandle<art::TFileService>());
{
  fUseCollectionOnly     = pset.get<bool>("UseCollectionOnly");
  fPFParticleModuleLabel = pset.get<art::InputTag> ("PFParticleModuleLabel");
  fHitModuleLabel        = pset.get<art::InputTag> ("HitModuleLabel");
}

void shower::SBNShowerAlg::OrderShowerHits(std::vector<art::Ptr<recob::Hit> >& hits, 
					   TVector3& ShowerStartPosition,
					   TVector3& ShowerDirection 
					   ){
  
  std::map<double, art::Ptr<recob::Hit> > OrderedHits;
  art::Ptr<recob::Hit> startHit = hits.front();

  //Get the wireID 
  const geo::WireID startWireID = startHit->WireID();
  
  //Get the plane 
  const geo::PlaneID planeid = startWireID.asPlaneID();

  //Get the pitch 
  double pitch = fGeom->WirePitch(planeid);

  TVector2 Shower2DStartPosition = { 
    fGeom->WireCoordinate(ShowerStartPosition, startHit->WireID().planeID())*pitch, 
    ShowerStartPosition.X()
  };
  
  //Vector of the plane
  TVector3 PlaneDirection = fGeom->Plane(planeid).GetIncreasingWireDirection();

  //get the shower 2D direction 
  TVector2 Shower2DDirection = { 
    ShowerDirection.Dot(PlaneDirection),
    ShowerDirection.X()
  };
    
  
  Shower2DDirection = Shower2DDirection.Unit();

  for(auto const& hit: hits){ 
    
    //Get the wireID
    const geo::WireID WireID = hit->WireID();

    if (WireID.asPlaneID() != startWireID.asPlaneID()) {
      break;
    }
    
    //Get the hit Vector.
    TVector2 hitcoord = HitCoordinates(hit);
    
    //Order the hits based on the projection
    TVector2 pos = hitcoord - Shower2DStartPosition;
    OrderedHits[pos*Shower2DDirection] = hit; 
  }

  //Transform the shower. 
  std::vector<art::Ptr<recob::Hit> > showerHits;
  std::transform(OrderedHits.begin(), OrderedHits.end(), std::back_inserter(showerHits), [](std::pair<double,art::Ptr<recob::Hit> > const& hit) { return hit.second; });

  //Sometimes get the order wrong. Depends on direction compared to the plane Correct for it here:
  art::Ptr<recob::Hit> frontHit = showerHits.front();
  art::Ptr<recob::Hit> backHit  = showerHits.back();
  
  //Get the hit Vector.  
  TVector2 fronthitcoord = HitCoordinates(frontHit);
  TVector2 frontpos = fronthitcoord - Shower2DStartPosition;

  
  //Get the hit Vector.                                                                            
  TVector2 backhitcoord  = HitCoordinates(backHit);
  TVector2 backpos = backhitcoord - Shower2DStartPosition;

  double frontproj = frontpos*Shower2DDirection;
  double backproj  = backpos*Shower2DDirection;
  if (TMath::Abs(backproj) < TMath::Abs(frontproj)){
   std::reverse(showerHits.begin(),showerHits.end());   
  }

  hits = showerHits;
  return;
}

void shower::SBNShowerAlg::OrderShowerSpacePoints( std::vector<art::Ptr<recob::SpacePoint> >& 
						   showersps, TVector3& vertex, 
						   TVector3& direction){

  std::map<double,art::Ptr<recob::SpacePoint> > OrderedSpacePoints;

  //Loop over the spacepoints and get the pojected distance from the vertex.                       
  for(auto const& sp: showersps){

    // Get the projection of the space point along the direction 
    double len = SpacePointProjection(sp, vertex, direction);

    //Add to the list                                                                              
    OrderedSpacePoints[len] = sp;
  }

  //Return an ordered list.                                                                        
  showersps.clear();
  for(auto const& sp: OrderedSpacePoints){
    showersps.push_back(sp.second);
  }
  return;
}

TVector3 shower::SBNShowerAlg::ShowerCentre(std::vector<art::Ptr<recob::SpacePoint> >& 
					    showerspcs, art::FindManyP<recob::Hit>& fmh){
  float totalCharge=0;
  TVector3 centre =  shower::SBNShowerAlg::ShowerCentre(showerspcs,fmh,totalCharge);
  return centre;

}


TVector3 shower::SBNShowerAlg::ShowerCentre(std::vector<art::Ptr<recob::SpacePoint> >& showersps,
					    art::FindManyP<recob::Hit>& fmh, float& totalCharge){

  TVector3 pos, chargePoint = TVector3(0,0,0);

  //Loop over the spacepoints and get the charge weighted center.
  for(auto const& sp: showersps){

    //Get the position of the spacepoint
    pos = SpacePointPosition(sp);

    //Get the associated hits
    std::vector<art::Ptr<recob::Hit> > hits = fmh.at(sp.key());

    //Average the charge unless sepcified.
    float charge  = 0;
    float charge2 = 0;
    for(auto const& hit: hits){

      if(fUseCollectionOnly){
	if(hit->SignalType() == geo::kCollection){
	  charge = hit->Integral();
	  break;
	}
      }

      //Check if any of the points are not withing 2 sigma.
      if(!fUseCollectionOnly){
	charge += hit->Integral();
	charge2 += hit->Integral() * hit->Integral();
      }
    }

    if(!fUseCollectionOnly){
      //Calculate the unbiased standard deviation and mean.
      float mean = charge/((float) hits.size());

      float rms = 1;
      
      if(hits.size() > 1){
	rms  = TMath::Sqrt((charge2 - charge*charge)/((float)(hits.size()-1)));
      }

      charge = 0;
      for(auto const& hit: hits){
	if(hit->Integral() > (mean - 2*rms) && hit->Integral() < (mean + 2*rms))
	  charge += hit->Integral();
      }
    }

    chargePoint += charge * pos;
    totalCharge += charge;

    if(charge == 0){
      mf::LogWarning("ShowerStartPosition") << 
	"Averaged charge, within 2 sigma, for a spacepoint is zero, Maybe this not a good method";
    }
  }
    
  double intotalcharge = 1/totalCharge;
  TVector3 centre = chargePoint *  intotalcharge;
  return centre;
 
}


TVector3 shower::SBNShowerAlg::SpacePointPosition(const art::Ptr<recob::SpacePoint>& sp){

  const Double32_t* sp_xyz = sp->XYZ();
  TVector3 sp_postiion = {sp_xyz[0], sp_xyz[1], sp_xyz[2]};
  return sp_postiion;
}

double shower::SBNShowerAlg::SpacePointCharge(art::Ptr<recob::SpacePoint> sp,
					      art::FindManyP<recob::Hit>& fmh){

  double Charge = 0;

  //Average over the charge even though there is only one 
  std::vector<art::Ptr<recob::Hit> > hits = fmh.at(sp.key());
  for(auto const& hit: hits){
    Charge += hit->Integral();
  }

  Charge /= hits.size();

  return Charge;
}


TVector2 shower::SBNShowerAlg::HitCoordinates(art::Ptr<recob::Hit> const& hit) {

  //Get the pitch
  const geo::WireID  WireID = hit->WireID();
  const geo::PlaneID planeid = WireID.asPlaneID();
  double pitch = fGeom->WirePitch(planeid);
  
  return TVector2(WireID.Wire*pitch ,fDetProp->ConvertTicksToX(hit->PeakTime(),planeid));
}

double shower::SBNShowerAlg::SpacePointProjection(const art::Ptr<recob::SpacePoint>&sp,
						 TVector3& vertex, TVector3& direction){
  
  // Get the position of the spacepoint                                                           
  TVector3 pos = shower::SBNShowerAlg::SpacePointPosition(sp) - vertex;
  
  // Get the the projected length                                                                 
  double projLen = pos.Dot(direction);

  return projLen;
}

double shower::SBNShowerAlg::SpacePointPerpendiular(const art::Ptr<recob::SpacePoint>&sp,
						    TVector3& vertex, TVector3& direction, 
						    double proj){

  // Get the position of the spacepoint                                                            
  TVector3 pos = shower::SBNShowerAlg::SpacePointPosition(sp) - vertex;

  // Take away the projection * distance to find the perpendicular vector

  pos = pos - proj * direction;
  
  // Get the the projected length                                                
  double perpLen = pos.Mag();
  
  return perpLen;
}


void shower::SBNShowerAlg::TrackValidationPlotter(const art::Ptr<recob::PFParticle>& pfparticle,
						  art::Event& Event,
						  reco::shower::ShowerPropertyHolder& ShowerPropHolder){

  //Function for drawing reco showers to check direction and initial track selection

  // Get run info to make unique canvas names
  int run    = Event.run();
  int subRun = Event.subRun();
  int event  = Event.event();
  int PFPID  = pfparticle->Self();

  // Create the canvas
  TString canvasName = Form("canvas_%i_%i_%i_%i",run,subRun,event,PFPID);
  TCanvas* canvas = tfs->make<TCanvas>(canvasName, canvasName);
    
  // Initialise variables
  float x;
  float y;
  float z;

  // Get a bunch of associations (again)
  // N.B. this is a horribly inefficient way of doing things but as this is only
  // going to be used to debug I don't care, I would rather have generality in this case

  art::Handle<std::vector<recob::PFParticle> > pfpHandle;
  if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
    throw cet::exception("Shower3DTrackFinderEMShower") << "Could not get the pandora pf particles\
. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
    return;
  }

  // Get the spacepoint - PFParticle assn                                                          
  art::FindManyP<recob::SpacePoint> fmspp(pfpHandle, Event, fPFParticleModuleLabel);
  if (!fmspp.isValid()){
    throw cet::exception("Shower3DTrackFinder") << "Trying to get the spacepoint and failed. Somet\
hing is not configured correctly. Stopping ";
    return;
  }

  // Get the SpacePoints                                                                           
  std::vector<art::Ptr<recob::SpacePoint> > spacePoints = fmspp.at(pfparticle.key());

  // Get the hit handle                                                                           
  art::Handle<std::vector<recob::Hit> > hitHandle;
  if (!Event.getByLabel(fHitModuleLabel, hitHandle)){
    throw cet::exception("Shower3DTrackFinder") << "Could not configure the hit handle. Something is configured incorrectly. Stopping";
    return;
  }

  // Get the hits associated with the space points                                                 
  art::FindManyP<recob::SpacePoint> fmsph(hitHandle, Event, fPFParticleModuleLabel);
  if(!fmsph.isValid()){
    throw cet::exception("Shower3DTrackFinderPosition") << "Spacepoint and hit association not val\
id. Stopping.";
    return;
  }

  //We cannot progress with no spacepoints.                                                        
  if(spacePoints.size() == 0){
    //throw cet::exception("Shower3DTrackFinder") << "No Space Points. Stopping.";
    return;
  }

  if(!ShowerPropHolder.CheckShowerStartPosition()){
    mf::LogError("Shower3DTrackFinder") << "Start position not set, returning "<< std::endl;
    return;
  }
  if(!ShowerPropHolder.CheckShowerDirection()){
    mf::LogError("Shower3DTrackFinder") << "Direction not set, returning "<< std::endl;
    return;
  }
  if(!ShowerPropHolder.CheckInitialTrackHits()){
    mf::LogError("Shower3DTrackFinder") << "TrackHits not set, returning "<< std::endl;
    return;
  }

  // Get info from shower property holder
  TVector3 showerStartPosition = ShowerPropHolder.GetShowerStartPosition();
  TVector3 showerDirection     = ShowerPropHolder.GetShowerDirection();
  std::vector<art::Ptr<recob::Hit> > trackHits=ShowerPropHolder.GetInitialTrackHits();
  
  // Create 3D point at vertex, chosed to be origin for ease of use of display
  double startXYZ[3] = {0,0,0};
  TPolyMarker3D* startPoly = new TPolyMarker3D(1,startXYZ);
  
  // get the space points associated to the initial track hits
  std::vector<art::Ptr<recob::SpacePoint> > trackSpacePoints;
  for (auto hit : trackHits){
    const std::vector<art::Ptr<recob::SpacePoint> > sps = fmsph.at(hit.key());
    const art::Ptr<recob::SpacePoint> sp = sps.front();
    trackSpacePoints.push_back(sp);
  }

  // Get the min and max projections along the direction to know how long to draw 
  // the direction line
  double minProj=0;
  double maxProj=0;

  //initialise counter point
  int point = 0;
  // Make 3D points for each spacepoint in the shower
  TPolyMarker3D* allPoly = new TPolyMarker3D(spacePoints.size());
  for (auto spacePoint : spacePoints){
    TVector3 pos = shower::SBNShowerAlg::SpacePointPosition(spacePoint) - showerStartPosition;

    x = pos.X();
    y = pos.Y();
    z = pos.Z();
    allPoly->SetPoint(point,x,y,z);
    ++point;
    
    // Calculate the projection of (point-startpoint) along the direction
    double proj = shower::SBNShowerAlg::SpacePointProjection(spacePoint, showerStartPosition, 
							     showerDirection);
    if (proj>maxProj) {
      maxProj = proj;
    } else if (proj<minProj) {
      minProj = proj ;
    }

  } // loop over spacepoints
  
  // Create TPolyLine3D arrays
  double xDirPoints[3] = {minProj*showerDirection.X(), 0, maxProj*showerDirection.X()};
  double yDirPoints[3] = {minProj*showerDirection.Y(), 0, maxProj*showerDirection.Y()};
  double zDirPoints[3] = {minProj*showerDirection.Z(), 0, maxProj*showerDirection.Z()};

  TPolyLine3D* dirPoly = new TPolyLine3D(3,xDirPoints,yDirPoints,zDirPoints);			 

  point = 0; // re-initialise counter
  TPolyMarker3D* trackPoly = new TPolyMarker3D(trackSpacePoints.size());
  for (auto spacePoint : trackSpacePoints){
    TVector3 pos = shower::SBNShowerAlg::SpacePointPosition(spacePoint) - showerStartPosition;
    x = pos.X();
    y = pos.Y();
    z = pos.Z();
    trackPoly->SetPoint(point,x,y,z);    
    ++point;
  } // loop over track spacepoints

  // TODO: make this a fcl parameter
  bool fDrawAllPFPs = true;

  // If we want to draw all of the PFParticles in the event
  if (fDrawAllPFPs){
    //Get the PFParticles                                                                           
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    std::vector<art::Ptr<recob::PFParticle> > pfps;
    if (Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      art::fill_ptr_vector(pfps, pfpHandle);
    }
    else {
      throw cet::exception("SBNShower") << "pfps not loaded." << std::endl;
    }
    // initialse counters
    // Split into tracks and showers to make it clearer what pandora is doing
    int pfpTrackCounter = 0;
    int pfpShowerCounter = 0;
    
    // initial loop over pfps to find nuber of spacepoints for tracks and showers
    for(auto const& pfp: pfps){
      std::vector<art::Ptr<recob::SpacePoint> > sps = fmspp.at(pfp.key());
      if (pfp->PdgCode()==11) { pfpShowerCounter += sps.size();
      } else if (pfp->PdgCode()==13) { pfpTrackCounter += sps.size(); }
    }

    TPolyMarker3D* pfpPolyTrack = new TPolyMarker3D(pfpTrackCounter);
    TPolyMarker3D* pfpPolyShower = new TPolyMarker3D(pfpShowerCounter);

    // initialise counters
    int trackPoints  = 0;
    int showerPoints = 0;

    for(auto const& pfp: pfps){
      std::vector<art::Ptr<recob::SpacePoint> > sps = fmspp.at(pfp.key());
      int pdg = pfp->PdgCode(); // Track or shower
      for (auto sp : sps){
	TVector3 pos = shower::SBNShowerAlg::SpacePointPosition(sp) - showerStartPosition;
	x = pos.X();
	y = pos.Y();
	z = pos.Z();
	if (pdg==11){
	  pfpPolyShower->SetPoint(showerPoints,x,y,z);
	  ++showerPoints;
	} else if (pdg==13){
	  pfpPolyTrack->SetPoint(trackPoints,x,y,z);
	  ++trackPoints;
	}
      } // loop over sps
    } // loop over pfps
    pfpPolyShower->SetMarkerStyle(20);
    pfpPolyShower->SetMarkerColor(4);
    pfpPolyShower->Draw();
    pfpPolyTrack->SetMarkerStyle(20);
    pfpPolyTrack->SetMarkerColor(6);
    pfpPolyTrack->Draw();
    
  } // if (fDrawAllPFPs)

  // Draw all of the things
  allPoly->SetMarkerStyle(20);
  allPoly->Draw();
  trackPoly->SetMarkerStyle(20);
  trackPoly->SetMarkerColor(2);
  trackPoly->Draw();
  startPoly->SetMarkerStyle(21);
  startPoly->SetMarkerSize(2);
  startPoly->SetMarkerColor(3);
  startPoly->Draw();
  dirPoly->SetLineWidth(1);
  dirPoly->SetLineColor(6);
  dirPoly->Draw();

  // Save the canvas. Don't usually need this when using TFileService but this in the alg
  // not a module and didn't work without this so im going with it. 
  canvas->Write();
}


