//############################################################################
//### Name:        Shower3DTrackFinder                                     ###
//### Author:      Ed Tyley                                                ###
//### Date:        14.06.19                                                ###
//### Description: Tool for finding the initial shower track using 3D      ###
//###              spacepoints within a cylinder along the shower          ###
//###              direction. fcl parameters define cylinder dimensions    ###
//############################################################################
#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Core/EDProducer.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft Includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/SBNShowerAlg.h"
#include "larreco/RecoAlg/ProjectionMatchingAlg.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

//C++ Includes
#include <iostream>
#include <math.h>

//Root Includes
#include "TVector3.h"

namespace ShowerRecoTools{

  class Shower3DTrackFinder:IShowerTool {
  public:

    Shower3DTrackFinder(const fhicl::ParameterSet& pset);

    ~Shower3DTrackFinder();

    //Generic Track Finder
    int CalculateProperty(const art::Ptr<recob::PFParticle>& pfparticle,
			  art::Event& Event,
			  reco::shower::ShowerPropertyHolder& ShowerPropHolder
			  ) override;

  private:

    void configure(const fhicl::ParameterSet& pset) override;

    std::vector<art::Ptr<recob::Hit> > FindInitialTrackHits(std::vector<art::Ptr<recob::Hit> >& hits, TVector3& ShowerStartPosition, TVector3& ShowerDirection);

    std::vector<art::Ptr<recob::SpacePoint> > FindTrackSpacePoints(std::vector<art::Ptr<recob::SpacePoint> >& spacePoints, TVector3& showerStartPosition,TVector3& showerDirection);

    Int_t WeightedFit(const Int_t n, const Double_t *x, const Double_t *y, const Double_t *w,  Double_t *parm);

    // Define standard art tool interface

    shower::SBNShowerAlg       fSBNShowerAlg;
    pma::ProjectionMatchingAlg fProjectionMatchingAlg;

    float fMaxProjectionDist;    // Maximum projection along shower direction, length of cylinder
    float fMaxPerpendicularDist; // Maximum perpendicular distance, radius of cylinder
    bool fForwardHitsOnly;       // Only take hits downstream of shower vertex (projection>0)
    bool fThreeDTrackFinding;    // Fcl to allow comparisons between 2D and 3D

    // For 2D method for comparison
    unsigned int               fNfitpass;
    std::vector<unsigned int>  fNfithits;
    std::vector<double>        fToler;
    bool                       fApplyChargeWeight;
    bool                       fDebugEVD;
    art::InputTag              fPFParticleModuleLabel;
    detinfo::DetectorProperties const* fDetProp;
    art::ServiceHandle<geo::Geometry> fGeom;
  };


  Shower3DTrackFinder::Shower3DTrackFinder(const fhicl::ParameterSet& pset)
    : fSBNShowerAlg(pset.get<fhicl::ParameterSet>("SBNShowerAlg")),
      fProjectionMatchingAlg(pset.get<fhicl::ParameterSet>("ProjectionMatchingAlg")),
      fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())
  {
    configure(pset);
  }

  Shower3DTrackFinder::~Shower3DTrackFinder()
  {
  }

  void Shower3DTrackFinder::configure(const fhicl::ParameterSet& pset)
  {
    fApplyChargeWeight     = pset.get<bool>                      ("ApplyChargeWeight");
    fDebugEVD              = pset.get<bool>                      ("DebugEVD");  
    fNfitpass              = pset.get<unsigned int>              ("Nfitpass");
    fNfithits              = pset.get<std::vector<unsigned int> >("Nfithits");
    fToler                 = pset.get<std::vector<double> >      ("Toler");
    fPFParticleModuleLabel = pset.get<art::InputTag>             ("PFParticleModuleLabel");

    fMaxProjectionDist     = pset.get<float> ("MaxProjectionDist");
    fMaxPerpendicularDist  = pset.get<float> ("MaxPerpendicularDist");
    fForwardHitsOnly       = pset.get<bool>  ("ForwardHitsOnly");
    fThreeDTrackFinding    = pset.get<bool>  ("ThreeDTrackFinding");


    if (fNfitpass!=fNfithits.size() ||
	fNfitpass!=fToler.size()) {
      throw art::Exception(art::errors::Configuration)
	<< "Shower3DTrackFinderEMShower: fNfithits and fToler need to have size fNfitpass";
    }
  }


  int Shower3DTrackFinder::CalculateProperty(const art::Ptr<recob::PFParticle>& pfparticle,
					   art::Event& Event,
					   reco::shower::ShowerPropertyHolder& ShowerPropHolder
					   ){
    std::cout <<"#########################################\n"<<
      "hello world track finder 3d\n" <<"#########################################\n"<< std::endl;

    //This is all based on the shower vertex being known. If it is not lets not do the track
    if(!ShowerPropHolder.CheckShowerStartPosition()){
      mf::LogError("Shower3DTrackFinder") << "Start position not set, returning "<< std::endl;
      return 1;
    }
    if(!ShowerPropHolder.CheckShowerDirection()){
      mf::LogError("Shower3DTrackFinder") << "Direction not set, returning "<< std::endl;
      return 1;
    }

    TVector3 ShowerStartPosition = ShowerPropHolder.GetShowerStartPosition();
    TVector3 ShowerDirection     = ShowerPropHolder.GetShowerDirection();

    // Get the assocated pfParicle Handle
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("Shower3DTrackFinderEMShower") << "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    // Get the spacepoint - PFParticle assn
    art::FindManyP<recob::SpacePoint> fmspp(pfpHandle, Event, fPFParticleModuleLabel);
    if (!fmspp.isValid()){
      throw cet::exception("Shower3DTrackFinder") << "Trying to get the spacepoint and failed. Something is not configured correctly. Stopping ";
      return 1;
    }

    // Get the spacepoints
    art::Handle<std::vector<recob::SpacePoint> > spHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, spHandle)){
      throw cet::exception("Shower3DTrackFinder") << "Could not configure the spacepoint handle. Something is configured incorrectly. Stopping";
      return 1;
    }

    // Get the hits associated with the space points
    art::FindManyP<recob::Hit> fmhsp(spHandle, Event, fPFParticleModuleLabel);
    if(!fmhsp.isValid()){
      throw cet::exception("Shower3DTrackFinderPosition") << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }



    // Get the SpacePoints
    std::vector<art::Ptr<recob::SpacePoint> > spacePoints = fmspp.at(pfparticle.key());

    //We cannot progress with no spacepoints.
    if(spacePoints.size() == 0){
       throw cet::exception("Shower3DTrackFinder") << "No Space Points. Stopping.";
      return 1;
    }
    
    std::cout<<"Spacepoints: "<<spacePoints.size()<<std::endl;

    // Order the spacepoints
    fSBNShowerAlg.OrderShowerSpacePoints(spacePoints,ShowerStartPosition,ShowerDirection);

    // If we are selecting the track hits in 3D do it before spliting into planes
    std::vector<art::Ptr<recob::SpacePoint> > trackSpacePoints;
    if (fThreeDTrackFinding) {
      trackSpacePoints = FindTrackSpacePoints(spacePoints,ShowerStartPosition,ShowerDirection);
    };

    spacePoints = trackSpacePoints; //TODO: think of a better way (maybe remove with 2D)

    // Get the hits associated to the space points and seperate them by planes
    std::map<geo::PlaneID, std::vector<art::Ptr<recob::Hit> > > plane_hits;
    for(auto const& spacePoint: spacePoints){

      //Get the hit
      const std::vector<art::Ptr<recob::Hit> > hits = fmhsp.at(spacePoint.key());
      // Should I be using FindOneP instead of FindManyP here?
      if( hits.size() !=1 ){
	mf::LogError("Shower3DTrackFinder") << "Wrong number of hits associated to spacepoint: "
					    << hits.size() << std::endl;
	return 1;
      }

      art::Ptr<recob::Hit> hit = hits.front();

      //Get the view.
      //geo::PlaneID plane = cluster->Plane();

      geo::WireID wire = hit->WireID();
      geo::PlaneID plane = wire.asPlaneID();
      plane_hits[plane].push_back(hit);

    }

    


    std::map<geo::PlaneID, std::vector<art::Ptr<recob::Hit> > > plane_trackhits;
    if (fThreeDTrackFinding) {
      plane_trackhits = plane_hits;
    } else {
      //Loop over the clusters and order the hits and get the initial track hits in that plane
      for(auto const& plane_hits_iter: plane_hits){
	//Get the Plane
	geo::PlaneID plane = plane_hits_iter.first;
	//Get the hits
	std::vector<art::Ptr<recob::Hit> > hits = plane_hits_iter.second;

      //Find the initial track hits
	std::vector<art::Ptr<recob::Hit> > trackhits = FindInitialTrackHits(hits,ShowerStartPosition,ShowerDirection);

	plane_trackhits[plane].insert(plane_trackhits[plane].end(),trackhits.begin(),trackhits.end());
      }
    }

    //Decide which plane to use
    int maxhits       = 1; // set to 1 as we require at least 2 track hits per plane
    geo::TPCID vtxTPC = fGeom->FindTPCAtPosition(ShowerStartPosition);
    geo::PlaneID maxplane; // Note this contains both tpc and cryostat information

    for(auto const& plane : plane_trackhits){
      std::vector<art::Ptr<recob::Hit> >  trackhits = plane.second;
      geo::TPCID maxTPC = (plane.first).asTPCID();
      if( maxTPC == vtxTPC){
	if((int) trackhits.size() > maxhits ){
	  maxplane = plane.first;
	  maxhits  = trackhits.size();
	}
      }
    }

    if( maxhits == 1 || !maxplane){
      mf::LogError("Shower3DTrackFinder") << "Max Plane not set " << std::endl;
      return 1;
    }

    int nextmaxhits  = 1;
    geo::PlaneID nextmaxplane;

    for(auto const& plane : plane_trackhits){
      //Check clusters are not in same plane
      if( (plane.first) == maxplane){continue;}
      //Need to make sure clusters are in same tpc
      geo::TPCID nextmaxTPC = (plane.first).asTPCID();
      if( nextmaxTPC == vtxTPC){
	std::vector<art::Ptr<recob::Hit> > trackhits = plane.second;
	if((int) trackhits.size() > nextmaxhits){
	  nextmaxplane = plane.first;
	  nextmaxhits  = trackhits.size();
	}
      }
    }


    if( nextmaxhits == 1 || !nextmaxplane){
      mf::LogError("Shower3DTrackFinder") << "Next Max Plane not set " << std::endl;
      return 1;
    }

    std::vector<art::Ptr<recob::Hit> > maxPlaneHits = (plane_trackhits.at(maxplane));
    std::vector<art::Ptr<recob::Hit> > nextmaxPlaneHits = (plane_trackhits.at(nextmaxplane));

    if(maxPlaneHits.size() < 2 && nextmaxPlaneHits.size() < 2){
      mf::LogWarning("Shower3DTrackFinder") << "Not Enough Hits" << std::endl;
      return 0;
    }

    //Build the 3D track
    pma::Track3D* pmatrack = fProjectionMatchingAlg.buildSegment(maxPlaneHits, nextmaxPlaneHits, ShowerStartPosition);

    if(!pmatrack){
      mf::LogWarning("Shower3DTrackFinder") << "No PMA track made " << std::endl;
      return 1;
    }

    //Get the spacepoints
    std::vector<TVector3> spts;
    for (size_t i = 0; i<pmatrack->size(); ++i){
      if ((*pmatrack)[i]->IsEnabled()){
	TVector3 p3d = (*pmatrack)[i]->Point3D();
	spts.push_back(p3d);
      }
    }

    if(spts.size() < 2){
      mf::LogWarning("Shower3DTrackFinder") << "Not Enough Spacepoints" << std::endl;
      return 1;
    }

    //Make the track. Adapted from PandoraShowerCreation
    recob::tracking::Positions_t xyz;
    recob::tracking::Momenta_t pxpypz;
    recob::TrackTrajectory::Flags_t flags;

    TVector3 spt1 = spts[0];

    for(unsigned int sp=0; sp<spts.size(); ++sp){

      TVector3 spt = spts[sp];

      if(sp < spts.size() -1){
	spt1 = spts[sp+1];
      }
      else{
	spt1 = spts[sp];
      }

      xyz.emplace_back(recob::tracking::Point_t(spt.X(), spt.Y(), spt.Z()));

      TVector3 dir = -(spt - spt1);

      pxpypz.emplace_back(recob::tracking::Vector_t(dir.X(), dir.Y(), dir.Z()));

      if (std::fabs(spt.X()-util::kBogusF)<std::numeric_limits<float>::epsilon() &&
	    std::fabs(spt.X()-util::kBogusF)<std::numeric_limits<float>::epsilon() &&
	  std::fabs(spt.X()-util::kBogusF)<std::numeric_limits<float>::epsilon())
	{
	  flags.emplace_back(recob::TrajectoryPointFlags(recob::TrajectoryPointFlags::InvalidHitIndex, recob::TrajectoryPointFlagTraits::NoPoint));
	}
      else {
	flags.emplace_back(recob::TrajectoryPointFlags());
      }
      spt1 = spt;
    }

    std::vector<art::Ptr<recob::Hit> > TrackHits;
    for(auto const& trackhits_p: plane_trackhits){
      TrackHits.insert(TrackHits.end(),trackhits_p.second.begin(),trackhits_p.second.end());
    }

    //Actually make the thing.
    recob::Track track = recob::Track(recob::TrackTrajectory(std::move(xyz), std::move(pxpypz), std::move(flags), false),
				      util::kBogusI, util::kBogusF, util::kBogusI, recob::tracking::SMatrixSym55(),
				      recob::tracking::SMatrixSym55(), pfparticle.key());

    ShowerPropHolder.SetInitialTrack(track);
    ShowerPropHolder.SetInitialTrackHits(TrackHits);

    if (fDebugEVD){
	std::cout<<"Do DebugEVD"<<std::endl;
	fSBNShowerAlg.DebugEVD(pfparticle,Event,ShowerPropHolder);
    }

    std::cout <<"#########################################\n"<<
      "track finder 3d done\n" <<"#########################################\n"<< std::endl;

    return 0;
  }

  std::vector<art::Ptr<recob::SpacePoint> > Shower3DTrackFinder::FindTrackSpacePoints(std::vector<art::Ptr<recob::SpacePoint> >& spacePoints, TVector3& showerStartPosition,TVector3& showerDirection){

    std::vector<art::Ptr<recob::SpacePoint> > trackSpacePoints;

    for (auto spacePoint : spacePoints){
      double proj = fSBNShowerAlg.SpacePointProjection(spacePoint, showerStartPosition, showerDirection);
      double perp = fSBNShowerAlg.SpacePointPerpendiular(spacePoint, showerStartPosition, showerDirection, proj);
      
      //std::cout<<"Proj: "<<proj<<", Perp: "<<perp<<std::endl;
      if (fForwardHitsOnly){
    	if (proj>0 && proj<fMaxProjectionDist && TMath::Abs(perp)<fMaxPerpendicularDist){
          trackSpacePoints.push_back(spacePoint);
	    }
      } else {
        if (TMath::Abs(proj)<fMaxProjectionDist && TMath::Abs(perp)<fMaxPerpendicularDist){
	    trackSpacePoints.push_back(spacePoint);
	    }
      }
    }

    return trackSpacePoints;
  }


  //Function to calculate the what are the initial tracks hits. Adapted from EMShower FindInitialTrackHits
  std::vector<art::Ptr<recob::Hit> > Shower3DTrackFinder::FindInitialTrackHits(std::vector<art::Ptr<recob::Hit> >& hits, TVector3& ShowerStartPosition, TVector3& ShowerDirection){

    std::vector<art::Ptr<recob::Hit> > trackHits;

    double parm[2];
    int fitok = 0;
    std::vector<double> wfit;
    std::vector<double> tfit;
    std::vector<double> cfit;

    //std::cout<<"fNfitpass "<<fNfitpass<<std::endl;

    for (size_t i = 0; i<fNfitpass; ++i){


      // Fit a straight line through hits
      unsigned int nhits = 0;
      for (auto &hit: hits){

	//Not sure I am a fan of doing things in wire tick space. What if id doesn't not iterate properly or the
	//two planes in each TPC are not symmetric.
	TVector2 coord = fSBNShowerAlg.HitCoordinates(hit);

	//#####################################################
	//check the magnitude
	// TVector2 coord_hit = coord - Shower2DStartPosition;
	// 	double proj = TMath::Abs(coord_hit.X()*Shower2DDirection.X() +
	// 			 coord_hit.Y()*Shower2DDirection.Y());
	// TVector2 perp = coord_hit - proj*Shower2DDirection;
	// double  len_perp = perp.Mod();

	//#####################################################

	//std::cout<<i<<" "<<hit->WireID()<<" "<<hit->PeakTime()<<std::endl;

	if (i==0||(std::abs((coord.Y()-(parm[0]+coord.X()*parm[1]))*cos(atan(parm[1])))<fToler[i-1])||fitok==1){
	  //	  std::cout << "SBN Shower passed" << std::endl;
	  ++nhits;
	  if (nhits==fNfithits[i]+1) break;
	  wfit.push_back(coord.X());
	  tfit.push_back(coord.Y());

	  if(fApplyChargeWeight) { cfit.push_back(hit->Integral());
	  } else { cfit.push_back(1.); };
	  if (i==fNfitpass-1) {
	    trackHits.push_back(hit);
	  }
	}
      }

      if (i<fNfitpass-1&&wfit.size()){
	fitok = WeightedFit(wfit.size(), &wfit[0], &tfit[0], &cfit[0], &parm[0]);
      }

      wfit.clear();
      tfit.clear();
      cfit.clear();
    }
    return trackHits;
  }

  //Stolen from EMShowerAlg, a regression fitting function
  Int_t Shower3DTrackFinder::WeightedFit(const Int_t n, const Double_t *x, const Double_t *y, const Double_t *w,  Double_t *parm){

  Double_t sumx=0.;
  Double_t sumx2=0.;
  Double_t sumy=0.;
  Double_t sumy2=0.;
  Double_t sumxy=0.;
  Double_t sumw=0.;
  Double_t eparm[2];

  parm[0]  = 0.;
  parm[1]  = 0.;
  eparm[0] = 0.;
  eparm[1] = 0.;

  for (Int_t i=0; i<n; i++) {
    sumx += x[i]*w[i];
    sumx2 += x[i]*x[i]*w[i];
    sumy += y[i]*w[i];
    sumy2 += y[i]*y[i]*w[i];
    sumxy += x[i]*y[i]*w[i];
    sumw += w[i];
  }

  if (sumx2*sumw-sumx*sumx==0.) return 1;
  if (sumx2-sumx*sumx/sumw==0.) return 1;

  parm[0] = (sumy*sumx2-sumx*sumxy)/(sumx2*sumw-sumx*sumx);
  parm[1] = (sumxy-sumx*sumy/sumw)/(sumx2-sumx*sumx/sumw);

  eparm[0] = sumx2*(sumx2*sumw-sumx*sumx);
  eparm[1] = (sumx2-sumx*sumx/sumw);

  if (eparm[0]<0. || eparm[1]<0.) return 1;

  eparm[0] = sqrt(eparm[0])/(sumx2*sumw-sumx*sumx);
  eparm[1] = sqrt(eparm[1])/(sumx2-sumx*sumx/sumw);

  return 0;

  }





//#####################################################
//check the magnitude
// TVector2 coord_hit = coord - Shower2DStartPosition;
// double proj = TMath::Abs(coord_hit.X()*Shower2DDirection.X() +
//			  coord_hit.Y()*Shower2DDirection.Y());
// TVector2 perp = coord_hit - proj*Shower2DDirection;
// double  len_perp = perp.Mod();
//#####################################################

}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::Shower3DTrackFinder)




