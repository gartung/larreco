//############################################################################
//### Name:        ShowerTrackFinder                                       ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the initial shower track using a rms   ###
//###              based method to define when the shower starts to        ###
//###              shower. This methd is derived from the EMShower_module  ###
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
#include "art/Persistency/Common/PtrMaker.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft Includes 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "larreco/RecoAlg/SBNShowerAlg.h"
#include "larreco/RecoAlg/ProjectionMatchingAlg.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "lardata/Utilities/AssociationUtil.h"

//C++ Includes 
#include <iostream>
#include <math.h>

namespace ShowerRecoTools{

  class ShowerTrackFinder:IShowerTool {
  public:

    ShowerTrackFinder(const fhicl::ParameterSet& pset);
    
    ~ShowerTrackFinder(); 
    
    //Generic Track Finder
    int findMetric(const art::Ptr<recob::PFParticle>& pfparticle,
		   art::Event& Event,
		   reco::shower::ShowerPropertyHolder& ShowerPropHolder
		   ) override;

  private:  
    
    void configure(const fhicl::ParameterSet& pset) override;

    std::vector<art::Ptr<recob::Hit> > FindInitialTrackHits(std::vector<art::Ptr<recob::Hit> >& hits, TVector3& ShowerStartPosition, TVector3& ShowerDirection);

    Int_t WeightedFit(const Int_t n, const Double_t *x, const Double_t *y, const Double_t *w,  Double_t *parm);

    // Define standard art tool interface

    shower::SBNShowerAlg fSBNShowerAlg;
    pma::ProjectionMatchingAlg fProjectionMatchingAlg;
    unsigned int         fNfitpass;
    std::vector<unsigned int>     fNfithits;
    std::vector<double>  fToler;
    bool fApplyChargeWeight;
    art::InputTag fPFParticleModuleLabel;

    // dEdx params
    art::ServiceHandle<geo::Geometry> fGeom;
    calo::CalorimetryAlg fCalorimetryAlg;
    double fdEdxTrackLength;
    bool fMaxHitPlane;
    bool fMissFirstPoint;
  };
  
  
  ShowerTrackFinder::ShowerTrackFinder(const fhicl::ParameterSet& pset)
    : fSBNShowerAlg(pset.get<fhicl::ParameterSet>("SBNShowerAlg")),
      fProjectionMatchingAlg(pset.get<fhicl::ParameterSet>("ProjectionMatchingAlg")),
      fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
  {

    //    produces<art::Assns<recob::Track, recob::Hit> >();
    configure(pset);
  }
  
  ShowerTrackFinder::~ShowerTrackFinder()
  {
  }
  
  void ShowerTrackFinder::configure(const fhicl::ParameterSet& pset)
  {
    fApplyChargeWeight     = pset.get<bool>                      ("ApplyChargeWeight");
    fNfitpass              = pset.get<unsigned int>              ("Nfitpass");
    fNfithits              = pset.get<std::vector<unsigned int> >("Nfithits");
    fToler                 = pset.get<std::vector<double> >      ("Toler");
    fPFParticleModuleLabel = pset.get<art::InputTag>             ("PFParticleModuleLabel");
    fdEdxTrackLength       = pset.get<double>                    ("dEdxTrackLength");
    fMaxHitPlane           = pset.get<bool>                      ("MaxHitPlane");
    fMissFirstPoint        = pset.get<bool>                      ("MissFirstPoint");


    if (fNfitpass!=fNfithits.size() ||
	fNfitpass!=fToler.size()) {
      throw art::Exception(art::errors::Configuration)
	<< "ShowerTrackFinderEMShower: fNfithits and fToler need to have size fNfitpass";
    }
  }
  
  
  int ShowerTrackFinder::findMetric(const art::Ptr<recob::PFParticle>& pfparticle,
					    art::Event& Event,
					    reco::shower::ShowerPropertyHolder& ShowerPropHolder
					    ){
    std::cout <<"#########################################\n"<<
      "hello world track finder\n" <<"#########################################\n"<< std::endl;

    //This is all based on the shower vertex being known. If it is not lets not do the track 
    if(!ShowerPropHolder.CheckShowerStartPosition()){
      mf::LogError("ShowerTrackFinder") << "Start position not set, returning "<< std::endl;
      return 0;
    }
    if(!ShowerPropHolder.CheckShowerDirection()){
      mf::LogError("ShowerTrackFinder") << "Direction not set, returning "<< std::endl;
      return 0;
    }
    
    //    auto assnhit = std::make_unique<art::Assns<recob::Track, recob::Hit>>();

    TVector3 ShowerStartPosition = ShowerPropHolder.GetShowerStartPosition();
    TVector3 ShowerDirection     = ShowerPropHolder.GetShowerDirection();

    // Get the assocated pfParicle vertex PFParticles
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("ShowerTrackFinderEMShower") << "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }
    
    //Get the clusters
    art::Handle<std::vector<recob::Cluster> > clusHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, clusHandle)){
      throw cet::exception("ShowerTrackFinderEMShower") << "Could not get the pandora clusters. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }
    art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleModuleLabel);
    std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(pfparticle.key());
    

    if(clusters.size()<2){
      mf::LogError("ShowerTrackFinder") << "Not enough clusters: "<<clusters.size() << std::endl;
      // throw cet::exception("ShowerTrackFinderEMShower") << "Not enough clusters: "
      // 						<<clusters.size(<<" ");
      return 1;
    }

    //Get the hit association 
    art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleModuleLabel);

    std::map<geo::PlaneID, std::vector<art::Ptr<recob::Hit> > > plane_clusters;
    std::cout<<"Num Clusters "<<clusters.size()<<std::endl;;    
    //Loop over the clusters in the plane and get the hits 
    for(auto const& cluster: clusters){
      
      //Get the hits 
      std::vector<art::Ptr<recob::Hit> > hits = fmhc.at(cluster.key());
      
      //Get the view.
      //geo::PlaneID plane = cluster->Plane();
      
      for (auto hit : hits) {
	geo::WireID wire = hit->WireID();
	geo::PlaneID plane = wire.asPlaneID();
	plane_clusters[plane].push_back(hit);
      }
      
      // Was having issues with clusters having hits in multiple planes breaking PMA
      // So switched to the method above. May want to switch back when using PandoraTrack
      //plane_clusters[plane].insert(plane_clusters[plane].end(),hits.begin(),hits.end());
    }

    std::map<geo::PlaneID, std::vector<art::Ptr<recob::Hit> > > plane_trackhits;
    //std::vector<std::vector<art::Ptr<recob::Hit> > > trackhits;
    //Loop over the clusters and order the hits and get the initial track hits in that plane
    for(auto const& cluster: plane_clusters){
      //Get the Plane
      geo::PlaneID plane = cluster.first;
      //Get the hits 
      std::vector<art::Ptr<recob::Hit> > hits = cluster.second;    

      //Order the hits 
      fSBNShowerAlg.OrderShowerHits(hits,ShowerStartPosition,ShowerDirection);


      //Find the initial track hits 
      std::vector<art::Ptr<recob::Hit> > trackhits = FindInitialTrackHits(hits,ShowerStartPosition,ShowerDirection);
      //view_trackhits.push_back(trackhits_plane);

      std::cout<<"Plane: "<<plane.toString()<<" with hits: "<<hits.size()
	       <<"  and track hits:"<<trackhits.size()<<std::endl;

      /*
      for (auto hit : trackhits){
	geo::View_t view = hit->View();
	geo::WireID wire = hit->WireID();
	double angle = fGeom->WireAngleToVertical(view,wire.asPlaneID());
	
	std::cout<<wire.toString()<<" with angle: "<<angle<<" View: "<<view;
       
        if (view==geo::kU) { std::cout<<" kU"<<std::endl; };
        if (view==geo::kV) { std::cout<<" kV"<<std::endl; };
        if (view==geo::kW) { std::cout<<" kW"<<std::endl; };
      }
      */
      plane_trackhits[plane].insert(plane_trackhits[plane].end(),trackhits.begin(),trackhits.end());
    }

    //Decide which plane to use 
    //int maxplane = -1;
    int maxhits       = 1; // set to 1 as we require at least 2 track hits per plane
    geo::TPCID vtxTPC = fGeom->FindTPCAtPosition(ShowerStartPosition); 
    geo::PlaneID maxplane; // Note this contains both tpc and cryostat information

    std::cout<<"Shower Vertex: X:"<<ShowerStartPosition.X()<<" Y: "<<ShowerStartPosition.Y()
	     <<" Z: "<<ShowerStartPosition.Z()<<" in TPC: "<<vtxTPC.toString()<<std::endl;


    for(auto const& plane : plane_trackhits){
      std::vector<art::Ptr<recob::Hit> >  trackhits = plane.second;    
      geo::TPCID maxTPC = (plane.first).asTPCID();
      if( maxTPC == vtxTPC){
	if((int) trackhits.size() > maxhits ){
	  maxplane = plane.first;
	  maxhits  = trackhits.size();
	  //maxTPC   = (plane.first).asTPCID(); // convert geo::PlaneID to geo::TPCID
	}
      }
    } 
    
    if( maxhits == 1 || !maxplane){
      std::cout<<"Next Max Plane not set "<<std::endl;
      mf::LogError("ShowerTrackFinder") << "Max Plane not set " << std::endl;
      //throw cet::exception("ShowerTrackFinderEMShower") << "Max Plane not set ";
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
      mf::LogError("ShowerTrackFinder") << "Next Max Plane not set " << std::endl;
      //throw cet::exception("ShowerTrackFinderEMShower") << "Next Max Plane not set ";
      return 1;
    }
    
    // Trying to debug PMA error where it claims there are not hits in 2 planes
    std::cout<<"Max Plane: "<<maxplane.toString()<<" with "<<maxhits<<std::endl;
    std::cout<<"Next MaxPlane: "<<nextmaxplane.toString()<<" with "<<nextmaxhits<<std::endl;

    std::vector<art::Ptr<recob::Hit> > maxPlaneHits = (plane_trackhits.at(maxplane));
    std::vector<art::Ptr<recob::Hit> > nextmaxPlaneHits = (plane_trackhits.at(nextmaxplane));

    // art::Ptr<recob::Hit> maxPlaneHit = maxPlaneHits.at(0);
    // art::Ptr<recob::Hit> nextmaxPlaneHit = nextmaxPlaneHits.at(0);
    
    // std::cout<<"planes: "<<maxPlaneHit->View()<<" and "<<nextmaxPlaneHit->View() <<std::endl;
    // std::cout<<"Max Plane"<<std::endl;
    // for (auto hit : maxPlaneHits){
    //   geo::View_t view = hit->View(); 
    //   std::cout<<(hit->WireID()).toString()<<" View: "<<view;
    //   if (view==geo::kU) { std::cout<<" kU"<<std::endl; };
    //   if (view==geo::kV) { std::cout<<" kV"<<std::endl; };
    //   if (view==geo::kW) { std::cout<<" kW"<<std::endl; };
    // }

    // std::cout<<"Next Max Plane"<<std::endl;
    // for (auto hit : nextmaxPlaneHits){
    //   geo::View_t view = hit->View(); 
    //   std::cout<<(hit->WireID()).toString()<<" View: "<<view;
    //   if (view==geo::kU) { std::cout<<" kU"<<std::endl; };
    //   if (view==geo::kV) { std::cout<<" kV"<<std::endl; };
    //   if (view==geo::kW) { std::cout<<" kW"<<std::endl; };
    // }
    
    
    //std::cout<<"planes: "<<(maxPlaneHit->WireID()).toString()<<" and "
    //         <<(nextmaxPlaneHit->WireID()).toString() <<std::endl;
    
    //Build the 3D track
    pma::Track3D* pmatrack = fProjectionMatchingAlg.buildSegment(maxPlaneHits, nextmaxPlaneHits, ShowerStartPosition);

    //  pma::Track3D* pmatrack = fProjectionMatchingAlg.buildSegment(track1, track2, trackStart);
    if(!pmatrack){
      // TODO fix this and put exception back in
      mf::LogError("ShowerTrackFinder") << "No PMA track made " << std::endl;
      //throw cet::exception("ShowerTrackFinderEMShower") << "No PMA track made ";
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
      
      TVector3 dir = spt - spt1;
      
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

    //Actually make the thing.
    recob::Track track = recob::Track(recob::TrackTrajectory(std::move(xyz), std::move(pxpypz), std::move(flags), false),
				      util::kBogusI, util::kBogusF, util::kBogusI, recob::tracking::SMatrixSym55(), 
				      recob::tracking::SMatrixSym55(), pfparticle.key());

    //Not correct need a different way.
    // const art::PtrMaker<recob::Track> makeTrackPtr(Event);
    // art::Ptr<recob::Track> initialtrack(makeTrackPtr(pfparticle.key()));

    // //Save the track assns
    // for (auto const& PlaneIter: trackhits)
    //   for(auto const& hitPtr: PlaneIter)
    // 	assnhit->addSingle(initialtrack, hitPtr);

    // Event.put(std::move(assnhit));
    ShowerPropHolder.SetInitialTrack(track); 

    std::cout<<"Track made"<<std::endl;

    // Shower dEdx calculation
    // To be moved when we have figured out how to transfer the hits associated with the track
    // TODO: when moved out in some protection to make sure there is an initial track    

    std::cout << "hello world dEdx (in track finder)" << std::endl;

    std::vector<double> dEdxVec;
    //std::map<geo::PlaneID, std::vector<art::Ptr<recob::Hit> > > trackhits;
    std::vector<std::vector<art::Ptr<recob::Hit> > > trackHits;
    TVector3 showerDir = ShowerPropHolder.GetShowerDirection();
    
    unsigned int numPlanes = fGeom->Nplanes();

    //dEdx.resize(numPlanes);
    trackHits.resize(numPlanes);

    //std::cout<<"Test: numPlanes: "<<numPlanes<<std::endl;

    // TODO replace trackhits look with loop over associated hits with track
    for(auto const& plane : plane_trackhits){
      std::cout<<(plane.first).toString()<<std::endl;
      std::vector<art::Ptr<recob::Hit> > trackhits = plane.second;
      for (unsigned int hitIt=0; hitIt<trackhits.size(); ++hitIt) {
	art::Ptr<recob::Hit> hit = trackhits.at(hitIt);
	geo::PlaneID hitWire = hit->WireID();
	geo::TPCID TPC = hitWire.asTPCID();
	if (TPC==vtxTPC){
	  (trackHits.at(hit->View())).push_back(hit);
	}
      }
    }

    int bestPlane = -999;
    int bestHitsPlane = 0;
    int bestPlaneHits = 0;
    double minPitch = 999;
    
    std::vector<art::Ptr<recob::Hit> > trackPlaneHits;
    for (unsigned int plane=0; plane<numPlanes; ++plane) {
      trackPlaneHits = trackHits.at(plane);
      std::cout<<"Plane "<<plane<<" with trackhits "<<trackPlaneHits.size()<<std::endl; 
      if (trackPlaneHits.size()){
	double fdEdx = -999;
	double totQ  = 0;
	double avgT  = 0;
	double pitch = 0;
	double wirepitch = fGeom->WirePitch(trackPlaneHits.at(0)->WireID().planeID());
	//TODO: check the numbers below
	double angleToVert = fGeom->WireAngleToVertical(fGeom->Plane(plane).View(),trackPlaneHits[0]->WireID().planeID()) - 0.5*TMath::Pi();
	double cosgamma = std::abs(sin(angleToVert)*showerDir.Y()+cos(angleToVert)*showerDir.Z());
	if (cosgamma>0) pitch = wirepitch/cosgamma;
	if (pitch){
	  //std::cout<<"pitch is good: "<<pitch<<std::endl;
	  
	  
	  int nhits = 0;
	  std::vector<float> vQ;
	  int w0 = trackPlaneHits.at(0)->WireID().Wire;
	  std::cout<<"Wire 0: "<<w0<<std::endl;
	  for (auto const& hit: trackPlaneHits){
	    int w1 = hit->WireID().Wire;	    
	    if (fMissFirstPoint){
	      if (w0==w1){
		continue;
	      }
	    }
	    //std::cout<<w1<<std::endl;
	    if (std::abs((w1-w0)*pitch)<fdEdxTrackLength){ //TODO: put back in
	      //if (std::abs((w1-w0)*pitch)<fdEdxTrackLength){ 
	      //should we be sorting from the vertex or first hit?
	      //std::cout<<"Wire: "<<w1<<" with dQdx: "<<hit->Integral()/pitch<<" and peak time: "<<hit->PeakTime()<<std::endl;
	      vQ.push_back(hit->Integral());
	      totQ += hit->Integral();
	      avgT+= hit->PeakTime();
	      ++nhits;
	    }
	  }
	  
	  if (totQ) {
	    if (pitch<minPitch){
	      //std::cout<<"New best plane: "<<plane<<" with pitch: "<<pitch<<std::endl;
	      minPitch  = pitch;
	      bestPlane = plane;
	    }

	    double dQdx = TMath::Median(vQ.size(), &vQ[0])/pitch;
	    //int expHits = (int)fdEdxTrackLength/pitch;	    
	    fdEdx = fCalorimetryAlg.dEdx_AREA(dQdx, avgT/nhits, trackPlaneHits[0]->WireID().Plane);
	    
	    //std::cout<<"dQdx: "<<dQdx<<" dEdx: "<<fdEdx<<" avgT/nhits: "<< avgT/nhits <<" numHits: "<<nhits<<" expected: "<<expHits<<" pitch: "<<pitch<<" plane: "<<trackPlaneHits[0]->WireID().Plane<<std::endl;

	    if (isinf(fdEdx)) { //TODO add error message logger
	      std::cout<<"Its inf"<<std::endl;
	      fdEdx=-999;
	    };
	   
	    if (nhits > bestPlaneHits || ((nhits==bestPlaneHits) && (pitch<minPitch))){
	       bestHitsPlane = plane;
	       bestPlaneHits = nhits;
	    }
	  }
	  //std::cout<<"dEdx: "<<fdEdx<<std::endl;
	  dEdxVec.push_back(fdEdx);  
	} else { // if not (pitch)
	  dEdxVec.push_back(-999);
	};
      } else { // if not (trackPlaneHits.size())
	dEdxVec.push_back(-999);
      }
      trackPlaneHits.clear();
    } //end loop over planes 

    
    ShowerPropHolder.SetShowerdEdx(dEdxVec);

    if (fMaxHitPlane){
      bestPlane=bestHitsPlane;
    }

    if (bestPlane==999){
      //mf::LogError("ShowerTrackFinder") << "No best plane set "<< std::endl;
      throw cet::exception("ShowerTrackFinderEMShower") << "No best plane set";
      return 1;
    } else {
      ShowerPropHolder.SetBestPlane(bestPlane);
      //ShowerPropHolder.SetBestPlane(bestHitsPlane);
    }
    
    std::cout<<"Best Plane: "<<bestPlane<<" and plane with most hits: "<<bestHitsPlane<<std::endl;


    std::cout <<"#########################################\n"<<
      "track finder done\n" <<"#########################################\n"<< std::endl;
    
    return 0;
  }
  
  //Function to calculate the what are the initial tracks hits. Adapted from EMShower FindInitialTrackHits
  std::vector<art::Ptr<recob::Hit> > ShowerTrackFinder::FindInitialTrackHits(std::vector<art::Ptr<recob::Hit> >& hits, TVector3& ShowerStartPosition, TVector3& ShowerDirection){
    
    std::vector<art::Ptr<recob::Hit> > trackHits;
    /*
    std::cout<<"ShowerHits: "<<hits.size()<<std::endl;
    for (auto hit : hits){
      TVector2 hitcoord = { (double) hit->WireID().Wire, hit->PeakTime()};
      std::cout<<"hit  : "<<hitcoord.X()<<" "<<hitcoord.Y()<<std::endl;
    }
    */
    double parm[2];
    int fitok = 0;
    std::vector<double> wfit;
    std::vector<double> tfit;
    std::vector<double> cfit;
    
    //std::cout<<"fNfitpass "<<fNfitpass<<std::endl;

    for (size_t i = 0; i<fNfitpass; ++i){
      
      //if (i>0) std::cout<<"loop: i "<<i<<" and fToler[i-1]: "<< fToler[i-1]<<" and fNfithits[i]+1 "<<fNfithits[i]+1<<std::endl;
      
      // Fit a straight line through hits
      unsigned int nhits = 0;
      for (auto &hit: hits){
	
	//Not sure I am a fan of doing things in wire tick space. What if id doesn't not iterate properly or the 
	//two planes in each TPC are not symmetric. 
	TVector2 coord = fSBNShowerAlg.HitCoordinates(hit);
	
	if (i==0||(std::abs((coord.Y()-(parm[0]+coord.X()*parm[1]))*cos(atan(parm[1])))<fToler[i-1])||fitok==1){
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
    /*
    std::cout<<"TrackHits: "<<trackHits.size()<<std::endl;
    for (auto hit : trackHits){
      TVector2 hitcoord = { (double) hit->WireID().Wire, hit->PeakTime()};
      std::cout<<"hit  : "<<hitcoord.X()<<" "<<hitcoord.Y()<<std::endl;
    }
    */
    return trackHits;
  }

  //Stolen from EMShowerAlg, a regression fitting function 
  Int_t ShowerTrackFinder::WeightedFit(const Int_t n, const Double_t *x, const Double_t *y, const Double_t *w,  Double_t *parm){

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


}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackFinder)
