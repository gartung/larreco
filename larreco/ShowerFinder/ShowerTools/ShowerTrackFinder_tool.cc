//############################################################################
//### Name:        ShowerTrackFinder                                       ###
//### Author:      Dominic Barker (dominic.barker@sheffield.ac.uk          ###
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
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft Includes 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
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

  class ShowerTrackFinder:IShowerTool {
  public:

    ShowerTrackFinder(const fhicl::ParameterSet& pset);
    
    ~ShowerTrackFinder(); 
    
    //Generic Track Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
			 art::Event& Event,
			 reco::shower::ShowerElementHolder& ShowerEleHolder
			 ) override;

  private:  
    
    void configure(const fhicl::ParameterSet& pset) override;

    void InitialiseProducers() override;

    //Function to add the assoctions
    int AddAssociations(art::Event& Event,
			reco::shower::ShowerElementHolder& ShowerEleHolder) override;
    


    std::vector<art::Ptr<recob::Hit> > FindInitialTrackHits(std::vector<art::Ptr<recob::Hit> >& hits, TVector3& ShowerStartPosition, TVector3& ShowerDirection);

    Int_t WeightedFit(const Int_t n, const Double_t *x, const Double_t *y, const Double_t *w,  Double_t *parm);

    // Define standard art tool interface

    shower::SBNShowerAlg       fSBNShowerAlg;
    pma::ProjectionMatchingAlg fProjectionMatchingAlg;
    unsigned int               fNfitpass;
    std::vector<unsigned int>  fNfithits;
    std::vector<double>        fToler;
    bool                       fApplyChargeWeight;
    art::InputTag              fPFParticleModuleLabel;
    detinfo::DetectorProperties const* fDetProp;
    art::ServiceHandle<geo::Geometry> fGeom;
  };
  
  
  ShowerTrackFinder::ShowerTrackFinder(const fhicl::ParameterSet& pset)
    : fSBNShowerAlg(pset.get<fhicl::ParameterSet>("SBNShowerAlg")),
      fProjectionMatchingAlg(pset.get<fhicl::ParameterSet>("ProjectionMatchingAlg")),
      fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())
  {
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

    if (fNfitpass!=fNfithits.size() ||
	fNfitpass!=fToler.size()) {
      throw art::Exception(art::errors::Configuration)
	<< "ShowerTrackFinderEMShower: fNfithits and fToler need to have size fNfitpass";
    }
  }

  void ShowerTrackFinder::InitialiseProducers(){
    if(producerPtr == NULL){
      mf::LogWarning("ShowerStartPosition") << "The producer ptr has not been set" << std::endl;
      return;
    }

    InitialiseProduct<std::vector<recob::Track> >("InitialTrack"); 
    InitialiseProduct<art::Assns<recob::Shower, recob::Track > >("ShowerTrackAssn");
    InitialiseProduct<art::Assns<recob::Track, recob::Hit > >("ShowerTrackHitAssn");


  }

  
  
  int ShowerTrackFinder::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
					  art::Event& Event,
					  reco::shower::ShowerElementHolder& ShowerEleHolder
					   ){
    std::cout <<"#########################################\n"<<
      "hello world track finder\n" <<"#########################################\n"<< std::endl;

    //This is all based on the shower vertex being known. If it is not lets not do the track 
    if(!ShowerEleHolder.CheckElement("ShowerStartPosition")){
      mf::LogError("ShowerTrackFinder") << "Start position not set, returning "<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement("ShowerDirection")){
      mf::LogError("ShowerTrackFinder") << "Direction not set, returning "<< std::endl;
      return 1;
    }
    
    TVector3 ShowerStartPosition = {-999,-999,-999};
    ShowerEleHolder.GetElement("ShowerStartPosition",ShowerStartPosition);
    
    TVector3 ShowerDirection     = {-999,-999,-999};
    ShowerEleHolder.GetElement("ShowerDirection",ShowerDirection);

    // Get the assocated pfParicle vertex PFParticles
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("ShowerTrackFinderEMShower") << "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      recob::Track track;
      ShowerEleHolder.SetElement(track,"InitialTrack");
      return 1;
    }
    
    //Get the clusters
    art::Handle<std::vector<recob::Cluster> > clusHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, clusHandle)){
      throw cet::exception("ShowerTrackFinderEMShower") << "Could not get the pandora clusters. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      recob::Track track;
      ShowerEleHolder.SetElement(track,"InitialTrack");
      return 1;
    }
    art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleModuleLabel);
    std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(pfparticle.key());
    

    if(clusters.size()<2){
      mf::LogError("ShowerTrackFinder") << "Not enough clusters: "<<clusters.size() << std::endl;
      // throw cet::exception("ShowerTrackFinderEMShower") << "Not enough clusters: "
      // <<clusters.size(<<" ");
      return 1;
    }

    //Get the hit association 
    art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleModuleLabel);
    std::map<geo::PlaneID, std::vector<art::Ptr<recob::Hit> > > plane_clusters;
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
      
      plane_trackhits[plane].insert(plane_trackhits[plane].end(),trackhits.begin(),trackhits.end());
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
      mf::LogError("ShowerTrackFinder") << "Max Plane not set " << std::endl;
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
      return 1;
    }
    
    std::vector<art::Ptr<recob::Hit> > maxPlaneHits = (plane_trackhits.at(maxplane));
    std::vector<art::Ptr<recob::Hit> > nextmaxPlaneHits = (plane_trackhits.at(nextmaxplane));

    if(maxPlaneHits.size() < 2 && nextmaxPlaneHits.size() < 2){
      mf::LogError("ShowerTrackFinder") << "Not enough hit to make track " << std::endl;
      return 0;
    }

    //Build the 3D track
    pma::Track3D* pmatrack = fProjectionMatchingAlg.buildSegment(maxPlaneHits, nextmaxPlaneHits, ShowerStartPosition);

    if(!pmatrack){
      mf::LogError("ShowerTrackFinder") << "Failed fit " << std::endl;
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
      mf::LogWarning("ShowerTrackFinder") << "Not Enough Spacepoints" << std::endl;
      recob::Track track;
      ShowerEleHolder.SetElement(track,"InitialTrack");
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

    ShowerEleHolder.SetElement(track,"InitialTrack"); 
    ShowerEleHolder.SetElement(TrackHits, "InitialTrackHits");

    std::cout <<"#########################################\n"<<
      "track finder done\n" <<"#########################################\n"<< std::endl;
    
    return 0;
  }
  
  //Function to calculate the what are the initial tracks hits. Adapted from EMShower FindInitialTrackHits
  std::vector<art::Ptr<recob::Hit> > ShowerTrackFinder::FindInitialTrackHits(std::vector<art::Ptr<recob::Hit> >& hits, TVector3& ShowerStartPosition, TVector3& ShowerDirection){

    //########################################
    // art::Ptr<recob::Hit> startHit = hits.front();
    
    // //Get the wireID 
    // const geo::WireID startWireID = startHit->WireID();
    
    // //Get the TPCID
    // const geo::TPCID tpcid = startWireID.asTPCID();  
    

    // //Get the plane
    // const geo::PlaneID planeid = startWireID.asPlaneID();

    // //Get the pitch
    // double pitch = fGeom->WirePitch(planeid);

    // //Get the projection vectors for the start position in 2D  
    // TVector2 Shower2DStartPosition = { 
    //   fGeom->WireCoordinate(ShowerStartPosition, startHit->WireID().planeID()), 
    //   fDetProp->ConvertXToTicks(ShowerStartPosition.X(),  startHit->WireID().planeID())
    // };
    
        
    // //Vector of the plane
    // TVector3 PlaneDirection = fGeom->Plane(planeid).GetIncreasingWireDirection();

    // double XTicksCoefficient    = 0.001 * fDetProp->DriftVelocity(fDetProp->Efield(), fDetProp->Temperature()) * fDetProp->SamplingRate();
    // if(fGeom->TPC(tpcid).DriftDirection() == geo::kNeg){ XTicksCoefficient*=-1;}

    // //get the shower 2D direction 
    // TVector2 Shower2DDirection = { 
    //   ShowerDirection.Dot(PlaneDirection)/pitch, 
    //   ShowerDirection.X()/XTicksCoefficient
    // };

    // Shower2DDirection = Shower2DDirection.Unit();

  //#######################################################
    
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
	// double proj = TMath::Abs(coord_hit.X()*Shower2DDirection.X() + coord_hit.Y()*Shower2DDirection.Y());
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

  int ShowerTrackFinder::AddAssociations(art::Event& Event,
					 reco::shower::ShowerElementHolder& ShowerEleHolder
					 ){
    
    //Check the track has been set 
    if(!ShowerEleHolder.CheckElement("InitialTrack")){
      mf::LogError("ShowerTrackFinderAddAssn") << "Track not set so the assocation can not be made  "<< std::endl;
      return 1;
    }

    //Get the size of the ptr as it is.
    int trackptrsize = GetVectorPtrSize("InitialTrack");

    const art::Ptr<recob::Track> trackptr = GetProducedElementPtr<recob::Track>("InitialTrack", ShowerEleHolder,trackptrsize-1);
    const art::Ptr<recob::Shower> showerptr = GetProducedElementPtr<recob::Shower>("shower", ShowerEleHolder);

    AddSingle<art::Assns<recob::Shower, recob::Track> >(showerptr,trackptr,"ShowerTrackAssn");

    std::vector<art::Ptr<recob::Hit> > TrackHits;
    ShowerEleHolder.GetElement("InitialTrackHits",TrackHits);   

    for(auto const& TrackHit: TrackHits){
      AddSingle<art::Assns<recob::Track, recob::Hit> >(trackptr,TrackHit,"ShowerTrackHitAssn");
    }
  
  return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackFinder)
