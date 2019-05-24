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
  };
  
  
  ShowerTrackFinder::ShowerTrackFinder(const fhicl::ParameterSet& pset)
    : fSBNShowerAlg(pset.get<fhicl::ParameterSet>("SBNShowerAlg")),
      fProjectionMatchingAlg(pset.get<fhicl::ParameterSet>("ProjectionMatchingAlg"))

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

    //This is all based on the shower vertex being known. If it is not lets not do the track 
    if(!ShowerPropHolder.CheckShowerStartPosition()){return 0;}
    if(!ShowerPropHolder.CheckShowerDirection()){return 0;}

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
    
    //Get the hit association 
    art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleModuleLabel);

    std::map<geo::View_t, std::vector<art::Ptr<recob::Hit> > > clusters_view;

    //Loop over the clusters in the plane and get the hits 
    for(auto const& cluster: clusters){
      
      //Get the hits 
      std::vector<art::Ptr<recob::Hit> > hits = fmhc.at(cluster.key());
      
      //Get the view.
      geo::View_t view = cluster->View();
 
      clusters_view[view].insert(clusters_view[view].end(),hits.begin(),hits.end());
    }


    std::vector<std::vector<art::Ptr<recob::Hit> > > trackhits;
    //Loop over the clusters and order the hits and get the initial track hits in that plane
    for(auto const& cluster: clusters_view){

      //Get the hits 
      std::vector<art::Ptr<recob::Hit> > hits = cluster.second;
      
      //Order the hits 
      fSBNShowerAlg.OrderShowerHits(hits,ShowerStartPosition,ShowerDirection);

      //Find the initial track hits 
      std::vector<art::Ptr<recob::Hit> > trackhits_plane = FindInitialTrackHits(hits,ShowerStartPosition,ShowerDirection);
      trackhits.push_back(trackhits_plane);
    }

    //Decide which plane to use 
    int maxplane = -1;
    int maxhits  = -1;
    for(unsigned int plane=0; plane<trackhits.size(); ++plane){
      if((int) trackhits[plane].size() > maxhits){
	maxplane = plane;
      }
    }
     
    if(maxplane == -1 || maxhits == -1){
            throw cet::exception("ShowerTrackFinderEMShower") << "Max Plane not set";
	    return 1;
    }


    int nextmaxplane = -1;
    maxhits      = -1;
    for(unsigned int plane=0; plane<trackhits.size(); ++plane){
      
      if((int) plane == maxplane){continue;}
      
      if((int) trackhits[plane].size() > maxhits){
	maxplane = plane;
      }
    }

    if(maxplane == -1 || maxhits == -1){
            throw cet::exception("ShowerTrackFinderEMShower") << "Next Max Plane not set";
	    return 1;
    }

    //Build the 3D track
    pma::Track3D* pmatrack = fProjectionMatchingAlg.buildSegment(trackhits[maxplane], trackhits[nextmaxplane]);

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
    return 0;
  }
  
  //Function to calculate the what are the initial tracks hits. Adapted from EMShower FindInitialTrackHits
  std::vector<art::Ptr<recob::Hit> > ShowerTrackFinder::FindInitialTrackHits(std::vector<art::Ptr<recob::Hit> >& hits, TVector3& ShowerStartPosition, TVector3& ShowerDirection){
    
    std::vector<art::Ptr<recob::Hit> > trackHits;

    double parm[2];
    int fitok = 0;
    std::vector<double> wfit;
    std::vector<double> tfit;
    std::vector<double> cfit;
    
    for (size_t i = 0; i<fNfitpass; ++i){
      
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

	  if(fApplyChargeWeight) cfit.push_back(hit->Integral());
	  else cfit.push_back(1.);
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


}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackFinder)
