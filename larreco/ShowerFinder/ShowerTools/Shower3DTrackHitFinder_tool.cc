//############################################################################
//### Name:        Shower3DTrackHitFinder                                     ###
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

	class Shower3DTrackHitFinder:IShowerTool {
		public:

			Shower3DTrackHitFinder(const fhicl::ParameterSet& pset);

			~Shower3DTrackHitFinder();

			//Generic Track Finder
			int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
					art::Event& Event,
					reco::shower::ShowerElementHolder& ShowerEleHolder
					) override;

		private:

			void configure(const fhicl::ParameterSet& pset) override;

			void InitialiseProducers() override;


			std::vector<art::Ptr<recob::SpacePoint> > FindTrackSpacePoints(std::vector<art::Ptr<recob::SpacePoint> >& spacePoints, TVector3& showerStartPosition,TVector3& showerDirection);

			// Define standard art tool interface
			shower::SBNShowerAlg       fSBNShowerAlg;
			pma::ProjectionMatchingAlg fProjectionMatchingAlg;

			float fMaxProjectionDist;    // Maximum projection along shower direction, length of cylinder
			float fMaxPerpendicularDist; // Maximum perpendicular distance, radius of cylinder
			bool  fForwardHitsOnly;       // Only take hits downstream of shower vertex (projection>0)
			bool  fDebugEVD;

			// For 2D method for comparison
			art::InputTag              fPFParticleModuleLabel;
			detinfo::DetectorProperties const* fDetProp;
			art::ServiceHandle<geo::Geometry> fGeom;
	};


	Shower3DTrackHitFinder::Shower3DTrackHitFinder(const fhicl::ParameterSet& pset)
		: fSBNShowerAlg(pset.get<fhicl::ParameterSet>("SBNShowerAlg")),
			fProjectionMatchingAlg(pset.get<fhicl::ParameterSet>("ProjectionMatchingAlg")),
			fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())
	{
		configure(pset);
	}

	Shower3DTrackHitFinder::~Shower3DTrackHitFinder()
	{
	}

	void Shower3DTrackHitFinder::configure(const fhicl::ParameterSet& pset)
	{
		fPFParticleModuleLabel = pset.get<art::InputTag>             ("PFParticleModuleLabel");

		fMaxProjectionDist     = pset.get<float> ("MaxProjectionDist");
		fMaxPerpendicularDist  = pset.get<float> ("MaxPerpendicularDist");
		fForwardHitsOnly       = pset.get<bool>  ("ForwardHitsOnly");
		fDebugEVD              = pset.get<bool>                      ("DebugEVD");  
	}

	void Shower3DTrackHitFinder::InitialiseProducers(){
		if(producerPtr == NULL){
			mf::LogWarning("Shower3DTrackHitFinder") << "The producer ptr has not been set" << std::endl;
			return;
		}

		InitialiseProduct<std::vector<recob::Track> >                   ("InitialTrack");
		InitialiseProduct<art::Assns<recob::Shower, recob::Track > >    ("ShowerTrackAssn");
		InitialiseProduct<art::Assns<recob::Track, recob::Hit > >       ("ShowerTrackHitAssn");
		InitialiseProduct<art::Assns<recob::Track, recob::SpacePoint> > ("ShowerTrackSpacePointAssn");
	}


	int Shower3DTrackHitFinder::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
			art::Event& Event,
			reco::shower::ShowerElementHolder& ShowerEleHolder
			){
		std::cout <<"#########################################\n"<<
			"hello world track finder 3d\n" <<"#########################################\n"<< std::endl;

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

		// Get the assocated pfParicle Handle
		art::Handle<std::vector<recob::PFParticle> > pfpHandle;
		if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
			throw cet::exception("Shower3DTrackHitFinderEMShower") << "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
			return 1;
		}

		// Get the spacepoint - PFParticle assn
		art::FindManyP<recob::SpacePoint> fmspp(pfpHandle, Event, fPFParticleModuleLabel);
		if (!fmspp.isValid()){
			throw cet::exception("Shower3DTrackHitFinder") << "Trying to get the spacepoint and failed. Something is not configured correctly. Stopping ";
			return 1;
		}

		// Get the spacepoints
		art::Handle<std::vector<recob::SpacePoint> > spHandle;
		if (!Event.getByLabel(fPFParticleModuleLabel, spHandle)){
			throw cet::exception("Shower3DTrackHitFinder") << "Could not configure the spacepoint handle. Something is configured incorrectly. Stopping";
			return 1;
		}

		// Get the hits associated with the space points
		art::FindOneP<recob::Hit> fohsp(spHandle, Event, fPFParticleModuleLabel);
		if(!fohsp.isValid()){
			throw cet::exception("Shower3DTrackHitFinderPosition") << "Spacepoint and hit association not valid. Stopping.";
			return 1;
		}


		// Get the SpacePoints
		std::vector<art::Ptr<recob::SpacePoint> > spacePoints = fmspp.at(pfparticle.key());

		//We cannot progress with no spacepoints.
		if(spacePoints.size() == 0){
			mf::LogError("Shower3DTrackHitFinder") << "No space points, returning "<< std::endl;
			return 1;
		}

		//std::cout<<"Spacepoints: "<<spacePoints.size()<<std::endl;

		// Order the spacepoints
		fSBNShowerAlg.OrderShowerSpacePoints(spacePoints,ShowerStartPosition,ShowerDirection);

		// If we are selecting the track hits in 3D do it before spliting into planes
		std::vector<art::Ptr<recob::SpacePoint> > trackSpacePoints;
		trackSpacePoints = FindTrackSpacePoints(spacePoints,ShowerStartPosition,ShowerDirection);

		// Get the hits associated to the space points and seperate them by planes
		std::vector<art::Ptr<recob::Hit> > trackHits;
		for(auto const& spacePoint: trackSpacePoints){
			//Get the hits
			const art::Ptr<recob::Hit> hit = fohsp.at(spacePoint.key());
			trackHits.push_back(hit);
		}

		if (trackHits.size()<2){
			mf::LogWarning("Shower3DTrackHitFinder") << "Not Enough Hits to make Track "<<std::endl;
			return 1;
		}

		ShowerEleHolder.SetElement(trackHits, "InitialTrackHits");
		ShowerEleHolder.SetElement(trackSpacePoints,"InitialTrackSpacePoints");

		if (fDebugEVD){
			std::cout<<"Do DebugEVD"<<std::endl;
			//fSBNShowerAlg.DebugEVD(pfparticle,Event,ShowerPropHolder);
		}

		std::cout <<"#########################################\n"<<
			"track finder 3d done\n" <<"#########################################\n"<< std::endl;

		return 0;
	}

	std::vector<art::Ptr<recob::SpacePoint> > Shower3DTrackHitFinder::FindTrackSpacePoints(std::vector<art::Ptr<recob::SpacePoint> >& spacePoints, TVector3& showerStartPosition,TVector3& showerDirection){

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

}


DEFINE_ART_CLASS_TOOL(ShowerRecoTools::Shower3DTrackHitFinder)
