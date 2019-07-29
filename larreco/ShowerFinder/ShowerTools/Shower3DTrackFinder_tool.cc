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

	class Shower3DTrackFinder:IShowerTool {
		public:

			Shower3DTrackFinder(const fhicl::ParameterSet& pset);

			~Shower3DTrackFinder();

			//Generic Track Finder
			int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
					art::Event& Event,
					reco::shower::ShowerElementHolder& ShowerEleHolder
					) override;

		private:

			void configure(const fhicl::ParameterSet& pset) override;

			void InitialiseProducers() override;

			int AddAssociations(art::Event& Event,
					reco::shower::ShowerElementHolder& ShowerEleHolder) override;


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
		fPFParticleModuleLabel = pset.get<art::InputTag>             ("PFParticleModuleLabel");

		fMaxProjectionDist     = pset.get<float> ("MaxProjectionDist");
		fMaxPerpendicularDist  = pset.get<float> ("MaxPerpendicularDist");
		fForwardHitsOnly       = pset.get<bool>  ("ForwardHitsOnly");
		fDebugEVD              = pset.get<bool>                      ("DebugEVD");  
	}

	void Shower3DTrackFinder::InitialiseProducers(){
		if(producerPtr == NULL){
			mf::LogWarning("Shower3DTrackFinder") << "The producer ptr has not been set" << std::endl;
			return;
		}

		InitialiseProduct<std::vector<recob::Track> >                   ("InitialTrack");
		InitialiseProduct<art::Assns<recob::Shower, recob::Track > >    ("ShowerTrackAssn");
		InitialiseProduct<art::Assns<recob::Track, recob::Hit > >       ("ShowerTrackHitAssn");
		InitialiseProduct<art::Assns<recob::Track, recob::SpacePoint> > ("ShowerTrackSpacePointAssn");
	}


	int Shower3DTrackFinder::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
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
		art::FindOneP<recob::Hit> fohsp(spHandle, Event, fPFParticleModuleLabel);
		if(!fohsp.isValid()){
			throw cet::exception("Shower3DTrackFinderPosition") << "Spacepoint and hit association not valid. Stopping.";
			return 1;
		}


		// Get the SpacePoints
		std::vector<art::Ptr<recob::SpacePoint> > spacePoints = fmspp.at(pfparticle.key());

		//We cannot progress with no spacepoints.
		if(spacePoints.size() == 0){
			mf::LogError("Shower3DTrackFinder") << "No space points, returning "<< std::endl;
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
			mf::LogWarning("Shower3DTrackFinder") << "Not Enough Hits to make Track "<<std::endl;
			return 1;
		}

		//Build the 3D track
		pma::Track3D* pmatrack = fProjectionMatchingAlg.buildSegment(trackHits, ShowerStartPosition);

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

		//Actually make the thing.
		recob::Track track = recob::Track(recob::TrackTrajectory(std::move(xyz), std::move(pxpypz), std::move(flags), false),
				util::kBogusI, util::kBogusF, util::kBogusI, recob::tracking::SMatrixSym55(),
				recob::tracking::SMatrixSym55(), pfparticle.key());

    ShowerEleHolder.SetElement(track,"InitialTrack");
		ShowerEleHolder.SetElement(trackHits, "InitialTrackHits");
		ShowerEleHolder.SetElement(trackSpacePoints,"InitialTrackSpacePoints");

		if (fDebugEVD){
			std::cout<<"Doing DebugEVD"<<std::endl;
			fSBNShowerAlg.DebugEVD(pfparticle,Event,ShowerEleHolder);
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


	int Shower3DTrackFinder::AddAssociations(art::Event& Event,
			reco::shower::ShowerElementHolder& ShowerEleHolder
			){

		if(!ShowerEleHolder.CheckElement("InitialTrack")){
			mf::LogError("ShowerTrackFinderAddAssn") << "Track not set so the assocation can not be made  "<< std::endl;
			return 1;
		}


		int trackptrsize = GetVectorPtrSize("InitialTrack");

		const art::Ptr<recob::Track> trackptr = GetProducedElementPtr<recob::Track>("InitialTrack", ShowerEleHolder,trackptrsize-1);
		const art::Ptr<recob::Shower> showerptr = GetProducedElementPtr<recob::Shower>("shower", ShowerEleHolder);
		AddSingle<art::Assns<recob::Shower, recob::Track> >(showerptr,trackptr,"ShowerTrackAssn");

		std::vector<art::Ptr<recob::Hit> > TrackHits;
		ShowerEleHolder.GetElement("InitialTrackHits",TrackHits);
		for(auto const& TrackHit: TrackHits){
			AddSingle<art::Assns<recob::Track, recob::Hit> >(trackptr,TrackHit,"ShowerTrackHitAssn");
		}

		std::vector<art::Ptr<recob::SpacePoint> > TrackSpacePoints;
		ShowerEleHolder.GetElement("InitialTrackSpacePoints",TrackSpacePoints);
		for(auto const& TrackSpacePoint: TrackSpacePoints){
			AddSingle<art::Assns<recob::Track, recob::SpacePoint> >(trackptr,TrackSpacePoint,"ShowerTrackSpacePointAssn");
		}
		return 0;
	}
}


DEFINE_ART_CLASS_TOOL(ShowerRecoTools::Shower3DTrackFinder)
