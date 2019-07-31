//############################################################################
//### Name:        ShowerStartPositionCheater                              ###
//### Author:      Ed Tyley                                                ###
//### Date:        16.07.19                                                ###
//### Description: Cheating tool using truth for shower direction          ###
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

//LArSoft Includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/RecoAlg/SBNShowerAlg.h"
#include "larreco/RecoAlg/SBNShowerCheatingAlg.h"
#include "larreco/RecoAlg/MCRecoUtils/ShowerUtils.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "nusimdata/SimulationBase/MCTruth.h"

//C++ Includes
#include <iostream>
#include <cmath>

//Root Includes
#include "TVector3.h"
#include "TMath.h"
#include "TPrincipal.h"
#include "TVector.h"
#include "TTree.h"

namespace ShowerRecoTools {

	class ShowerStartPositionCheater:IShowerTool {

		public:

			ShowerStartPositionCheater(const fhicl::ParameterSet& pset);

			~ShowerStartPositionCheater();

			//Generic Direction Finder
			int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
					art::Event& Event,
					reco::shower::ShowerElementHolder& ShowerEleHolder
					) override;

		private:

			// Define standard art tool interface
			void configure(const fhicl::ParameterSet& pset) override;

			//Algorithm functions
			shower::SBNShowerAlg         fSBNShowerAlg;
			shower::SBNShowerCheatingAlg fSBNShowerCheatingAlg;

			//FCL
			art::InputTag fPFParticleModuleLabel;
	};


	ShowerStartPositionCheater::ShowerStartPositionCheater(const fhicl::ParameterSet& pset)
		: fSBNShowerAlg(pset.get<fhicl::ParameterSet>("SBNShowerAlg")),
		fSBNShowerCheatingAlg(pset.get<fhicl::ParameterSet>("SBNShowerCheatingAlg"))
	{
		configure(pset);
	}

	ShowerStartPositionCheater::~ShowerStartPositionCheater()
	{
	}

	void ShowerStartPositionCheater::configure(const fhicl::ParameterSet& pset)
	{
		fPFParticleModuleLabel  = pset.get<art::InputTag>("PFParticleModuleLabel","");
	}

	int ShowerStartPositionCheater::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle, art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){

		std::cout <<"#########################################\n"<<
			"hello world start position cheater\n" <<"#########################################\n"<< std::endl;

		//Could store these in the shower element holder and just calculate once?
		std::map<int,const simb::MCParticle*> trueParticles = fSBNShowerCheatingAlg.GetTrueParticleMap();
		std::map<int,std::vector<int> > showersMothers = fSBNShowerCheatingAlg.GetTrueChain(trueParticles);

		//Get the hits from the shower:
		art::Handle<std::vector<recob::PFParticle> > pfpHandle;
		if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
			throw cet::exception("ShowerLinearEnergy") << "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
			return 1;
		}

		//Get the clusters
		art::Handle<std::vector<recob::Cluster> > clusHandle;
		if (!Event.getByLabel(fPFParticleModuleLabel, clusHandle)){
			throw cet::exception("ShowerLinearEnergy") << "Could not get the pandora clusters. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
			return 1;
		}

		art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleModuleLabel);
		std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(pfparticle.key());

		//Get the hit association
		art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleModuleLabel);

		std::vector<art::Ptr<recob::Hit> > showerHits;
		for(auto const& cluster: clusters){

			//Get the hits
			std::vector<art::Ptr<recob::Hit> > hits = fmhc.at(cluster.key());
			showerHits.insert(showerHits.end(),hits.begin(),hits.end());
		}

		//Get the true particle from the shower
		std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain(showersMothers,showerHits,2);

		if(ShowerTrackInfo.first==-99999) {
			std::cout<<"True Shower Not Found"<<std::endl;
			return 1;
		}

		const simb::MCParticle* trueParticle = trueParticles[ShowerTrackInfo.first];
		TVector3 trueStartPos = {trueParticle->Vx(),trueParticle->Vy(),trueParticle->Vz()};

		TVector3 trueStartPosErr = {-999,-999,-999};
		ShowerEleHolder.SetElement(trueStartPos,trueStartPosErr,"ShowerStartPosition");

		ShowerEleHolder.SetElement(trueParticle,"TrueParticle");

		std::cout <<"#########################################\n"<<
			"start position cheater Done\n" <<"#########################################\n"<< std::endl;
		return 0;
	}
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerStartPositionCheater)
