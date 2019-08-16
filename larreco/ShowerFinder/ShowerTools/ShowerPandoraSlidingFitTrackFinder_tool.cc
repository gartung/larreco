//#############################################################################
//### Name:        ShowerPandoraSlidingFitTrackFinder                       ###
//### Author:      Dominic Barker (dominic.barker@sheffield.ac.uk           ###
//### Date:        30.07.19                                                 ###
//### Description: Tool for finding the initial shower track using the      ###
//###              pandora sliding fit calculation. This method is derived  ###
//###              from the PandoraTrackCreationModule.cc                   ###
//#############################################################################
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

//LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

//C++ Includes
#include <iostream>
#include <math.h>

//Root Includes
#include "TVector3.h"

namespace ShowerRecoTools{

  class ShowerPandoraSlidingFitTrackFinder:IShowerTool {
    public:

      ShowerPandoraSlidingFitTrackFinder(const fhicl::ParameterSet& pset);

      ~ShowerPandoraSlidingFitTrackFinder();

      //Generic Track Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder) override;
    private:

      void InitialiseProducers() override;

      //Function to add the assoctions
      int AddAssociations(art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder) override;


      // Define standard art tool interface.
      art::ServiceHandle<geo::Geometry> fGeom;
      float fSlidingFitHalfWindow;
      float fMinTrajectoryPoints;
  };


  ShowerPandoraSlidingFitTrackFinder::ShowerPandoraSlidingFitTrackFinder(const fhicl::ParameterSet& pset)
  {
    fSlidingFitHalfWindow = pset.get<float> ("SlidingFitHalfWindow");
    fMinTrajectoryPoints  = pset.get<float> ("MinTrajectoryPoints");
  }

  ShowerPandoraSlidingFitTrackFinder::~ShowerPandoraSlidingFitTrackFinder()
  {
  }

  void ShowerPandoraSlidingFitTrackFinder::InitialiseProducers(){
    if(producerPtr == NULL){
      mf::LogWarning("ShowerPandoraSlidingFitTrackFinder") << "The producer ptr has not been set" << std::endl;
      return;
    }

    InitialiseProduct<std::vector<recob::Track> >("InitialTrack");
    InitialiseProduct<art::Assns<recob::Shower, recob::Track > >("ShowerTrackAssn");
    InitialiseProduct<art::Assns<recob::Track, recob::SpacePoint > >("ShowerTrackSpacepointAssn");

  }


  //This whole idea is stolen from PandoraTrackCreationModule so credit goes to the Pandora guys.
  int ShowerPandoraSlidingFitTrackFinder::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event,
      reco::shower::ShowerElementHolder& ShowerEleHolder
      ){
    //This is all based on the shower vertex being known. If it is not lets not do the track
    if(!ShowerEleHolder.CheckElement("ShowerStartPosition")){
      mf::LogError("ShowerPandoraSlidingFitTrackFinder") << "Start position not set, returning "<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement("ShowerDirection")){
      mf::LogError("ShowerPandoraSlidingFitTrackFinder") << "Direction not set, returning "<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement("InitialTrackSpacePoints")){
      mf::LogError("ShowerPandoraSlidingFitTrackFinder") << "Initial Spacepoints not set, returning "<< std::endl;
      return 1;
    }

    TVector3 ShowerStartPosition = {-999,-999,-999};
    ShowerEleHolder.GetElement("ShowerStartPosition",ShowerStartPosition);

    TVector3 ShowerDirection     = {-999,-999,-999};
    ShowerEleHolder.GetElement("ShowerDirection",ShowerDirection);

    std::vector<art::Ptr<recob::SpacePoint> > spacepoints;
    ShowerEleHolder.GetElement("InitialTrackSpacePoints",spacepoints);

    const pandora::CartesianVector vertexPosition(ShowerStartPosition.X(), ShowerStartPosition.Y(),
        ShowerStartPosition.Z());

    pandora::CartesianPointVector cartesianPointVector;
    for (const art::Ptr<recob::SpacePoint> spacePoint : spacepoints)
      cartesianPointVector.emplace_back(pandora::CartesianVector(spacePoint->XYZ()[0],
            spacePoint->XYZ()[1], spacePoint->XYZ()[2]));

    lar_content::LArTrackStateVector trackStateVector;
    pandora::IntVector indexVector;
    try{
      lar_content::LArPfoHelper::GetSlidingFitTrajectory(cartesianPointVector, vertexPosition,
          fSlidingFitHalfWindow, fGeom->WirePitch(geo::kW), trackStateVector, &indexVector);
    }
    catch (const pandora::StatusCodeException &){
      mf::LogWarning("ShowerPandoraSlidingFitTrackFinder") << "Unable to extract sliding fit trajectory" << std::endl;
      return 1;
    }
    if (trackStateVector.size() < fMinTrajectoryPoints){
      mf::LogWarning("ShowerPandoraSlidingFitTrackFinder") << "Insufficient input trajectory points to build track: " << trackStateVector.size();
      return 1;
    }

    if (trackStateVector.empty())
      throw cet::exception("ShowerPandoraSlidingFitTrackFinder") << "BuildTrack - No input trajectory points provided " << std::endl;

    recob::tracking::Positions_t xyz;
    recob::tracking::Momenta_t pxpypz;
    recob::TrackTrajectory::Flags_t flags;

    for (const lar_content::LArTrackState &trackState : trackStateVector)
    {
      xyz.emplace_back(recob::tracking::Point_t(trackState.GetPosition().GetX(),
            trackState.GetPosition().GetY(), trackState.GetPosition().GetZ()));
      pxpypz.emplace_back(recob::tracking::Vector_t(trackState.GetDirection().GetX(),
            trackState.GetDirection().GetY(), trackState.GetDirection().GetZ()));
      // Set flag NoPoint if point has bogus coordinates, otherwise use clean flag set
      if (std::fabs(trackState.GetPosition().GetX()-util::kBogusF)<std::numeric_limits<float>::epsilon() &&
          std::fabs(trackState.GetPosition().GetY()-util::kBogusF)<std::numeric_limits<float>::epsilon() &&
          std::fabs(trackState.GetPosition().GetZ()-util::kBogusF)<std::numeric_limits<float>::epsilon())
      {
        flags.emplace_back(recob::TrajectoryPointFlags(recob::TrajectoryPointFlags::InvalidHitIndex,
              recob::TrajectoryPointFlagTraits::NoPoint));
      } else {
        flags.emplace_back(recob::TrajectoryPointFlags());
      }
    }

    // note from gc: eventually we should produce a TrackTrajectory, not a Track with empty covariance matrix and bogus chi2, etc.
    recob::Track InitialTrack  =  recob::Track(recob::TrackTrajectory(std::move(xyz),
          std::move(pxpypz), std::move(flags), false),
        util::kBogusI, util::kBogusF, util::kBogusI, recob::tracking::SMatrixSym55(),
        recob::tracking::SMatrixSym55(), pfparticle.key());

    ShowerEleHolder.SetElement(InitialTrack,"InitialTrack");

    TVector3 Start = {InitialTrack.Start().X(), InitialTrack.Start().Y(), InitialTrack.Start().Z()};
    TVector3 End   = {InitialTrack.End().X(), InitialTrack.End().Y(),InitialTrack.End().Z()};
    float tracklength = (Start-End).Mag();

    ShowerEleHolder.SetElement(tracklength,"InitialTrackLength");


    return 0;
  }


  int ShowerPandoraSlidingFitTrackFinder::AddAssociations(art::Event& Event,
      reco::shower::ShowerElementHolder& ShowerEleHolder
      ){

    //Check the track has been set
    if(!ShowerEleHolder.CheckElement("InitialTrack")){
      mf::LogError("ShowerPandoraSlidingFitTrackFinderAddAssn") << "Track not set so the assocation can not be made  "<< std::endl;
      return 1;
    }

    //Get the size of the ptr as it is.
    int trackptrsize = GetVectorPtrSize("InitialTrack");

    const art::Ptr<recob::Track> trackptr = GetProducedElementPtr<recob::Track>("InitialTrack",
        ShowerEleHolder,trackptrsize-1);
    const art::Ptr<recob::Shower> showerptr = GetProducedElementPtr<recob::Shower>("shower",
        ShowerEleHolder);

    AddSingle<art::Assns<recob::Shower, recob::Track> >(showerptr,trackptr,"ShowerTrackAssn");

    std::vector<art::Ptr<recob::SpacePoint> >  TrackSpacepoints;
    ShowerEleHolder.GetElement("InitialTrackSpacePoints",TrackSpacepoints);

    for(auto const& TrackSpacepoint: TrackSpacepoints){
      AddSingle<art::Assns<recob::Track, recob::SpacePoint> >
        (trackptr,TrackSpacepoint,"ShowerTrackSpacepointAssn");
    }

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerPandoraSlidingFitTrackFinder)
