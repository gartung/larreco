//////////////////////////////////////////////////
//
// GeneralCalorimetry module based on Calorimetry module
// but does not assume anything about which view is collection vs induction
//
//  This algorithm is designed to perform the calorimetric reconstruction 
//  of the 3D reconstructed tracks
//
// brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <vector>
#include <string>
#include <algorithm>

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/ArtDataHelper/TrackUtils.h" // lar::util::TrackPitchInView()

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "canvas/Persistency/Common/FindMany.h"
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

///calorimetry
namespace calo {
   
    class GeneralCalorimetry : public art::EDProducer {
    
    public:
    
      explicit GeneralCalorimetry(fhicl::ParameterSet const& pset); 
    
      void reconfigure(fhicl::ParameterSet const& pset);
      void produce(art::Event& evt);

    private:
        
      std::string    fTrackModuleLabel; ///< module creating the track objects and assns to hits
      double         fADCToElectrons;   ///< filled using the detinfo::DetectorPropertiesService service
      geo::View_t    fCollectionView;   ///< view of the collection plane
      unsigned int   fCollectionPlane;  ///< plane of the collection plane
      art::ServiceHandle<geo::Geometry> fGeo;

      CalorimetryAlg caloAlg;


    }; // class GeneralCalorimetry
  
}

//-------------------------------------------------
calo::GeneralCalorimetry::GeneralCalorimetry(fhicl::ParameterSet const& pset)
  : fCollectionView(geo::kUnknown)
  , fCollectionPlane(0)
  , caloAlg(pset.get< fhicl::ParameterSet >("CaloAlg"))
{
  auto const* dp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  fADCToElectrons = 1./dp->ElectronsToADC();

  // determine the view of the collection plane
  // just look at one cryostat, the first TPC and loop over those 
  // planes
  geo::TPCID FirstTPC(0, 0);
  for(unsigned int p = 0; p < fGeo->Nplanes(FirstTPC); ++p){
    geo::PlaneID planeid{ FirstTPC, p };
    if(fGeo->SignalType(planeid) == geo::kCollection){
      fCollectionView = fGeo->View(planeid);
      fCollectionPlane = p;
    }
  }

  this->reconfigure(pset);

  produces< std::vector<anab::Calorimetry>              >();
  produces< art::Assns<recob::Track, anab::Calorimetry> >();
}

//------------------------------------------------------------------------------------//
void calo::GeneralCalorimetry::reconfigure(fhicl::ParameterSet const& pset)
{
  fTrackModuleLabel = pset.get< std::string >("TrackModuleLabel");
  return;
}

//------------------------------------------------------------------------------------//
void calo::GeneralCalorimetry::produce(art::Event& evt)
{ 
  art::Handle< std::vector<recob::Track> > trackHandle;
  evt.getByLabel(fTrackModuleLabel, trackHandle);
  std::vector< art::Ptr<recob::Track> > tracks;
  art::fill_ptr_vector(tracks, trackHandle);
  
  auto const* dp = lar::providerFrom<detinfo::DetectorPropertiesService>();

  //create anab::Calorimetry objects and make association with recob::Track
  std::unique_ptr< std::vector<anab::Calorimetry> > calorimetrycol(new std::vector<anab::Calorimetry>);
  std::unique_ptr< art::Assns<recob::Track, anab::Calorimetry> > assn(new art::Assns<recob::Track, anab::Calorimetry>);

  art::FindMany<recob::Hit> fmht(trackHandle, evt, fTrackModuleLabel);

  // loop over the tracks
  for(size_t t = 0; t < tracks.size(); ++t){

    art::Ptr<recob::Track> trk = tracks.at(t);

    double viewPitch = 0.;
    try{
      viewPitch = lar::util::TrackPitchInView(*trk, fCollectionView);
    }
    catch( cet::exception &e){
      mf::LogWarning("GeneralCalorimetry") << "caught exception " 
					   << e 
					   << "\n pitch now set to 0";
    }
    
    // loop over the track trajectory to get the kinetic energy,
    // residual range and the dQdx
    float kineticEnergy = 0.;
    std::vector<float> vdEdx;
    std::vector<float> vresRange;
    std::vector<float> vdQdx;
    std::vector<float> deadwire; //residual range for dead wires

    // Check that the number of trajectory points and number of dQdx for the
    // collection view are the same
    if( trk->NumberTrajectoryPoints() != trk->NumberdQdx(fCollectionView) )
      throw cet::exception("GeneralCalorimetry") << "inconsistent number of track trajectory "
						 << " and dQdx points\n";
      
    if (trk->NumberTrajectoryPoints()>2){
      for(size_t p = 1; p < trk->NumberTrajectoryPoints()-1; ++p){	
	if (!trk->DQdxAtPoint(p, fCollectionView)) continue;
	vresRange.push_back(trk->Length(p));
	vdQdx.push_back(trk->DQdxAtPoint(p, fCollectionView));
	//vdEdx.push_back(vdQdx.back() * fADCToElectrons/util::kGeVToElectrons*1000);
	vdEdx.push_back(caloAlg.dEdx_AMP(vdQdx.back(),
					 dp->ConvertXToTicks(trk->LocationAtPoint(p)[0],fCollectionPlane,0,0),
					 fCollectionPlane));
	kineticEnergy += vdEdx.back(); // \todo should this be converted from electrons to energy?
	std::cout<<vresRange.back()<<" "<<vdQdx.back()<<" "<<vdEdx.back()<<std::endl;
      }
      
      geo::PlaneID planeID(0,0,fCollectionPlane);
      calorimetrycol->push_back(anab::Calorimetry(kineticEnergy,
						  vdEdx,
						  vdQdx,
						  vresRange,
						  deadwire,
						  trk->Length(),
						  viewPitch,
						  planeID));
      
      util::CreateAssn(*this, evt, *calorimetrycol, trk, *assn);
      
    }
  }// end of loop over all tracks 

  if (calorimetrycol->size()){
    evt.put(std::move(calorimetrycol));
    evt.put(std::move(assn));
  }

  return;
}

namespace calo{

  DEFINE_ART_MODULE(GeneralCalorimetry)
  
} // end namespace 

