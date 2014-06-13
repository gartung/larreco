#ifndef CORNERFINDER_H
#define CORNERFINDER_H
/*!
 * Title:   CornerFinder class
 * Author:  wketchum@lanl.gov
 * Inputs:  recob::Wire (calibrated)
 * Outputs: recob::EndPoint2D
 *
 * Description:
 * This module, is designed to produce EndPoint2D's using the CornerFinderAlg.
 * It's inlikely this module would actually be used on its own for vertex finding, but
 * is meant to be used as an input into a different module, or as a test piece.
 */


//Basic includes
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

//Framework includes
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

//LArSoft includes
#include "RecoBase/Wire.h"
#include "RecoBase/EndPoint2D.h"
#include "RecoAlg/CornerFinderAlg.h"
#include "Geometry/Geometry.h"

//header info for CornerFinder class

namespace vertex {
   
  class CornerFinder :  public art::EDProducer {
    
  public:
    
    explicit CornerFinder(fhicl::ParameterSet const& pset); 
    virtual ~CornerFinder();        
    void reconfigure(fhicl::ParameterSet const& p);

    
    void produce(art::Event& evt);
    
  private:
    corner::CornerFinderAlg           fCornerAlg;
    art::ServiceHandle<geo::Geometry> fGeometryHandle;
    std::string                       fCalDataModuleLabel;
    bool                              fPrintFlag;

    void printEndpoints(std::vector<recob::EndPoint2D> const& corner_vector);

    struct scoreCompare{
      bool operator() (recob::EndPoint2D a, recob::EndPoint2D b)
      { return a.Strength() > b.Strength(); }
    } fScoreCompare;

  }; //class CornerFinder
  
  

  //-----------------------------------------------------------------------------
  CornerFinder::CornerFinder(fhicl::ParameterSet const& pset) {  
    this->reconfigure(pset);    
    produces< std::vector<recob::EndPoint2D> >();
  }
  
  //-----------------------------------------------------------------------------
  CornerFinder::~CornerFinder(){}

  //---------------------------------------------------------------------------
  void CornerFinder::reconfigure(fhicl::ParameterSet const& pset) {
    fCalDataModuleLabel                = pset.get< std::string  >("CalDataModuleLabel");
    fPrintFlag                         = pset.get< bool >("PrintFlag",false);
    fhicl::ParameterSet pset_CornerAlg = pset.get< fhicl::ParameterSet >("CornerAlgParamSet");

    fCornerAlg.setSparsify(pset_CornerAlg.get<bool>("Sparsify",true));
    fCornerAlg.setSparseReserveSize(pset_CornerAlg.get<unsigned int>("SparseReserveSize",1e6));

    fCornerAlg.setSmoothKernelParams( pset_CornerAlg.get<unsigned int>("SmoothNeighborhoodX",4),
				      pset_CornerAlg.get<unsigned int>("SmoothNeighborhoodY",4),
				      pset_CornerAlg.get<float>("SmoothWidthX",1.),
				      pset_CornerAlg.get<float>("SmoothWidthY",1.));

    if( (pset_CornerAlg.get<std::string>("ImageTransformAlg","Threshold")).compare("Threshold")==0 )
      fCornerAlg.setImageAlg_Threshold( pset_CornerAlg.get<float>("InductionADCThreshold",3),
					pset_CornerAlg.get<float>("ColectionADCThreshold",6),
					pset_CornerAlg.get<float>("InductionAboveThresholdValue",-1),
					pset_CornerAlg.get<float>("CollectionAboveThresholdValue",-1),
					pset_CornerAlg.get<float>("InductionBelowThresholdValue",0),
					pset_CornerAlg.get<float>("CollectionBelowThresholdValue",0));
    else if( (pset_CornerAlg.get<std::string>("ImageTransformAlg")).compare("Nothing")==0 )
      fCornerAlg.setImageAlg_Nothing();

    fCornerAlg.setStructureTensorWindow( pset_CornerAlg.get<std::string>("STWindowType","Gaussian"),
					 pset_CornerAlg.get<unsigned int>("STWindowNeighborhood",2));

    fCornerAlg.setCornerScoreOptions( pset_CornerAlg.get<std::string>("CornerScoreMethod","Harris"),
				      pset_CornerAlg.get<double>("CornerScoreMin",1e10),
				      pset_CornerAlg.get<unsigned int>("LocalMaxNeighborhood",2));

  }

  //-----------------------------------------------------------------------------
  void CornerFinder::produce(art::Event& evt){
  
    //We need do very little here, as it's all handled by the corner finder.
    const bool DEBUG_TEST = fPrintFlag; //turn on/off some messages

    //We need to grab out the wires.
    art::Handle< std::vector<recob::Wire> > wireHandle;
    evt.getByLabel(fCalDataModuleLabel,wireHandle);
    std::vector<recob::Wire> const& wireVec(*wireHandle);

    geo::Geometry const& my_geometry(*fGeometryHandle);

    //First, have it process the wires.
    fCornerAlg.GrabWires(wireVec,my_geometry);
    
    //now, make a vector of recob::EndPoint2Ds, and hand that to CornerAlg to fill out
    std::unique_ptr< std::vector<recob::EndPoint2D> > corner_vector(new std::vector<recob::EndPoint2D>);
    fCornerAlg.GetFeaturePoints(*corner_vector,my_geometry);

    std::sort (corner_vector->begin(),corner_vector->end(),fScoreCompare);
    
    //mf::LogInfo("CornerFinderModule") 
    std::cout  
      << "CornerFinderAlg finished, and returned " 
      << corner_vector->size() << " endpoints." << std::endl;
    
    if(DEBUG_TEST) printEndpoints(*corner_vector);
    
    //and now, put this on the event.
    evt.put(std::move(corner_vector));

    //Done!

  } // end of produce


  //-----------------------------------------------------------------------------
  void CornerFinder::printEndpoints(std::vector<recob::EndPoint2D> const& corner_vector){

    for(auto iter=corner_vector.begin(); iter!=corner_vector.end(); iter++){
      geo::WireID wid = iter->WireID();
      //mf::LogVerbatim("CornerFinderModule") 
      std::cout	
	<< "Endpoint found: (plane,wire,time,strength,charge)=(" 
	<< wid.Plane << "," 
	<< wid.Wire << "," 
	<< iter->DriftTime() << ","
	<< iter->Strength() << ","
	<< iter->Charge() << ","
	<< ")" << std::endl;
    }

  }


  DEFINE_ART_MODULE(CornerFinder)

}
#endif // CORNERFINDER_H
