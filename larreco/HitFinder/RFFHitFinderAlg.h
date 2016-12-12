#ifndef RFFHITFINDERALG_H
#define RFFHITFINDERALG_H

/*!
 * Title:   RFFHitFinderAlg Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description: 
 * Class that runs the RFF HitFinder. Implements an RFFHitFitter, and takes 
 * the result and stores it in recob::Hit objects. 
 *
 * Input:  recob::Wire
 * Output: recob::Hit
*/

#include "fhiclcpp/ParameterSet.h"

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"

#include "RFFHitFitter.h"
#include <thread>

namespace hit{

  class RFFHitFinderAlg{

  static constexpr float SQRT_TWO_PI = 2.506628;

  public:
    RFFHitFinderAlg(fhicl::ParameterSet const&);
    ~RFFHitFinderAlg();

    void SetFitterParamsVectors(geo::Geometry const&);
    void Run(std::vector<recob::Wire> const&,
	     std::vector<recob::Hit>&,
	     geo::Geometry const&);

    typedef struct FitterInput{
      geo::SigType_t      sigtype;
      geo::WireID         wireid;
      raw::ChannelID_t    channel;
      geo::View_t         view;
      raw::TDCtick_t      roi_start;
      raw::TDCtick_t      roi_end;
      float               match_thresh;
      unsigned int        merge_multi;
      float               amp_thresh;
      const std::vector<float>* roi_dataptr;      
    } FitterInput_t;


    
  private:

    std::vector<float> fMatchThresholdVec;
    std::vector<unsigned int> fMergeMultiplicityVec;
    std::vector<float> fAmpThresholdVec;

    //void SetFitterParams(unsigned int);

    static void FitterWorker(void*);
    static void SendDataToWorkers(std::vector<recob::Wire> const&,
				  geo::Geometry const&,
				  void*,
				  std::vector<float> const&,
				  std::vector<unsigned int > const&,
				  std::vector<float> const&);

    
    static size_t EmplaceHit(std::vector<recob::Hit>&,
			     RFFHitFitter &,
			     float const&,
			     raw::TDCtick_t const&, raw::TDCtick_t const&,
			     geo::SigType_t const&, geo::WireID const&,
			     raw::ChannelID_t const&, geo::View_t const&);

    
    //RFFHitFitter fFitter;
    void *fContext;
    void *fPairNROI;
    void *fReceiver_Hits;
    void *fReceiver_fin;
    void* fController;

    size_t fNWorkers;
    std::vector<std::thread> fWorkers;
    
  };

}


#endif
