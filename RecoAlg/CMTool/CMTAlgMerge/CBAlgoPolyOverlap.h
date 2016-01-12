/**
 * \file CMalgoPolyOverlap.h
 *
 * \ingroup CMTool
 * 
 * \brief Class def header for a class CBAlgoPolyOverlap
 *
 * @author David Caratelli
 */

/** \addtogroup CMTool

    @{*/
#ifndef RECOTOOL_CBALGOPOLYOVERLAP_H
#define RECOTOOL_CBALGOPOLYOVERLAP_H

#include <iostream>
#include "RecoAlg/CMTool/CMToolBase/CBoolAlgoBase.h"
#include "Utilities/GeometryUtilities.h"


namespace cmtool {
  /**
     \class CMalgoPolyContain
     Merge Polygons if the two overlap even partially
  */
  class CBAlgoPolyOverlap : public CBoolAlgoBase{
    
  public:
    
    /// Default constructor
    CBAlgoPolyOverlap();
    
    /// Default destructor
    virtual ~CBAlgoPolyOverlap(){};
 
    /**
       Core function: given the ClusterParamsAlg input, return whether a cluster should be
       merged or not.
    */
    virtual bool Bool(const ::cluster::cluster_params &cluster1,
                      const ::cluster::cluster_params &cluster2);


    void SetDebug(bool debug) { _debug = debug; }

    //both clusters must have > this # of hits to be considered for merging
    void SetMinNumHits(size_t nhits) { _min_hits = nhits; }

    void SetOverlapFraction(float);

    /// Method to re-configure the instance
    void reconfigure();

  private:
    
    bool _debug;
    size_t _min_hits;
    float _overlap_fraction;  // if used, the smaller cluster has to be
                              // more than this percent contained with 
                              // the larger cluster
  };
}

#endif
/** @} */ // end of doxygen group 

