/**
 * \file FillGeomParams.h
 *
 * \ingroup ClusterRecoUtil
 * 
 * \brief Class def header for a class FillGeomParams
 *
 * @author cadams
 */

/** \addtogroup ClusterRecoUtil

    @{*/
#ifndef FILLGEOMPARAMS_H
#define FILLGEOMPARAMS_H

#include <iostream>
#include "RecoAlg/ClusterRecoUtil/ParamsAlgBase.h"

namespace larlitecluster {

  /**
     \class FillGeomParams
     User defined class FillGeomParams ... these comments are used to generate
     doxygen documentation!
   */
  class FillGeomParams : public ParamsAlgBase {

  public:

    /// Default constructor
    FillGeomParams();

    /// Default destructor
    ~FillGeomParams(){}

    /**
       Calculates the following variables:

      */
    void do_params_fill(::cluster::cluster_params &);

  protected:

    float _percentage;
    size_t _Nbins;
  };

} // cluster


#endif
/** @} */ // end of doxygen group 

