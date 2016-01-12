#ifndef FILLPOLYGON_CXX
#define FILLPOLYGON_CXX

#include "FillPolygon.h"
#include "RecoAlg/ClusterRecoUtil/CRUException.h"
#include "Utilities/GeometryHelper.h"
#include <map>

namespace larlitecluster {

  void FillPolygon::do_params_fill(::cluster::cluster_params & cluster){

    // Geometry Utilities
    auto geomHelper = ::util::GeometryHelper::GetME();

    // get hit list
    auto const& hits = cluster.hit_vector;
    
    if(hits.size()) {
      std::vector<const Hit2D*> polygonEdges;
      geomHelper->SelectPolygonHitList(hits,polygonEdges,_frac);
      //now making Polygon Object
      std::pair<float,float> tmpvertex;
      //make Polygon Object as in mac/PolyOverlap.cc
      std::vector<std::pair<float,float> > vertices;
      for (unsigned int i=0; i < polygonEdges.size(); i++){
        tmpvertex = std::make_pair( polygonEdges.at(i)->w,
				    polygonEdges.at(i)->t );
        vertices.push_back( tmpvertex );
      }
      cluster.PolyObject = Polygon2D( vertices );
    }

    return;
  }
  
} // cluster

#endif
