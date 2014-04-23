////////////////////////////////////////////////////////////////////////
/// \file  CornerFinderAlgAlg.h
/// \brief algorithm to find feature 2D points
///
/// \author  
////////////////////////////////////////////////////////////////////////

#ifndef CORNERFINDERALG_H
#define CORNERFINDERALG_H

#include <vector>
#include <string>
#include "RecoBase/Wire.h"
#include "RecoBase/EndPoint2D.h"
#include "Geometry/Geometry.h"

#include "Eigen/Core"
#include "Eigen/SparseCore"

namespace trkf {
  class BezierTrack;
}

namespace corner {

   class CornerFinderAlg {
   
   public:
    
     CornerFinderAlg(); 
     ~CornerFinderAlg();        
     
     //this one creates the histograms we want to use
     void GrabWires( std::vector<recob::Wire> const&, geo::Geometry const&);                                      

     //and this runs the corner-finding
     void GetFeaturePoints(std::vector<recob::EndPoint2D>&, geo::Geometry const&);

     //setters for our configurations
     void setSparsify(bool option=true){ fSparsify = option; }
     void setSparseReserveSize(unsigned int size=1e6)
     { fSparseReserveSize = size; }
     void setSmoothKernelParams(unsigned int nX, unsigned int nY,
				float sigmaX, float sigmaY)
     { fSmoothNeighborhoodX = nX; fSmoothNeighborhoodY = nY;
       fSmoothWidthX = sigmaX;    fSmoothWidthY = sigmaY; }
     void setImageAlg_Threshold(float i_adc_t, float c_adc_t,
				float i_pass_v=-1, float c_pass_v=-1,
				float i_fail_v=0, float c_fail_v=0){ 
       fImageTransformAlg="Threshold";
       fInductionADCThreshold = i_adc_t;
       fInductionAboveThresholdValue = i_pass_v;
       fInductionBelowThresholdValue = i_fail_v;
       fCollectionADCThreshold = c_adc_t;
       fCollectionAboveThresholdValue = c_pass_v;
       fCollectionBelowThresholdValue = c_fail_v;
     }
     void setImageAlg_Nothing() { fImageTransformAlg="Nothing"; }
     void setStructureTensorWindow( std::string type, unsigned int nhood)
     { fSTWindowType = type; fSTWindowNeighborhood = nhood; }
     void setCornerScoreOptions(std::string method, double min, unsigned int local_nhood)
     { fCornerScoreMethod = method; fCornerScoreMin = min; fLocalMaxNeighborhood = local_nhood; }
     
    private:
     
     //internal configurable options
     bool          fSparsify;
     unsigned int  fSparseReserveSize;

     unsigned int  fSmoothNeighborhoodX;
     unsigned int  fSmoothNeighborhoodY;
     float         fSmoothWidthX;
     float         fSmoothWidthY;

     std::string   fImageTransformAlg;
     float         fInductionADCThreshold;
     float         fCollectionADCThreshold;
     float         fInductionBelowThresholdValue;
     float         fCollectionBelowThresholdValue;
     float         fInductionAboveThresholdValue;
     float         fCollectionAboveThresholdValue;

     std::string   fSTWindowType;
     unsigned int  fSTWindowNeighborhood;

     std::string   fCornerScoreMethod;
     double        fCornerScoreMin;
     unsigned int  fLocalMaxNeighborhood;

     //internal Eigen arrays
     std::vector< Eigen::ArrayXXf >                fWireArrays;
     std::vector< std::vector< geo::WireID > >     fWireIDs;

     void CleanCornerFinderAlg(); //optinal cleanup
     void InitializeGeometry(geo::Geometry const&);
   
     void FillSmoothKernel(Eigen::MatrixXf&);
     void FillWindowKernel(Eigen::MatrixXf&);

     void AttachFeaturePoints( Eigen::ArrayXXf&,
			       std::vector<geo::WireID> const&,
			       geo::View_t, 
			       std::vector<recob::EndPoint2D>&,
			       Eigen::MatrixXf const&,
			       Eigen::MatrixXf const&,
			       Eigen::MatrixXf const&,
			       Eigen::MatrixXf const&);

     void thresholdInput_InPlace(Eigen::ArrayXXf & wireArray,
				 float const threshold_val,
				 float const output_above_threshold=-1,
				 float const output_below_threshold=0);
     void thresholdInput_SparseImage(Eigen::ArrayXXf const&,
				     Eigen::SparseMatrix<float>&,
				     float const threshold_val,
				     float const output_above_threshold=-1);
     
     void construct_Derivative(Eigen::ArrayXXf const&,
			       Eigen::MatrixXf const&,
			       Eigen::MatrixXf const&,
			       Eigen::ArrayXXf &);

     void construct_SparseDerivative(Eigen::SparseMatrix<float> const&,
				     Eigen::MatrixXf const&,
				     Eigen::MatrixXf const&,
				     Eigen::SparseMatrix<float> &);

     void construct_CornerScore(Eigen::ArrayXXf const&, 
				Eigen::ArrayXXf const&,
				Eigen::MatrixXf const&,
				Eigen::ArrayXXd&);     

     void fill_CornerVector(Eigen::ArrayXXf const& wireArray,
			    Eigen::ArrayXXd const& cornerScoreArray,
			    std::vector<recob::EndPoint2D> & corner_vector,
			    std::vector<geo::WireID>, 
			    geo::View_t view);


     
     void construct_SparseCornerScore(Eigen::SparseMatrix<float> const&,
				      Eigen::SparseMatrix<float> const&,
				      Eigen::MatrixXf const&,
				      Eigen::SparseMatrix<double> &);
     
     void fill_SparseCornerVector(Eigen::ArrayXXf const&,
				  Eigen::SparseMatrix<double> const&,
				  std::vector<recob::EndPoint2D> &,
				  std::vector<geo::WireID>,
				  geo::View_t view);

     float Gaussian_2D(float x, float y, 
		       float mean_x, float mean_y,
		       float sigma_x=1, float sigma_y=1,
		       float amp = 1.);

   };//<---End of class CornerFinderAlg
   
   
}

#endif //CORNERFINDERALG_H
