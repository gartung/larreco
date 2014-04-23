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

namespace corner { //<---Not sure if this is the right namespace

   class CornerFinderAlg {
   
   public:
    
     explicit CornerFinderAlg(fhicl::ParameterSet const& pset); 
     virtual ~CornerFinderAlg();        
     
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
       fSmoothWidthX = sigmaX;    fSmoothSigmaY = sigmaY; }

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


     //internal Eigen arrays
     std::vector< Eigen::ArrayXXf >                fWireArrays;
     std::vector< std::vector< geo::WireID > >     fWireIDs;

     void CleanCornerFinderAlg(); //optinal cleanup
     void InitializeGeometry(geo::Geometry const&);
   
     void FillSmoothKernel(Eigen::MatrixXf&);

     void AttachFeaturePoints( Eigen::ArrayXXf&,
			       std::vector<geo::WireID> const&,
			       geo::View_t, 
			       std::vector<recob::EndPoint2D>&,
			       Eigen::MaxtrixXf const&,
			       Eigen::MaxtrixXf const&,
			       Eigen::MaxtrixXf const&);

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

     void construct_CornerScore(Eigen::ArrayXXf const& derivativeXArray, 
				Eigen::ArrayXXf const& derivativeYArray,
				Eigen::ArrayXXd & cornerScoreArray);     

     void fill_CornerVector(Eigen::ArrayXXf const& wireArray,
			    Eigen::ArrayXXd const& cornerScoreArray,
			    std::vector<recob::EndPoint2D> & corner_vector,
			    std::vector<geo::WireID> wireIDs, 
			    geo::View_t view);


     
     void construct_SparseCornerScore(Eigen::SparseMatrix<float> const&,
				      Eigen::SparseMatrix<float> const&,
				      Eigen::SparseMatrix<double> &);
     
     void fill_SparseCornerVector(Eigen::ArrayXXf const&,
				  Eigen::SparseMatrix<double> const&,
				  std::vector<recob::EndPoint2D> &,
				  std::vector<geo::WireID> const&,
				  geo::View_t view);

     float Gaussian_2D(float x, float y, 
		       float mean_x, float mean_y,
		       float sigma_x=1, float sigma_y=1,
		       float amp = 1.);

     unsigned int event_number;
     unsigned int run_number;
     
     void create_image_histo(TH2F const& h_wire_data, TH2F & h_conversion);
     void create_derivative_histograms(TH2F const& h_conversion, TH2F & h_derivative_x, TH2F & h_derivative_y);
     void create_cornerScore_histogram(TH2F const& h_derivative_x, TH2F const& h_derivative_y, TH2D & h_cornerScore);
     size_t perform_maximum_suppression(TH2D const& h_cornerScore, 
					std::vector<recob::EndPoint2D> & corner_vector,
					std::vector<geo::WireID> wireIDs, 
					geo::View_t view,
					TH2D & h_maxSuppress,
					int startx=0,
					int starty=0);
     
     size_t calculate_line_integral_score( TH2F const& h_wire_data, 
					   std::vector<recob::EndPoint2D> const & corner_vector, 
					   std::vector<recob::EndPoint2D> & corner_lineIntegralScore_vector,
					   TH2F & h_lineIntegralScore);
     
     void attach_feature_points(TH2F const& h_wire_data, 
				std::vector<geo::WireID> wireIDs, 
				geo::View_t view,
				std::vector<recob::EndPoint2D>&,
				int startx=0,int starty=0);
     void attach_feature_points_LineIntegralScore(TH2F const& h_wire_data, 
						  std::vector<geo::WireID> wireIDs, 
						  geo::View_t view,
						  std::vector<recob::EndPoint2D>&);
     
     
     void create_smaller_histos(geo::Geometry const&);
     void remove_duplicates(std::vector<recob::EndPoint2D>&);
     

   };//<---End of class CornerFinderAlg
   
   
}

#endif //CORNERFINDERALG_H
