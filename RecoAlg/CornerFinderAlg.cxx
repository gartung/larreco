////////////////////////////////////////////////////////////////////////
//
// CornerFinderAlg class
//
// wketchum@fnal.gov
//
// CornerFinder is meant to use image-processing techniques (mainly Harris-Stephens
// corner-finding) to find "corners" using the information from calibrated wires.
//  
//  ImageTransform options:  
//     Threshold --- only uses wire info above a certain value. Can use a 
//                   constant value for ticks above threshold, or use wire data.
//     Nothing   --- takes a copy of the input wire data.
//       
//  CornerScore_algorithm options:
//     Noble        --- determinant / (trace + 1e-5)
//     Harris       --- determinant - (trace)^2 *0.05
//     STWindowType --- Can be "Flat" or "Gaussian"; this is applied to elements 
//                      forming structure tensor
////////////////////////////////////////////////////////////////////////


#include <cmath>
#include "RecoAlg/CornerFinderAlg.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"

#include "RecoObjects/BezierTrack.h"

#include "Eigen/LU"

//-----------------------------------------------------------------------------
corner::CornerFinderAlg::CornerFinderAlg(){
  fSparsify = true;
  fSparseReserveSize = 1e6;
  fSmoothNeighborhoodX = 4;
  fSmoothNeighborhoodY = 4;
  fSmoothWidthX = 1.;
  fSmoothWidthY = 1.;

  fImageTransformAlg = "";
  fInductionADCThreshold = 6;
  fCollectionADCThreshold = 6;

  fSTWindowType = "Gaussian";
  fSTWindowNeighborhood = 2;

  fCornerScoreMethod = "Harris";
  fCornerScoreMin = 1e10;
  fLocalMaxNeighborhood = 2;

}

corner::CornerFinderAlg::~CornerFinderAlg(){}

//-----------------------------------------------------------------------------
void corner::CornerFinderAlg::CleanCornerFinderAlg(){}

//-----------------------------------------------------------------------------
void corner::CornerFinderAlg::GrabWires(std::vector<recob::Wire> const& wireVec,
					geo::Geometry const& my_geometry){

  InitializeGeometry(my_geometry);

  const unsigned int nTimeTicks = wireVec.at(0).NSignal();
  const unsigned int nPlanes = my_geometry.Nplanes();

  for( unsigned int i_plane=0; i_plane < nPlanes; i_plane++)
    fWireArrays.at(i_plane).resize(my_geometry.Nwires(i_plane),nTimeTicks);

  for(auto iwire : wireVec){
    
    geo::WireID this_wireID;
    try { this_wireID = (my_geometry.ChannelToWire(iwire.Channel())).at(0);}
    catch(cet::exception& excep) { LOG_ERROR("CornerFinderAlg") << "Bail out! No Possible Wires!\n"; }
    
    unsigned int i_plane = this_wireID.Plane;
    unsigned int i_wire  = this_wireID.Wire;
    
    fWireIDs.at(i_plane).at(i_wire) = this_wireID;
    
    //Note, we make a copy of the wire data here, wire by wire. 
    //We will overwrite fWireArrays throughout the algorithm, so it's not so wasteful.
    for(unsigned int tick=0; tick<iwire.Signal().size(); tick++)
      (fWireArrays.at(i_plane))(i_wire,tick) = iwire.Signal().at(tick); 
  }

}

//-----------------------------------------------------------------------------------
// This gives us a vecotr of EndPoint2D objects that correspond to possible corners
void corner::CornerFinderAlg::GetFeaturePoints(std::vector<recob::EndPoint2D> & corner_vector, 
					       geo::Geometry const& my_geometry){


  //declare some kernels used for derivatives.
  //these could be changed to be options declared from input parameters...
  Eigen::MatrixXf SMOOTH_KERNEL(fSmoothNeighborhoodX*2+1,fSmoothNeighborhoodY*2+1);
  FillSmoothKernel(SMOOTH_KERNEL);
  
  const Eigen::Matrix3f SOBEL_KERNEL_X = 
    (Eigen::Matrix3f() <<
     -0.25, 0.00, 0.25,
     -0.50, 0.00, 0.50,
     -0.25, 0.00, 0.50).finished();
  
  const Eigen::Matrix3f SOBEL_KERNEL_Y = 
    (Eigen::Matrix3f() <<
     -0.25, -0.50, -0.25,
      0.00,  0.00,  0.00,
      0.25,  0.50,  0.25).finished();

  Eigen::MatrixXf WINDOW_KERNEL(fSTWindowNeighborhood*2+1,fSTWindowNeighborhood*2+1);
  FillWindowKernel(WINDOW_KERNEL);


  for(auto pid : my_geometry.PlaneIDs())
    AttachFeaturePoints(fWireArrays.at(pid.Plane),
			fWireIDs.at(pid.Plane),
			my_geometry.View(pid),
			corner_vector,
			SOBEL_KERNEL_X,
			SOBEL_KERNEL_Y,
			SMOOTH_KERNEL,
			WINDOW_KERNEL);

}


//-----------------------------------------------------------------------------
void corner::CornerFinderAlg::InitializeGeometry(geo::Geometry const& my_geometry){

  CleanCornerFinderAlg();
  fWireArrays.resize(my_geometry.Nplanes());

  fWireIDs.resize(my_geometry.Nplanes());
  for(unsigned int i_plane=0; i_plane < my_geometry.Nplanes(); ++i_plane)
    fWireIDs.at(i_plane).resize(my_geometry.Nwires(i_plane));
}

//-----------------------------------------------------------------------------
void corner::CornerFinderAlg::FillSmoothKernel(Eigen::MatrixXf & SMOOTH_KERNEL){

  for(unsigned int i=0; i<SMOOTH_KERNEL.rows(); i++){
    for(unsigned int j=0; j<SMOOTH_KERNEL.cols(); j++){
      SMOOTH_KERNEL(i,j) = Gaussian_2D((float)i,             //xval
				       (float)j,             //yval
				       1,                    //amp
				       fSmoothNeighborhoodX, //xmean  
				       fSmoothNeighborhoodY, //ymean
				       fSmoothWidthX,        //xwidth
				       fSmoothWidthY);       //ywidth
    }
  }
  float integral = SMOOTH_KERNEL.sum();
  SMOOTH_KERNEL /= integral;
}

//-----------------------------------------------------------------------------
void corner::CornerFinderAlg::FillWindowKernel(Eigen::MatrixXf & WINDOW_KERNEL){

  if(fSTWindowType.compare("Gaussian")==0){
    FillSmoothKernel(WINDOW_KERNEL);
    return;
  }
  else if(fSTWindowType.compare("Flat")==0){
    int nRows = WINDOW_KERNEL.rows();
    int nCols = WINDOW_KERNEL.cols();
    float value = 1./(float)(nRows+nCols);
    WINDOW_KERNEL = Eigen::MatrixXf::Constant(nRows,nCols,value);
    return;
  }
  else{
      LOG_ERROR("CornerFinderAlg") << "No such window type option "
				   << fSTWindowType;
      return;
  }

}

//-----------------------------------------------------------------------------
float corner::CornerFinderAlg::Gaussian_2D(float x, float y, 
					   float mean_x, float mean_y,
					   float sigma_x, float sigma_y,
					   float amp){
  return amp* std::exp( -0.5*((x-mean_x)*(x-mean_x)/sigma_x/sigma_x + 
			      (y-mean_y)*(y-mean_y)/sigma_y/sigma_y));
}

//-----------------------------------------------------------------------------
// This puts on all the feature points in a given view, using a given data histogram
void corner::CornerFinderAlg::AttachFeaturePoints( Eigen::ArrayXXf & wireArray, 
						   std::vector<geo::WireID> const& wireIDs, 
						   geo::View_t view, 
						   std::vector<recob::EndPoint2D> & corner_vector,
						   Eigen::MatrixXf const& kernel_derivative_x,
						   Eigen::MatrixXf const& kernel_derivative_y,
						   Eigen::MatrixXf const& kernel_derivative_smooth,
						   Eigen::MatrixXf const& kernel_window){
  
  const unsigned int N_ROWS = wireArray.rows();
  const unsigned int N_COLS = wireArray.cols();

  
  if(fSparsify){

    //make image array
    Eigen::SparseMatrix<float> imageSMatrix(N_ROWS,N_COLS);
    imageSMatrix.reserve(fSparseReserveSize);

    if(fImageTransformAlg.compare("Threshold")==0){
      if(view==geo::View_t::kU || view==geo::View_t::kV)
	thresholdInput_SparseImage(wireArray,
				   imageSMatrix,
				   fInductionADCThreshold,
				   fInductionAboveThresholdValue);
      else if(view==geo::View_t::kZ)
	thresholdInput_SparseImage(wireArray,
				   imageSMatrix,
				   fCollectionADCThreshold,
				   fCollectionAboveThresholdValue);
    }
    else{
      LOG_ERROR("CornerFinderAlg") << "No such image transform option "
				   << fImageTransformAlg
				   << " for Sparse Matrices.";
      return;
    }

    imageSMatrix.makeCompressed();

    //declare some internal temporary arrays
    const unsigned int new_reserve_size = 2*imageSMatrix.nonZeros();

    Eigen::SparseMatrix<float>  derivativeXSMatrix(N_ROWS,N_COLS);
    Eigen::SparseMatrix<float>  derivativeYSMatrix(N_ROWS,N_COLS);
    Eigen::SparseMatrix<double> cornerScoreSMatrix(N_ROWS,N_COLS);
    
    derivativeXSMatrix.reserve(new_reserve_size);
    derivativeYSMatrix.reserve(new_reserve_size);
    cornerScoreSMatrix.reserve(new_reserve_size);
    
    //do the derivatives
    construct_SparseDerivative(imageSMatrix,kernel_derivative_x,kernel_derivative_smooth,derivativeXSMatrix);
    construct_SparseDerivative(imageSMatrix,kernel_derivative_y,kernel_derivative_smooth,derivativeYSMatrix);
    
    derivativeXSMatrix.makeCompressed();
    derivativeYSMatrix.makeCompressed();

    //get the corner score
    construct_SparseCornerScore(derivativeXSMatrix,derivativeYSMatrix,kernel_window,cornerScoreSMatrix);
    
    cornerScoreSMatrix.makeCompressed();

    //fill the corner vector
    fill_SparseCornerVector(wireArray,cornerScoreSMatrix,corner_vector,wireIDs,view);
    
  }//end if sparsify
  
  else{

    //make image array
    if(fImageTransformAlg.compare("Threshold")==0){
      if(view==geo::View_t::kU || view==geo::View_t::kV)
	thresholdInput_InPlace(wireArray,
			       fInductionADCThreshold,
			       fInductionAboveThresholdValue,
			       fInductionBelowThresholdValue);
      else if(view==geo::View_t::kZ)
	thresholdInput_InPlace(wireArray,
			       fCollectionADCThreshold,
			       fCollectionAboveThresholdValue,
			       fCollectionBelowThresholdValue);
    }
    else if(fImageTransformAlg.compare("Nothing")==0){
    }
    else{
      LOG_ERROR("CornerFinderAlg") << "No such image transform option "
				   << fImageTransformAlg
				   << " for full arrays.";
      return;
    }

    //set up temporary arrays
    Eigen::ArrayXXf derivativeXArray = Eigen::ArrayXXf::Zero(N_ROWS,N_COLS);
    Eigen::ArrayXXf derivativeYArray = Eigen::ArrayXXf::Zero(N_ROWS,N_COLS);
    Eigen::ArrayXXd cornerScoreArray = Eigen::ArrayXXd::Zero(N_ROWS,N_COLS);
    
    //make derivatives
    construct_Derivative(wireArray,kernel_derivative_x,kernel_derivative_smooth,derivativeXArray);
    construct_Derivative(wireArray,kernel_derivative_y,kernel_derivative_smooth,derivativeYArray);
        
    //get corner score
    construct_CornerScore(derivativeXArray,derivativeYArray,kernel_window,cornerScoreArray);
    
    //fill corner vector
    fill_CornerVector(wireArray,cornerScoreArray,corner_vector,wireIDs,view);

  }//end else (not sparse)

}//end AttachFeaturePoints




void corner::CornerFinderAlg::thresholdInput_InPlace(Eigen::ArrayXXf & wireArray,
						     float const threshold_val,
						     float const output_above_threshold,
						     float const output_below_threshold){

  const unsigned int nRows = wireArray.rows();
  const unsigned int nCols = wireArray.cols();

  for(unsigned int irow=0; irow<nRows; irow++){
    for(unsigned int icol=0; icol<nCols; icol++){

      if(wireArray(irow,icol) < threshold_val)
	wireArray(irow,icol) = output_below_threshold;

      else if(output_above_threshold > output_below_threshold)
	wireArray(irow,icol) = output_above_threshold;

    }//end cols
  }//end rows

}

void corner::CornerFinderAlg::thresholdInput_SparseImage(Eigen::ArrayXXf const& wireArray,
							 Eigen::SparseMatrix<float> & imageSMatrix,
							 float const threshold_val,
							 float const output_above_threshold){

  const unsigned int nRows = wireArray.rows();
  const unsigned int nCols = wireArray.cols();

  std::vector< Eigen::Triplet<float> > triplet_vector;

  for(unsigned int irow=0; irow<nRows; irow++){
    for(unsigned int icol=0; icol<nCols; icol++){

      if(wireArray(irow,icol) > threshold_val){
	if(output_above_threshold > 0)
	  triplet_vector.emplace_back(irow,icol,output_above_threshold);
	else
	  triplet_vector.emplace_back(irow,icol,wireArray(irow,icol));
      }//end if above threshold

    }//end cols
  }//end rows

  imageSMatrix.setFromTriplets(triplet_vector.begin(),triplet_vector.end());
}


void corner::CornerFinderAlg::construct_SparseDerivative(Eigen::SparseMatrix<float> const& imageSMatrix,
							 Eigen::MatrixXf const& KERNEL_DERIVATIVE,
							 Eigen::MatrixXf const& KERNEL_SMOOTH,
							 Eigen::SparseMatrix<float> & derivativeSMatrix){
  
  const int nRows = imageSMatrix.rows();
  const int nCols = imageSMatrix.cols();

  const int d_neighborhood = (KERNEL_DERIVATIVE.rows() - 1)/2;
  const int d_rows = KERNEL_DERIVATIVE.rows();
  const int d_cols = KERNEL_DERIVATIVE.cols();
  const int s_neighborhood = (KERNEL_SMOOTH.rows() - 1)/2;
  const int s_rows = KERNEL_SMOOTH.rows();
  const int s_cols = KERNEL_SMOOTH.cols();


  std::vector< Eigen::Triplet<float> > triplet_vector;

  for(int i_outer=0; i_outer < imageSMatrix.outerSize(); i_outer++){
    for(Eigen::SparseMatrix<float>::InnerIterator it(imageSMatrix,i_outer); it; ++it){

      if( it.row()<d_neighborhood || it.col()<d_neighborhood || 
	  it.row()>(nRows-1-d_neighborhood) || it.col()>(nCols-1-d_neighborhood)) continue;

      float sum = (imageSMatrix.block(it.row()-d_neighborhood,
				      it.col()-d_neighborhood,
				      d_rows,
				      d_cols).cwiseProduct(KERNEL_DERIVATIVE)).sum();

      if(std::abs(sum)>1e-3)
	triplet_vector.emplace_back(it.row(),it.col(),sum);

    }
  }

  Eigen::SparseMatrix<float> tmp_SMatrix(nRows,nCols);
  tmp_SMatrix.setFromTriplets(triplet_vector.begin(),triplet_vector.end());


  triplet_vector.clear();

  for(int i_outer=0; i_outer < tmp_SMatrix.outerSize(); i_outer++){
    for(Eigen::SparseMatrix<float>::InnerIterator it(tmp_SMatrix,i_outer); it; ++it){
      
      if( it.row()<s_neighborhood || it.col()<s_neighborhood || 
	  it.row()>(nRows-1-s_neighborhood) || it.col()>(nCols-1-s_neighborhood)) continue;
      
      float sum = (tmp_SMatrix.block(it.row()-s_neighborhood,
				     it.col()-s_neighborhood,
				     s_rows,
				     s_cols).cwiseProduct(KERNEL_SMOOTH)).sum();
      
      if(std::abs(sum)>1e-3)
	triplet_vector.emplace_back(it.row(),it.col(),sum);

    }
  }
  
  std::cout << "Triplet size is " << triplet_vector.size() << std::endl;
  derivativeSMatrix.setFromTriplets(triplet_vector.begin(),triplet_vector.end());

}


void corner::CornerFinderAlg::construct_Derivative(Eigen::ArrayXXf const& imageArray,
						   Eigen::MatrixXf const& KERNEL_DERIVATIVE,
						   Eigen::MatrixXf const& KERNEL_SMOOTH,
						   Eigen::ArrayXXf & derivativeXArray){

  const unsigned int nRows = imageArray.rows();
  const unsigned int nCols = imageArray.cols();

  const int d_neighborhood = (KERNEL_DERIVATIVE.rows() - 1)/2;
  const int d_rows = KERNEL_DERIVATIVE.rows();
  const int d_cols = KERNEL_DERIVATIVE.cols();
  const int s_neighborhood = (KERNEL_SMOOTH.rows() - 1)/2;
  const int s_rows = KERNEL_SMOOTH.rows();
  const int s_cols = KERNEL_SMOOTH.cols();

  Eigen::ArrayXXf tmp_Array = derivativeXArray;

  for(unsigned int irow=d_neighborhood; irow < nRows-d_neighborhood; irow++){
    for(unsigned int icol=d_neighborhood; icol < nCols-d_neighborhood; icol++){
      tmp_Array(irow,icol) = (imageArray.block(irow-d_neighborhood,
					       icol-d_neighborhood,
					       d_rows,
					       d_cols)*(KERNEL_DERIVATIVE.array()) ).sum();
    }
  }

  for(unsigned int irow=s_neighborhood; irow < tmp_Array.rows()-s_neighborhood; irow++){
    for(unsigned int icol=s_neighborhood; icol < tmp_Array.cols()-s_neighborhood; icol++){
      derivativeXArray(irow,icol) = (tmp_Array.block(irow-s_neighborhood,
						     icol-s_neighborhood,
						     s_rows,
						     s_cols)*(KERNEL_SMOOTH.array()) ).sum();
    }
  }
  
}


void corner::CornerFinderAlg::construct_SparseCornerScore(Eigen::SparseMatrix<float> const& derivativeXSMatrix, 
							  Eigen::SparseMatrix<float> const& derivativeYSMatrix,
							  Eigen::MatrixXf const& windowMatrix,
							  Eigen::SparseMatrix<double> & cornerScoreSMatrix){

  Eigen::Matrix2f structure_tensor = Eigen::Matrix2f::Zero();

  const int nRows = derivativeXSMatrix.rows();
  const int nCols = derivativeXSMatrix.cols();

  const int w_rows = windowMatrix.rows();
  const int w_cols = windowMatrix.cols();
  const int w_neighborhood = (w_rows-1)/2;

  std::vector< Eigen::Triplet<double> > triplet_vector;

  for(int i_outer=0; i_outer < derivativeXSMatrix.outerSize(); i_outer++){
    for(Eigen::SparseMatrix<float>::InnerIterator it(derivativeXSMatrix,i_outer); it; ++it){
      
      if( it.row()<w_neighborhood || it.col()<w_neighborhood || 
	  it.row()>(nRows-1-w_neighborhood) || it.col()>(nCols-1-w_neighborhood)) continue;

      Eigen::MatrixXf x_block = derivativeXSMatrix.block(it.row()-w_neighborhood,
							 it.col()-w_neighborhood,
							 w_rows,
							 w_cols);
      Eigen::MatrixXf y_block = derivativeYSMatrix.block(it.row()-w_neighborhood,
							 it.col()-w_neighborhood,
							 w_rows,
							 w_cols);

      //if(x_block.isZero(1e-3) && y_block.isZero(1e-3)) continue;


      structure_tensor(0,0) = (x_block.array()*x_block.array()*windowMatrix.array()).sum();
      structure_tensor(1,1) = (y_block.array()*y_block.array()*windowMatrix.array()).sum();
      structure_tensor(0,1) = (x_block.array()*y_block.array()*windowMatrix.array()).sum();
      structure_tensor(1,0) = structure_tensor(0,1);
      
      double score;
      if(fCornerScoreMethod.compare("Harris")==0)
	score = structure_tensor.determinant() - 0.05 * structure_tensor.trace()*structure_tensor.trace();
      else if(fCornerScoreMethod.compare("Noble")==0)
	score = 2*(structure_tensor.determinant() / (structure_tensor.trace()+1e-5));
      else{
	LOG_ERROR("CornerFinderAlg") << "No such corner score option "
				     << fCornerScoreMethod;
	return;
      }

      if(score>fCornerScoreMin)
	triplet_vector.emplace_back(it.row(),it.col(),score);

    }
  }

  std::cout << "Triplet size is " << triplet_vector.size() << std::endl;

  for(int i_outer=0; i_outer < derivativeYSMatrix.outerSize(); i_outer++){
    for(Eigen::SparseMatrix<float>::InnerIterator it(derivativeYSMatrix,i_outer); it; ++it){
      
      if( it.row()<w_neighborhood || it.col()<w_neighborhood || 
	  it.row()>(nRows-1-w_neighborhood) || it.col()>(nCols-1-w_neighborhood)) continue;

      Eigen::MatrixXf x_block = derivativeXSMatrix.block(it.row()-w_neighborhood,
							 it.col()-w_neighborhood,
							 w_rows,
							 w_cols);
      Eigen::MatrixXf y_block = derivativeYSMatrix.block(it.row()-w_neighborhood,
							 it.col()-w_neighborhood,
							 w_rows,
							 w_cols);

      //if(x_block.isZero(1e-3) && y_block.isZero(1e-3)) continue;


      structure_tensor(0,0) = (x_block.array()*x_block.array()*windowMatrix.array()).sum();
      structure_tensor(1,1) = (y_block.array()*y_block.array()*windowMatrix.array()).sum();
      structure_tensor(0,1) = (x_block.array()*y_block.array()*windowMatrix.array()).sum();
      structure_tensor(1,0) = structure_tensor(0,1);
      
      double score;
      if(fCornerScoreMethod.compare("Harris")==0)
	score = structure_tensor.determinant() - 0.05 * structure_tensor.trace()*structure_tensor.trace();
      else if(fCornerScoreMethod.compare("Noble")==0)
	score = 2*(structure_tensor.determinant() / (structure_tensor.trace()+1e-5));
      else{
	LOG_ERROR("CornerFinderAlg") << "No such corner score option "
				     << fCornerScoreMethod;
	return;
      }

      if(score>fCornerScoreMin)
	triplet_vector.emplace_back(it.row(),it.col(),score);

    }
  }

  std::cout << "Triplet size is " << triplet_vector.size() << std::endl;

  cornerScoreSMatrix.setFromTriplets(triplet_vector.begin(),triplet_vector.end());

}

void corner::CornerFinderAlg::construct_CornerScore(Eigen::ArrayXXf const& derivativeXArray, 
						    Eigen::ArrayXXf const& derivativeYArray,
						    Eigen::MatrixXf const& windowMatrix,
						    Eigen::ArrayXXd & cornerScoreArray){


  Eigen::Matrix2f structure_tensor = Eigen::Matrix2f::Zero();

  const int nRows = derivativeXArray.rows();
  const int nCols = derivativeXArray.cols();

  const int w_rows = windowMatrix.rows();
  const int w_cols = windowMatrix.cols();
  const int w_neighborhood = (w_rows-1)/2;

  for(int irow=w_neighborhood; irow < nRows-1-w_neighborhood; irow++){
    for(int icol=w_neighborhood; icol < nCols-1-w_neighborhood; icol++){

      Eigen::ArrayXXf x_block = derivativeXArray.block(irow-w_neighborhood,
						       icol-w_neighborhood,
						       w_rows,
						       w_cols);
      Eigen::ArrayXXf y_block = derivativeYArray.block(irow-w_neighborhood,
						       icol-w_neighborhood,
						       w_rows,
						       w_cols);
      
      if(x_block.isZero(1e-3) && y_block.isZero(1e-3)) continue;
      
      structure_tensor(0,0) = (x_block*x_block*windowMatrix.array()).sum();
      structure_tensor(1,1) = (y_block*y_block*windowMatrix.array()).sum();
      structure_tensor(0,1) = (x_block*y_block*windowMatrix.array()).sum();
      structure_tensor(1,0) = structure_tensor(0,1);

      double score;
      if(fCornerScoreMethod.compare("Harris")==0)
	score = structure_tensor.determinant() - 0.05 * structure_tensor.trace()*structure_tensor.trace();
      else if(fCornerScoreMethod.compare("Noble")==0)
	score = 2*(structure_tensor.determinant() / (structure_tensor.trace()+1e-5));
      else{
	LOG_ERROR("CornerFinderAlg") << "No such corner score option "
				     << fCornerScoreMethod;
	return;
      }
      cornerScoreArray(irow,icol) = score;
    }
  }

}

void corner::CornerFinderAlg::fill_SparseCornerVector(Eigen::ArrayXXf const& wireArray,
						      Eigen::SparseMatrix<double> const& cornerScoreSMatrix,
						      std::vector<recob::EndPoint2D> & corner_vector,
						      std::vector<geo::WireID> wireIDs, 
						      geo::View_t view){

  const int nRows = cornerScoreSMatrix.rows();
  const int nCols = cornerScoreSMatrix.cols();
  
  unsigned int max_row; unsigned int max_col;

  for(int i_outer=0; i_outer < cornerScoreSMatrix.outerSize(); i_outer++){
    for(Eigen::SparseMatrix<double>::InnerIterator it(cornerScoreSMatrix,i_outer); it; ++it){

      if( it.row()<(int)fLocalMaxNeighborhood || it.col()<(int)fLocalMaxNeighborhood || 
	  it.row()>(nRows-1-(int)fLocalMaxNeighborhood) || it.col()>(nCols-1-(int)fLocalMaxNeighborhood)) continue;

      Eigen::MatrixXd my_block(fLocalMaxNeighborhood*2+1,fLocalMaxNeighborhood*2+1);
      my_block = cornerScoreSMatrix.block(it.row()-fLocalMaxNeighborhood,
					  it.col()-fLocalMaxNeighborhood,
					  fLocalMaxNeighborhood*2+1,
					  fLocalMaxNeighborhood*2+1);
      double max = my_block.maxCoeff(&max_row,&max_col);
      
      if( max_row==fLocalMaxNeighborhood && max_col==fLocalMaxNeighborhood && max>fCornerScoreMin )
	corner_vector.emplace_back(it.col(),
				   wireIDs.at(it.row()),
				   max,
				   0,//id
				   view,
				   wireArray(it.row(),it.col())); //totalQ
      
      
    }
  }

}

void corner::CornerFinderAlg::fill_CornerVector(Eigen::ArrayXXf const& wireArray,
						Eigen::ArrayXXd const& cornerScoreArray,
						std::vector<recob::EndPoint2D> & corner_vector,
						std::vector<geo::WireID> wireIDs, 
						geo::View_t view){

  const unsigned int nRows = cornerScoreArray.rows();
  const unsigned int nCols = cornerScoreArray.cols();
  
  unsigned int max_row; unsigned int max_col;

  for(unsigned int irow=fLocalMaxNeighborhood; irow < nRows-fLocalMaxNeighborhood; irow++){
    for(unsigned int icol=fLocalMaxNeighborhood; icol < nCols-fLocalMaxNeighborhood; icol++){

      float max = (cornerScoreArray.block(irow-fLocalMaxNeighborhood,
					  icol-fLocalMaxNeighborhood,
					  fLocalMaxNeighborhood*2+1,
					  fLocalMaxNeighborhood*2+1)).maxCoeff(&max_row,&max_col);
      
      if(max_row==fLocalMaxNeighborhood && 
	 max_col==fLocalMaxNeighborhood && 
	 max>fCornerScoreMin)
	corner_vector.emplace_back(icol,
				   wireIDs.at(irow),
				   max,
				   0,//id
				   view,
				   wireArray(irow,icol)); //totalQ

    }
  }

}

