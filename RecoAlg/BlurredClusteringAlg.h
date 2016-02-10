////////////////////////////////////////////////////////////////////
// Implementation of the Blurred Clustering algorithm
//
// Converts a hit map into a 2D image of the hits before convoling
// with a Gaussian function to introduce a weighted blurring.
// Clustering proceeds on this blurred image to create more
// complete clusters.
//
// M Wallbank (m.wallbank@sheffield.ac.uk), May 2015
////////////////////////////////////////////////////////////////////

#ifndef BlurredClustering_h
#define BlurredClustering_h

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "Geometry/Geometry.h"
#include "SimulationBase/MCParticle.h"
#include "MCCheater/BackTracker.h"

// ROOT & C++
#include <TTree.h>
#include <TH2F.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TString.h>
#include <TMarker.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TF1.h>
#include <TLine.h>
#include <TPrincipal.h>
#include <TMath.h>
#include <TVector.h>
#include <TVectorD.h>
#include <TVector2.h>

#include <string>
#include <vector>
#include <map>
#include <sstream>


namespace cluster {
  class BlurredClusteringAlg;
}

class cluster::BlurredClusteringAlg {
public:

  BlurredClusteringAlg(fhicl::ParameterSet const& pset);
  virtual ~BlurredClusteringAlg();

  void reconfigure(fhicl::ParameterSet const&p);

  /// Create the PDF to save debug images
  void CreateDebugPDF(int run, int subrun, int event);

  /// Takes a vector of clusters (itself a vector of hits) and turns them into clusters using the initial hit selection
  void ConvertBinsToClusters(std::vector<std::vector<double> > const& image,
			     std::vector<std::vector<int> > const& allClusterBins,
			     std::vector<art::PtrVector<recob::Hit> >& clusters);

  /// Takes hit map and returns a 2D vector representing wire and tick, filled with the charge
  std::vector<std::vector<double> > ConvertRecobHitsToVector(std::vector<art::Ptr<recob::Hit> > const& hits);

  /// Find clusters in the histogram
  int FindClusters(std::vector<std::vector<double> > const& image, std::vector<std::vector<int> >& allcluster);

  /// Find the global wire position
  int GlobalWire(geo::WireID const& wireID);

  /// Applies Gaussian blur to image
  std::vector<std::vector<double> > GaussianBlur(std::vector<std::vector<double> > const& image);

  /// Minimum size of cluster to save
  unsigned int GetMinSize() { return fMinSize; }

  /// Converts a 2D vector in a histogram for the debug pdf
  TH2F* MakeHistogram(std::vector<std::vector<double> > const& image, TString name);

  /// Save the images for debugging
  /// This version takes the final clusters and overlays on the hit map
  void SaveImage(TH2F* image, std::vector<art::PtrVector<recob::Hit> > const& allClusters, int pad, int tpc, int plane);

  /// Save the images for debugging
  void SaveImage(TH2F* image, int pad, int tpc, int plane);

  /// Save the images for debugging
  /// This version takes a vector of bins and overlays the relevant bins on the hit map
  void SaveImage(TH2F* image, std::vector<std::vector<int> > const& allClusterBins, int pad, int tpc, int plane);

  std::map<int,std::map<int,art::Ptr<recob::Hit> > > fHitMap;

private:

  /// Converts a vector of bins into a hit selection - not all the hits in the bins vector are real hits
  art::PtrVector<recob::Hit> ConvertBinsToRecobHits(std::vector<std::vector<double> > const& image, std::vector<int> const& bins);

  /// Converts a bin into a recob::Hit (not all of these bins correspond to recob::Hits - some are fake hits created by the blurring)
  art::Ptr<recob::Hit> ConvertBinToRecobHit(std::vector<std::vector<double> > const& image, int bin);

  /// Convolves the Gaussian kernel with the image to blur
  std::vector<std::vector<double> > Convolve(std::vector<std::vector<double> > const& image, std::vector<double> const& kernel, int width, int height);

  /// Converts an xbin and a ybin to a global bin number                                                                                                                       
  int ConvertWireTickToBin(std::vector<std::vector<double> > const& image, int xbin, int ybin);

  /// Returns the charge stored in the global bin value                                                                                                                        
  double ConvertBinToCharge(std::vector<std::vector<double> > const& image, int bin);

  /// Dynamically find the blurring radii and Gaussian sigma in each dimension
  void FindBlurringParameters(int& blurwire, int& blurtick, int& sigmawire, int& sigmatick);

  /// Returns the hit time of a hit in a particular bin
  double GetTimeOfBin(std::vector<std::vector<double> > const& image, int bin);

  /// Determines the number of clustered neighbours of a hit
  unsigned int NumNeighbours(int nx, std::vector<bool> const& used, int bin);

  /// Determine if a hit is within a time threshold of any other hits in a cluster
  bool PassesTimeCut(std::vector<double> const& times, double time);

  // Parameters used in the Blurred Clustering algorithm
  int          fBlurWire;                 // blur radius for Gauss kernel in the wire direction
  int          fBlurTick;                 // blur radius for Gauss kernel in the tick direction
  double       fBlurSigma;                // sigma for Gaussian kernel
  int          fClusterWireDistance;      // how far to cluster from seed in wire direction
  int          fClusterTickDistance;      // how far to cluster from seed in tick direction
  unsigned int fMinMergeClusterSize;      // minimum size of a cluster to consider merging it to another
  double       fMergingThreshold;        // the PCA eigenvalue needed to consider two clusters a merge
  unsigned int fNeighboursThreshold;      // min. number of neighbors to add to cluster
  int          fMinNeighbours;            // minumum number of neighbors to keep in the cluster
  unsigned int fMinSize;                  // minimum size for cluster
  double       fMinSeed;                  // minimum seed after blurring needed before clustering proceeds
  double       fTimeThreshold;            // time threshold for clustering
  double       fChargeThreshold;          // charge threshold for clustering

  // Wire and tick information for histograms
  int fLowerHistTick, fUpperHistTick;
  int fLowerHistWire, fUpperHistWire;

  // Blurring stuff
  int fLastBlurWire;
  int fLastBlurTick;
  double fLastSigma;
  std::vector<double> fLastKernel;

  // For the debug pdf
  TCanvas *fDebugCanvas;
  std::string fDebugPDFName;

  // art service handles
  art::ServiceHandle<geo::Geometry> fGeom;
  art::ServiceHandle<util::DetectorProperties> fDetProp;

};

#endif
