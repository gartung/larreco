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

#ifndef BlurredCluster_h
#define BlurredCluster_h

// Framework includes
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "larcore/Geometry/Geometry.h"

// ROOT
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

// c++
#include <string>
#include <vector>
#include <map>
#include <sstream>


namespace cluster {
  class BlurredClusterAlg;
}

class cluster::BlurredClusterAlg {
public:

  BlurredClusterAlg(fhicl::ParameterSet const& pset);
  virtual ~BlurredClusterAlg();

  void reconfigure(fhicl::ParameterSet const&p);

  /// Create the PDF to save debug images
  void CreateDebugPDF(int run, int subrun, int event);

  /// Projects a 3D point onto a plane and returns the bin number this point corresponds to on the plane
  TVector2 Convert3DPointToPlaneBins(const TVector3& point, int plane);

  /// Takes a vector of clusters (itself a vector of hits) and turns them into clusters using the initial hit selection
  void ConvertBinsToClusters(std::vector<std::vector<double> > const& image,
			     std::vector<std::vector<int> > const& allClusterBins,
			     std::vector<art::PtrVector<recob::Hit> >& clusters);

  /// Takes hit map and returns a 2D vector representing wire and tick, filled with the charge
  std::vector<std::vector<double> > ConvertRecobHitsToVector(std::vector<art::Ptr<recob::Hit> > const& hits);

  /// Find clusters in the histogram
  std::vector<std::vector<int> > FindClusters(std::vector<std::vector<double> > const& image, const std::vector<TVector2>& vertices, bool reblur = false);

  /// Find the global wire position
  int GlobalWire(geo::WireID const& wireID);

  /// Applies Gaussian blur to image
  std::vector<std::vector<double> > GaussianBlur(std::vector<std::vector<double> > const& image, int bin = -1);

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

private:

  /// Return the coordinates of the specified bin in the hit map (units [cm])
  TVector2 ConvertBinTo2DPosition(const std::vector<std::vector<double> >& image, int bin);

  /// Converts a vector of bins into a hit selection - not all the hits in the bins vector are real hits
  art::PtrVector<recob::Hit> ConvertBinsToRecobHits(std::vector<std::vector<double> > const& image, std::vector<int> const& bins);

  /// Converts a bin into a recob::Hit (not all of these bins correspond to recob::Hits - some are fake hits created by the blurring)
  art::Ptr<recob::Hit> ConvertBinToRecobHit(std::vector<std::vector<double> > const& image, int bin);

  /// Converts a global bin number to an xbin and a ybin
  void ConvertBinToWireTickBins(int bin, const std::vector<std::vector<double> >& blurred, int& xbin, int& ybin);

  /// Converts an xbin and a ybin to a global bin number
  int ConvertWireTickBinsToBin(std::vector<std::vector<double> > const& image, int xbin, int ybin);

  /// Returns the charge stored in the global bin value
  double ConvertBinToCharge(std::vector<std::vector<double> > const& image, int bin);

  /// Count how many dead wires there are in the blurring region for a particular hit
  /// Returns a pair of counters representing how many dead wires there are below and above the hit respectively
  std::pair<int,int> DeadWireCount(int wire_bin, int width);

  /// Dynamically find the blurring radii and Gaussian sigma in each dimension
  /// Use the full hit map to determine the rough direction of the showers (good for single particle events)
  void FindBlurringParameters(int& blurwire, int& blurtick, int& sigmawire, int& sigmatick);

  /// Dynamically find the blurring radii and Gaussian sigma in each dimension
  /// Look over a rangle of angles about a specific point to determine the rough direction of the shower
  /// Intended to be used in multi shower events as part of 'reblurring'
  void FindBlurringParameters(int& blurwire, int& blurtick, int& sigmawire, int& sigmatick, TVector2 point);

  /// Return the coordinates of this hit in global wire/tick space
  TVector2 HitCoordinates(art::Ptr<recob::Hit> const& hit);

  /// Return the coordinates of this hit in units of cm
  TVector2 HitPosition(art::Ptr<recob::Hit> const& hit);

  /// Return the coordinates of this hit in units of cm
  TVector2 HitPosition(TVector2 const& pos, geo::PlaneID planeID);

  /// Makes all the kernels which could be required given the tuned parameters
  void MakeKernels();

  /// Determines the number of clustered neighbours of a hit
  unsigned int NumNeighbours(int nx, std::vector<bool> const& used, int bin);

  /// Projects a 3D point (units [cm]) onto a 2D plane
  /// Returns 2D point (units [wire/tick])
  TVector2 Project3DPointOntoPlane(TVector3 const& point, int plane, int cryostat = 0);

  /// Returns bool indiciating whether or not the specified bin corresponds to a real hit
  bool RealHit(std::vector<std::vector<double> > const& image, int bin);

  bool fDebug;
  std::string fDetector;

  // Parameters used in the Blurred Clustering algorithm
  int          fBlurWire;                 // blur radius for Gauss kernel in the wire direction
  int          fBlurTick;                 // blur radius for Gauss kernel in the tick direction
  double       fSigmaWire;                // sigma for Gaussian kernel in the wire direction
  double       fSigmaTick;                // sigma for Gaussian kernel in the tick direction
  int          fMaxTickWidthBlur;         // maximum distance to blur a hit based on its natural width in time
  double       fShowerDirectionWidth;     // the maximum rectangle width to consider when trying to determine the direction of a shower in reblurring
  int          fClusterWireDistance;      // how far to cluster from seed in wire direction
  int          fClusterTickDistance;      // how far to cluster from seed in tick direction
  unsigned int fMinMergeClusterSize;      // minimum size of a cluster to consider merging it to another
  unsigned int fNeighboursThreshold;      // min. number of neighbors to add to cluster
  int          fMinNeighbours;            // minumum number of neighbors to keep in the cluster
  unsigned int fMinSize;                  // minimum size for cluster
  double       fMinSeed;                  // minimum seed after blurring needed before clustering proceeds
  double       fChargeThreshold;          // charge threshold for clustering

  // Blurring stuff
  int fKernelWidth, fKernelHeight;
  std::vector<std::vector<std::vector<double> > > fAllKernels;

  // Hit containers
  std::vector<std::vector<art::Ptr<recob::Hit> > > fHitMap;

  // Other useful information
  std::vector<bool> fDeadWires;

  int fLowerTick, fUpperTick;
  int fLowerWire, fUpperWire;

  // For the debug pdf
  TCanvas* fDebugCanvas;
  std::string fDebugPDFName;

  // art service handles
  art::ServiceHandle<geo::Geometry> fGeom;
  detinfo::DetectorProperties const* fDetProp;
  lariov::ChannelStatusProvider const& fChanStatus = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();

};

#endif
