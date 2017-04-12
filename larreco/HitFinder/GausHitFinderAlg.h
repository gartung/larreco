#ifndef GAUS_HIT_FINDER_ALG_H
#define GAUS_HIT_FINDER_ALG_H

/*
 * Gaus Hit Finder Alg
 *
 * Author: Jonathan Asaadi
 *         Corey Adams
 *
 * The code in this module takes much of the functionality from the original
 * hit finder module and packages it into an algorithm class.  This was done
 * to allow access to those functions from second pass hit finders.
 *
 * JA wrote this software.
 * CA repackaged it with minor organizational changes, no content changes
 *
 *
 *
 */

// Framework includes
#include "canvas/Persistency/Common/FindOneP.h"
#include "art/Framework/Principal/Event.h"


// LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "HitFilterAlg.h"


#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Optional/TFileService.h"


#include "TGraphErrors.h"
#include "TH1D.h"
#include "TDecompSVD.h"
#include "TMath.h"
#include "TF1.h"
#include "TTree.h"
#include "TStopwatch.h"

namespace hit {

class GausHitFinderAlg
{
public:
  GausHitFinderAlg(fhicl::ParameterSet const&);
  ~GausHitFinderAlg(){}

  void RunOnWire(art::Ptr<recob::Wire> const&,
                 std::vector<recob::Hit>&,
                 geo::Geometry const&);

  void RunOnROI(const lar::sparse_vector<float>::datarange_t &,
                std::vector<recob::Hit> &,
                art::Ptr<recob::Wire>,
                geo::Geometry const&);

  void beginJob();
  void endJob(){}

private:

  using TimeValsVec      = std::vector<std::tuple<int, int, int>>;
  using PeakTimeWidVec   = std::vector<std::pair<int, int>>;
  using MergedTimeWidVec = std::vector<std::tuple<int, int, PeakTimeWidVec>>;
  using ParameterVec     = std::vector<std::pair<double, double>>; //< parameter/error vec
  using ROIRegion        = lar::sparse_vector<float>::datarange_t;

  double              threshold           = 0.;  // minimum signal size for id'ing a hit
  double              fitWidth            = 0.;  // hit fit width initial value
  double              minWidth            = 0 ;  // hit minimum width

  std::vector<double> fMinSig;                   ///<signal height threshold
  std::vector<double> fInitWidth;                ///<Initial width for fit
  std::vector<double> fMinWidth;                 ///<Minimum hit width
  std::vector<int>    fLongMaxHits;              ///<Maximum number hits on a really long pulse train
  std::vector<int>    fLongPulseWidth;           ///<Sets width of hits used to describe long pulses

  size_t              fMaxMultiHit;              ///<maximum hits for multi fit
  int                 fAreaMethod;               ///<Type of area calculation
  std::vector<double> fAreaNorms;                ///<factors for converting area to same units as peak height
  bool              fTryNplus1Fits;            ///<whether we will (true) or won't (false) try n+1 fits
  double              fChi2NDFRetry;             ///<Value at which a second n+1 Fit will be tried
  double              fChi2NDF;                  ///maximum Chisquared / NDF allowed for a hit to be saved
  size_t              fNumBinsToAverage;         ///< If bin averaging for peak finding, number bins to average

  std::unique_ptr<HitFilterAlg> fHitFilterAlg;   ///algorithm used to filter out noise hits

  TH1F* fFirstChi2;
  TH1F* fChi2;



  void findCandidatePeaks(std::vector<float>::const_iterator startItr,
                          std::vector<float>::const_iterator stopItr,
                          TimeValsVec&                       timeValsVec,
                          float&                             roiThreshold,
                          int                                firstTick) const;

  void mergeCandidatePeaks(const std::vector<float>&,
                           TimeValsVec&,
                           MergedTimeWidVec&) const;


  void FitGaussians(const std::vector<float>& SignalVector,
                    const PeakTimeWidVec&     PeakVals,
                    int                       StartTime,
                    int                       EndTime,
                    double                    ampScaleFctr,
                    ParameterVec&             paramVec,
                    double&                   chi2PerNDF,
                    int&                      NDF);


  void doBinAverage(const std::vector<float>& inputVec,
                    std::vector<float>&       outputVec,
                    size_t                    binsToAverage) const;

  void reBin(const std::vector<float>& inputVec,
             std::vector<float>&       outputVec,
             size_t                    nBinsToCombine) const;



  void FillOutHitParameterVector(const std::vector<double>& input,
                                 std::vector<double>& output);




};

}

#endif // GAUS_HIT_FINDER_ALG_H
