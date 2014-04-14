////////////////////////////////////////////////////////////////////////
// ClusterCrawlerAlg.h
//
// ClusterCrawlerAlg class
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef CCHITFINDERALG_H
#define CCHITFINDERALG_H

#include "TMath.h"

#include <vector>

#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 

#include "Geometry/Geometry.h"
#include "RecoBase/Hit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"

//namespace recob { class Hit; }

namespace cluster {

  class CCHitFinderAlg {
  
  public:
    
    struct CCHit {
      float Charge;
      float ChargeErr;
      float Amplitude;
      float AmplitudeErr;
      float Time;
      float TimeErr;
      float RMS;
      float RMSErr;
      float ChiDOF;
      art::Ptr<recob::Wire> Wire;
      unsigned short WireNum;
      unsigned short numHits;
      unsigned short LoHitID;
      float LoTime;   // defines the Lo(Hi) time region of the hit
      float HiTime;   // or hit multiplet
      short InClus;
    };
    std::vector< CCHit > allhits;
    
    // struct for passing hit fitting cuts to ClusterCrawler
    struct HitCuts {
      float MinSigInd;
      float MinSigCol;
      float MinRMSInd;
      float MinRMSCol;
      float ChiSplit;
      std::vector<float> ChiNorms;
      std::vector<float> TimeOffsets;
      std::vector<float> ChgNorms;
    };
    HitCuts hitcuts;

    CCHitFinderAlg(fhicl::ParameterSet const& pset);
    virtual ~CCHitFinderAlg();

    void reconfigure(fhicl::ParameterSet const& pset);

    void RunCCHitFinder(art::Event & evt);
    
  private:
    
    std::string     fCalDataModuleLabel;
    float fMinSigInd;     ///<Induction signal height threshold 
    float fMinSigCol;     ///<Collection signal height threshold 
    float fMinRMSInd;      ///<Initial rms for induction fit
    float fMinRMSCol;      ///<Initial rms for collection fit
    unsigned short fMaxBumps; // make a crude hit if > MaxBumps are found in the RAT
    unsigned short fMaxXtraHits; // max num of hits in Region Above Threshold
    float fChiSplit;      ///<Estimated noise error on the Signal
    float ChgNorm;     // Area norm for the wire we are working on

    std::vector<float> fChiNorms;
    std::vector<float> fTimeOffsets;
    std::vector<float> fChgNorms;

    uint32_t theChannel;
    unsigned short theWireNum;
    unsigned short thePlane;
    float minRMS;
    float minSig;
    float chinorm;
    float timeoff;
    const float Sqrt2Pi = 2.5066;
    const float SqrtPi  = 1.7725;


//    bool prt;
    
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;
    art::ServiceHandle<util::DetectorProperties> detprop;

    // fit n Gaussians possibly with bounds setting (parmin, parmax)
    void FitNG(unsigned short nGaus, unsigned short npt, float *ticks,
       float *signl);
    // parameters, errors, lower limit, upper limits for FitNG
    std::vector<double> par;
    std::vector<double> parerr;
    std::vector<double> parmin;
    std::vector<double> parmax;
    float chidof;
    std::vector<unsigned short> bumps;
    
    // make a cruddy hit if fitting fails
    void MakeCrudeHit(unsigned short npt, float *ticks, float *signl);
    // store the hits
    void StoreHits(unsigned short TStart, unsigned short npt, 
      art::Ptr<recob::Wire>& theWire);
/*
    // study hit finding and fitting
    void StudyHits(unsigned short flag, unsigned short npt = 0,
      float *ticks = 0, float *signl = 0, unsigned short tstart = 0);
    std::vector<int> bumpCnt;
    std::vector<int> RATCnt;
    std::vector<float> bumpChi;
    std::vector<float> bumpRMS;
    std::vector<int> hitCnt;
    std::vector<float> hitRMS;
    // use to determine the slope of protons
    std::vector<float> loWire;
    std::vector<float> loTime;
    std::vector<float> hiWire;
    std::vector<float> hiTime;
    bool SelRAT; // set true if a Region Above Threshold should be studied
*/

  }; // class CCHitFinderAlg

} // cchit

#endif // ifndef CCHITFINDERALG_H
