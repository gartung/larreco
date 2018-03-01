////////////////////////////////////////////////////////////////////////
//
//
// TrajClusterAlg
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERALG_H
#define TRAJCLUSTERALG_H

#include "larreco/RecoAlg/TCAlg/Utils.h"
#include "larreco/RecoAlg/TCAlg/TCTruth.h"
#include "larreco/RecoAlg/TCAlg/TCCR.h"

// C/C++ standard libraries
#include <array>
#include <vector>
#include <utility> // std::pair<>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/Exception.h"

// LArSoft libraries
#include "larreco/RecoAlg/LinFitAlg.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"

namespace tca {
  
  class TrajClusterAlg {
    public:
    
 
    /// @{
    /// @name Data structures for the reconstruction results


    TrajClusterAlg(fhicl::ParameterSet const& pset);

    virtual void reconfigure(fhicl::ParameterSet const& pset);

    void RunTrajClusterAlg(const art::Event & evt);

    void DefineShTree(TTree* t);
    
    void DefineCRTree(TTree* t);

    TjStuff const& GetTJS() const {return tjs;}

//    std::vector<short> const& GetinClus() const {return tjs.inClus; }
    
    /// Returns (and loses) the art::Ptr collection of previously reconstructed hits (e.g. gaushit)
    std::vector<recob::Hit> YieldHits();
    art::InputTag const& GetHitFinderModuleLabel() { return fHitFinderModuleLabel; }
    
    /// Returns a constant reference to the clusters found
    std::vector<ClusterStore> const& GetClusters() const { return tjs.tcl; }
    
    /// Returns a constant reference to the 2D end points found
    std::vector<VtxStore> const& GetEndPoints() const { return tjs.vtx; }
    
    /// Returns a constant reference to the 3D vertices found
    std::vector<Vtx3Store> const& GetVertices() const { return tjs.vtx3; }
    
    // Get the index list of matchVec entries that have PFParticle info defined
    std::vector<PFPStruct> const& GetPFParticles() const { return tjs.pfps; }
    // Get the cluster index of a trajectory ID
    unsigned int GetTjClusterIndex(unsigned int TjID) { return tjs.allTraj[TjID - 1].ClusterIndex; }
    // Get a ShowerStuct3D entry
    unsigned short GetShowerStructSize() { return tjs.showers.size(); };
    ShowerStruct3D const& GetShowerStruct(unsigned short ish) { return tjs.showers[ish]; };
    
    std::vector<unsigned int> const& GetAlgModCount() const {return fAlgModCount; }
    std::vector<std::string> const& GetAlgBitNames() const {return AlgBitNames; }
    
    /// Deletes all the results
    void ClearResults();
    
    TjStuff tjs;
    HistStuff hist;
    TruthMatcher tm{tjs};
    
    private:
    
    art::InputTag fHitFinderModuleLabel; ///< label of module producing input hits
    art::InputTag fHitTruthModuleLabel; ///< label of module producing MCParticle -> hit associations
    
    short fMode;            ///  StepCrawl mode (0 = turn off)
    std::vector<unsigned short> fMinPtsFit; ///< Reconstruct in two passes
    std::vector<unsigned short> fMinPts;    ///< min number of Pts required to make a trajectory
    std::vector<unsigned short> fMaxAngleCode;   ///< max allowed angle code for each pass
    std::vector<short> fMinMCSMom;   ///< Min MCSMom for each pass
    float fMultHitSep;      ///< preferentially "merge" hits with < this separation
    float fMaxChi;
    std::vector<float> fQualityCuts; ///< Min points/wire, min consecutive pts after a gap
    float fMaxWireSkipNoSignal;    ///< max number of wires to skip w/o a signal on them
    float fMaxWireSkipWithSignal;  ///< max number of wires to skip with a signal on them
    float fProjectionErrFactor;
    bool fMakeNewHits;
    float fJTMaxHitSep2;  /// Max hit separation for making junk trajectories. < 0 to turn off
    
    bool fTagAllTraj;              ///< tag clusters as shower-like or track-like
    float fMaxTrajSep;     ///< max trajectory point separation for making showers
    bool fStudyMode;       ///< study cuts
    bool fUseOldBackTracker {false};

    float fHitErrFac;   ///< hit time error = fHitErrFac * hit RMS used for cluster fit
    float fMinAmp;      ///< min amplitude required for declaring a wire signal is present
    float fVLAStepSize;
    
    float fLAClusSlopeCut;
    unsigned short fAllowNoHitWire;
    std::vector<float> fChkStopCuts; ///< [Min Chg ratio, Chg slope pull cut, Chg fit chi cut]
    
    bool fIsRealData;


    TH2F *fMCSMom_TruMom_e;
    TH2F *fMCSMom_TruMom_mu;
    TH2F *fMCSMom_TruMom_pi;
    TH2F *fMCSMom_TruMom_p;

    TH2F *fMCSMomEP_TruMom_e;
    
    // Reco-MC vertex position difference
    TH1F* fNuVtx_dx;
    TH1F* fNuVtx_dy;
    TH1F* fNuVtx_dz;
    TH1F* fNuVtx_Score;
    TProfile* fNuVtx_Enu_Score_p;
    
    // Vertex score for 2D and 3D vertices
    TH1F* fVx2_Score;
    TH1F* fVx3_Score;
    
    // Reco-MC stopping wire difference for different MC Particles
    TH1F* fdWire[5];
    // EP vs KE for different MC Particles
    TProfile* fEP_T[5];

    // SHOWER VARIABLE TREE
    TTree* showertree;

    // Cosmic Removal Variable Tree
    TTree* crtree;

    // number of primary particles in the event
    unsigned short nTruPrimary;
    float fNeutrinoEnergy;
    float fSourceParticleEnergy; //< in MeV
    // number of reconstructable primary particles in the event
    unsigned short nTruPrimaryOK;
    // number of reconstructable neutrino vertices in ALL events
    unsigned short nTruPrimaryVtxOK;
    // number of reconstructable neutrino vertices in ALL events that were reconstructed
    unsigned short nTruPrimaryVtxReco;

 
    bool prt;
    bool mrgPrt;
    bool didPrt;
    int TJPrt; // Set to the WorkID of a trajectory that is being debugged
    bool fDebugMode;
    
    trkf::LinFitAlg fLinFitAlg;
    calo::CalorimetryAlg fCaloAlg;

    int fWorkID;


    std::string fhitsModuleLabel;
    
    // variables for step crawling - updated at each TP
    // tracking flags
    bool fGoodTraj;         // the working trajectory is good and should be stored
    bool fTryWithNextPass;  // Try with next pass settings
    bool fMaskedLastTP;
    bool fQuitAlg;          // A significant error occurred. Delete everything and return
    
    std::vector<unsigned int> fAlgModCount;

//    bool fKalmanFilterFit;
//    trkf::TrackStatePropagator prop;
//    trkf::TrackKalmanFitter kalmanFitter;

    static bool SortByMultiplet(TCHit const& a, TCHit const& b);

//    short watchInTraj;
    // runs the TrajCluster algorithm on one plane specified by the calling routine
    void RunStepCrawl();
    void ReconstructAllTraj(CTP_t inCTP);
    // Main stepping/crawling routine
    void StepCrawl(Trajectory& tj);
    // Add hits on the trajectory point ipt that are close to the trajectory point Pos
    void AddHits(Trajectory& tj, unsigned short ipt, bool& sigOK);
    // Large Angle version
    void AddLAHits(Trajectory& tj, unsigned short ipt, bool& sigOK);
    // Try to use unused nearby hits in all trajectories after stepping is done
    void UseUnusedHits();
    // Finds junk trajectories using unassigned hits
    void FindJunkTraj(CTP_t inCTP);
    // Finds junk trajectories using unassigned hits
    bool MakeJunkTraj(std::vector<unsigned int> tHits, bool prt);
    // Step through TPs starting at the end and moving to the beginning
    void ReversePropagate(Trajectory& tj);
    // Start a trajectory going from fromHit to (toWire, toTick)
    bool StartTraj(Trajectory& tj, const unsigned int& fromHit, const unsigned int& toHit, const unsigned short& pass);
    bool StartTraj(Trajectory& tj, const float& fromWire, const float& fromTick, const float& toWire, const float& toTick, const CTP_t& tCTP, const unsigned short& pass);
    void GetHitMultiplet(unsigned int theHit, std::vector<unsigned int>& hitsInMultiplet);
    void GetHitMultiplet(unsigned int theHit, std::vector<unsigned int>& hitsInMultiplet, unsigned short& localIndex);
    // Returns fHits[iht]->RMS() * fScaleF * fHitErrFac * fHits[iht]->Multiplicity();
    float HitTimeErr(const unsigned int iht);
    // Estimates the error^2 of the time using all hits in hitVec
    float HitsTimeErr2(const std::vector<unsigned int>& hitVec);
    // defines HitPos, HitPosErr2 and Chg for the used hits in the trajectory point
    void ChkStopEndPts(Trajectory& tj, bool prt);
    void DefineHitPos(TrajPoint& tp);
    // Decide which hits to use to determine the trajectory point
    // fit, charge, etc. This is done by setting UseHit true and
    // setting inTraj < 0
     void FindUseHits(Trajectory& tj, unsigned short ipt, float maxDelta, bool useChg);
    // Print debug output if hit iht exists in work or allTraj
    void FindHit(std::string someText, unsigned int iht);
    // Check allTraj -> inTraj associations
    void ChkInTraj(std::string someText);
     // Make clusters from all trajectories in allTraj
    void MakeAllTrajClusters();
    void FindMissedVxTjs(const geo::TPCID& tpcid);
    void FindVtxTjs(CTP_t inCTP);
    void FindVtxTraj(VtxStore& theVtx);
    // Check the quality of the trajectory and possibly trim it
    void CheckTraj(Trajectory& tj);
     // Truncates the trajectory if a soft kink is found in it
    void FindSoftKink(Trajectory& tj);
    // Fill gaps in the trajectory
    void FillGaps(Trajectory& tj);
    // Check for many unused hits and try to use them
    void CheckHiMultUnusedHits(Trajectory& tj);
    void CheckHiMultEndHits(Trajectory& tj);
    // Check for high values of Delta at the beginning of the trajectory
    void HiEndDelta(Trajectory& tj);
    // Check for a TJ that is close to the Large Angle cut
    void CheckNearLA();
    // Updates the last added trajectory point fit, average hit rms, etc.
    void UpdateTraj(Trajectory& tj);
   // Estimate the Delta RMS of the TPs on the end of tj.
    void UpdateDeltaRMS(Trajectory& tj);
    void MaskBadTPs(Trajectory& tj, float const& maxChi);
    // The hits in the TP at the end of the trajectory were masked off. Decide whether to continue stepping with the
    // current configuration or whether to stop and possibly try with the next pass settings
    bool MaskedHitsOK(Trajectory& tj);
    // Any re-sizing should have been done by the calling routine. This code updates the Pass and adjusts the number of
    // fitted points to get FitCHi < 2
    bool StopIfBadFits(Trajectory& tj);
    void PrepareForNextPass(Trajectory& tj);
    // Does a local fit of just-added TPs to identify a kink while stepping.
    // Truncates the vector and returns true if one is found.
    void GottaKink(Trajectory& tj, unsigned short& killPts);
    // Update the parameters at the beginning of the trajectory
    void FixTrajBegin(Trajectory& tj);
    void FixTrajBegin(Trajectory& tj, unsigned short atPt);
    void FixTrajEnd(Trajectory& tj, unsigned short atPt);
    bool IsGhost(std::vector<unsigned int>& tHits, unsigned short& ofTraj);
    bool IsGhost(Trajectory& tj);
    void CheckTrajEnd();
    void EndMerge(CTP_t inCTP, bool lastPass);
    // Erases delHit and makes corrections to inTraj, allTraj and WireHitRange
    bool EraseHit(const unsigned int& delHit);
    // Creates a hit in tjs.fHits using the supplied information. Returns UINT_MAX if there is failure.
    // Returns the index of the newly created hit
    void DefineHit(TCHit& tcHit, CTP_t& hitCTP, unsigned int& hitWire);
    unsigned int CreateHit(TCHit tcHit);
    // Merges all of the hits used in each TP into one hit
    void MergeTPHits();
    void MaskTrajEndPoints(Trajectory& tj, unsigned short nPts);
    // Sets the StopsAtEnd bits for the trajectory
    void ChkStop(Trajectory& tj);
    void SplitTrajCrossingVertices();
    // Check the Michel electron topology, lastGoodPt is the last point of muon
    bool ChkMichel(Trajectory& tj, unsigned short& lastGoodPt);
    // TY: Split high charge hits near the trajectory end
    void ChkHiChgHits(CTP_t inCTP);
    void SplitHiChgHits(Trajectory& tj);

    void GetHitCollection(const art::Event& evt);

  }; // class TrajClusterAlg

} // namespace cluster

#endif // ifndef TRAJCLUSTERALG_H
