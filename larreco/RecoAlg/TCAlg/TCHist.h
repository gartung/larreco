////////////////////////////////////////////////////////////////////////
//
//
// TCAlg debug struct
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERALGHISTSTRUCT_H
#define TRAJCLUSTERALGHISTSTRUCT_H

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

namespace tca {
  
  struct HistStuff {
    void CreateHists(art::TFileService& tfs);
    
    TH1F *fUnMatchedHitFrac;
    
    // True kinetic energy (MeV)
    TH1F *fTruT[5];
    
    TH2F *fMCSMom_TruMom_e;
    TH2F *fMCSMom_TruMom_mu;
    TH2F *fMCSMom_TruMom_pi;
    TH2F *fMCSMom_TruMom_p;
    
    TH2F *fMCSMomEP_TruMom_e;
    
    // Reco-MC vertex position difference
    TH1F* fNuVtx_dx;
    TH1F* fNuVtx_dy;
    TH1F* fNuVtx_dz;

    // Vertex score for 2D vertices that are near the neutrino interaction vertex
    TH1F* fNuVx3Score;
    TH1F* fNuVx2Score;
    TH1F* fNuVx3ScoreDiff;
    TH1F* fVxTopoMat;
    TH1F* fVxTopoNoMat;
    // Vertex score for 2D and 3D vertices
    TH1F* fVx2Score;
    TH1F* fVx3Score;
    
    // Reco-MC stopping wire difference for different MC Particles
    TH1F* fdWire[5];
    // EP vs KE for different MC Particles
    TProfile* fEP_T[5];
    // Trajectory proton likelihood
    TH1F* fProtonLike[5];
    TProfile* fProtonLike_T[5];
    
    // PFParticle PDGCode vs true PDG code
    TH2F* PDGCode_reco_true;
    
    
  };
} // namespace tca

#endif // ifndef TRAJCLUSTERALGHISTSTRUCT_H
