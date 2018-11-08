////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       ShowerValidation                                                                          //
// Modlue Type: Analyser                                                                                  //
// File         ShowerValidation_module.cc                                                                //
// Author:      Dominic Barker dominic.barker@sheffield.ac.uk                                             //
//                                                                                                        //
// Usage:       The Validation module takes an array of shower modules and perfoms truth matching on      //
//              reconstructed EM showers. It calculates the truth information from geant information      //
//              such as the diretion and starting position and energy deposited. Metrics are then defined //
//              to compare the efficiency of the shower reconstruction modules. These are presented       //
//              as histograms and in terms of energy. Energy plots look at the metrics for specific       //
//              energies between +- fEnergyWidth of the values described in the fcl table                 //
//                                                                                                        //
// Updates:     31.10.2018  Clustering Validation Added                                                   //
//                                                                                                        //
////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Framework includes                                                                
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//Larsoft includes 
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/MCRecoUtils/RecoUtils.h"
#include "larreco/RecoAlg/MCRecoUtils/ShowerUtils.h"

//Root Includes
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

//C++ Includes 
#include <vector>
#include <iostream>


namespace ana {
  class ShowerValidation;
}

class ana::ShowerValidation : public art::EDAnalyzer {
public:

  ShowerValidation(const fhicl::ParameterSet& pset);

  void analyze(const art::Event& evt);
  void endJob();
  void beginJob();

  void ClusterValidation(std::vector<art::Ptr<recob::Cluster> >& clusters, const art::Event& evt, art::Handle<std::vector<recob::Cluster> >& clusterHandle,  std::map<int,std::vector<int> >& ShowerMotherTrackIDs,std::map<int,float>& MCTrack_Energy_map, std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > > & MCTrack_hit_map, int& TrueShowerID, float& simenergy, std::string & fShowerModuleLabel);


private:

  //fcl parameters 
  std::string fGenieGenModuleLabel;
  std::string fLArGeantModuleLabel;
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;

  bool fUseBiggestShower;
  int  fVerbose;
  float fEnergyWidth;

  std::vector<float>       fEnergies;
  std::vector<std::string> fShowerModuleLabels;
  
  std::map<std::string,TH1F*> ShowerDirection_X_HistMap;
  std::map<std::string,TH1F*> ShowerDirection_Y_HistMap;
  std::map<std::string,TH1F*> ShowerDirection_Z_HistMap;
  std::map<std::string,TH1F*> ShowerStart_X_HistMap;
  std::map<std::string,TH1F*> ShowerStart_Y_HistMap;
  std::map<std::string,TH1F*> ShowerStart_Z_HistMap;
  std::map<std::string,TH1F*> ShowerLength_HistMap;
  std::map<std::string,TH1F*> ShowerEnergyDiff_HistMap;
  std::map<std::string,TH1F*> ShowerdEdx_HistMap;
  std::map<std::string,TH1F*> EventSeggy_HistMap;
  std::map<std::string,TH1F*> ShowerCompleteness_HistMap;
  std::map<std::string,TH1F*> ShowerPurity_HistMap;
  std::map<std::string,TH1F*> ShowerEnergy_HistMap;
  std::map<std::string,TH1F*> ShowerHitNum_HistMap;
  std::map<std::string,TH1F*> ShowerTotalEnergyDiff_HistMap;
  std::map<std::string,TH1F*> HitEnergy_HistMap;
  std::map<std::string,TH1F*> ShowerMag_HistMap;
  std::map<std::string,TH1F*> ShowerDirectionDiff_HistMap;
  std::map<std::string,TH1F*> ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap;
  std::map<std::string,TH1F*> TrueEnergy_HistMap;

  std::map<std::string,std::map<geo::PlaneID,TH1F*> > ClusterProjectionMatchedEnergy_HistMap;
  std::map<std::string,std::map<geo::PlaneID,TH1F*> > ClusterCompletenessEnergy_HistMap;
  std::map<std::string,std::map<geo::PlaneID,TH1F*> > ClusterPurityEnergy_HistMap;
  std::map<std::string,std::map<geo::PlaneID,TH1F*> > ClusterCompletnessHits_HistMap;
  std::map<std::string,std::map<geo::PlaneID,TH1F*> > ClusterPurityHits_HistMap;
  std::map<std::string,std::map<geo::PlaneID,TH1F*> > ClusterCompPurityEnergy_HistMap;
  std::map<std::string,std::map<geo::PlaneID,TH1F*> > ClusterCompPurityHits_HistMap;

  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerStart_X_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerStart_Y_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerStart_Z_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerLength_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerEnergyDiff_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerdEdx_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_EventSeggy_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerCompleteness_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerPurity_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerEnergy_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerHitNum_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerTotalEnergyDiff_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_HitEnergy_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerMag_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerDirectionDiff_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_TrueEnergy_HistMap;
  
  std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > > Energies_ClusterProjectionMatchedEnergy_HistMap;
  std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > > Energies_ClusterCompletenessEnergy_HistMap;
  std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > > Energies_ClusterPurityEnergy_HistMap;
  std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > > Energies_ClusterCompletnessHits_HistMap;
  std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > > Energies_ClusterPurityHits_HistMap;
  std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > > Energies_ClusterCompPurityEnergy_HistMap;
  std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > > Energies_ClusterCompPurityHits_HistMap;
  

  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerStart_X_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerStart_Y_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerStart_Z_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerLength_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerEnergyDiff_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerdEdx_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_EventSeggy_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerCompleteness_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerPurity_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerEnergy_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerHitNum_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerTotalEnergyDiff_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_HitEnergy_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerMag_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerDirectionDiff_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_TrueEnergy_GraphMap;
  
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_Mean_ClusterProjectionMatchedEnergy_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_Mean_ClusterCompletenessEnergy_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_Mean_ClusterPurityEnergy_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_Mean_ClusterCompletnessHits_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_Mean_ClusterPurityHits_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_Mean_ClusterCompPurityEnergy_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_Mean_ClusterCompPurityHits_GraphMap;


  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerStart_X_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerStart_Y_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerStart_Z_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerLength_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerEnergyDiff_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerdEdx_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_EventSeggy_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerCompleteness_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerPurity_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerEnergy_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerHitNum_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerTotalEnergyDiff_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_HitEnergy_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerMag_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerDirectionDiff_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_TrueEnergy_GraphMap;

  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_RMS_ClusterProjectionMatchedEnergy_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_RMS_ClusterCompletenessEnergy_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_RMS_ClusterPurityEnergy_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_RMS_ClusterCompletnessHits_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_RMS_ClusterPurityHits_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_RMS_ClusterCompPurityEnergy_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_RMS_ClusterCompPurityHits_GraphMap;


  TMultiGraph* Energies_Mean_ShowerStart_X_Multi;
  TMultiGraph* Energies_Mean_ShowerStart_Y_Multi;
  TMultiGraph* Energies_Mean_ShowerStart_Z_Multi;
  TMultiGraph* Energies_Mean_ShowerLength_Multi;
  TMultiGraph* Energies_Mean_ShowerEnergyDiff_Multi;
  TMultiGraph* Energies_Mean_ShowerdEdx_Multi;
  TMultiGraph* Energies_Mean_EventSeggy_Multi;
  TMultiGraph* Energies_Mean_ShowerCompleteness_Multi;
  TMultiGraph* Energies_Mean_ShowerPurity_Multi;
  TMultiGraph* Energies_Mean_ShowerEnergy_Multi;
  TMultiGraph* Energies_Mean_ShowerHitNum_Multi;
  TMultiGraph* Energies_Mean_ShowerTotalEnergyDiff_Multi;
  TMultiGraph* Energies_Mean_HitEnergy_Multi;
  TMultiGraph* Energies_Mean_ShowerMag_Multi;
  TMultiGraph* Energies_Mean_ShowerDirectionDiff_Multi;
  TMultiGraph* Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_Multi;
  TMultiGraph* Energies_Mean_TrueEnergy_Multi;

  std::map<geo::PlaneID,TMultiGraph*> Energies_Mean_ClusterProjectionMatchedEnergy_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_Mean_ClusterCompletenessEnergy_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_Mean_ClusterPurityEnergy_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_Mean_ClusterCompletnessHits_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_Mean_ClusterPurityHits_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_Mean_ClusterCompPurityEnergy_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_Mean_ClusterCompPurityHits_Multi;
  

  TMultiGraph* Energies_RMS_ShowerStart_X_Multi;
  TMultiGraph* Energies_RMS_ShowerStart_Y_Multi;
  TMultiGraph* Energies_RMS_ShowerStart_Z_Multi;
  TMultiGraph* Energies_RMS_ShowerLength_Multi;
  TMultiGraph* Energies_RMS_ShowerEnergyDiff_Multi;
  TMultiGraph* Energies_RMS_ShowerdEdx_Multi;
  TMultiGraph* Energies_RMS_EventSeggy_Multi;
  TMultiGraph* Energies_RMS_ShowerCompleteness_Multi;
  TMultiGraph* Energies_RMS_ShowerPurity_Multi;
  TMultiGraph* Energies_RMS_ShowerEnergy_Multi;
  TMultiGraph* Energies_RMS_ShowerHitNum_Multi;
  TMultiGraph* Energies_RMS_ShowerTotalEnergyDiff_Multi;
  TMultiGraph* Energies_RMS_HitEnergy_Multi;
  TMultiGraph* Energies_RMS_ShowerMag_Multi;
  TMultiGraph* Energies_RMS_ShowerDirectionDiff_Multi;
  TMultiGraph* Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_Multi;
  TMultiGraph* Energies_RMS_TrueEnergy_Multi;

  std::map<geo::PlaneID,TMultiGraph*> Energies_RMS_ClusterProjectionMatchedEnergy_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_RMS_ClusterCompletenessEnergy_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_RMS_ClusterPurityEnergy_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_RMS_ClusterCompletnessHits_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_RMS_ClusterPurityHits_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_RMS_ClusterCompPurityEnergy_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_RMS_ClusterCompPurityHits_Multi;


  TCanvas* Energies_Mean_ShowerStart_X_canvasMulti;
  TCanvas* Energies_Mean_ShowerStart_Y_canvasMulti;
  TCanvas* Energies_Mean_ShowerStart_Z_canvasMulti;
  TCanvas* Energies_Mean_ShowerLength_canvasMulti;
  TCanvas* Energies_Mean_ShowerEnergyDiff_canvasMulti;
  TCanvas* Energies_Mean_ShowerTotalEnergyDiff_canvasMulti;
  TCanvas* Energies_Mean_ShowerdEdx_canvasMulti;
  TCanvas* Energies_Mean_EventSeggy_canvasMulti;
  TCanvas* Energies_Mean_ShowerCompleteness_canvasMulti;
  TCanvas* Energies_Mean_ShowerPurity_canvasMulti;
  TCanvas* Energies_Mean_ShowerEnergy_canvasMulti;
  TCanvas* Energies_Mean_ShowerHitNum_canvasMulti;
  TCanvas* Energies_Mean_HitEnergy_canvasMulti;
  TCanvas* Energies_Mean_ShowerDirectionDiff_canvasMulti;
  TCanvas* Energies_Mean_ShowerMag_canvasMulti;
  TCanvas* Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMulti;
  TCanvas* Energies_Mean_TrueEnergy_canvasMulti;

  std::map<geo::PlaneID,TCanvas*> Energies_Mean_ClusterProjectionMatchedEnergy_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_Mean_ClusterCompletenessEnergy_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_Mean_ClusterPurityEnergy_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_Mean_ClusterCompletnessHits_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_Mean_ClusterPurityHits_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_Mean_ClusterCompPurityEnergy_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_Mean_ClusterCompPurityHits_canvasMulti;


  TCanvas* Energies_RMS_ShowerStart_X_canvasMulti;
  TCanvas* Energies_RMS_ShowerStart_Y_canvasMulti;
  TCanvas* Energies_RMS_ShowerStart_Z_canvasMulti;
  TCanvas* Energies_RMS_ShowerLength_canvasMulti;
  TCanvas* Energies_RMS_ShowerEnergyDiff_canvasMulti;
  TCanvas* Energies_RMS_ShowerTotalEnergyDiff_canvasMulti;
  TCanvas* Energies_RMS_ShowerdEdx_canvasMulti;
  TCanvas* Energies_RMS_EventSeggy_canvasMulti;
  TCanvas* Energies_RMS_ShowerCompleteness_canvasMulti;
  TCanvas* Energies_RMS_ShowerPurity_canvasMulti;
  TCanvas* Energies_RMS_ShowerEnergy_canvasMulti;
  TCanvas* Energies_RMS_ShowerHitNum_canvasMulti;
  TCanvas* Energies_RMS_HitEnergy_canvasMulti;
  TCanvas* Energies_RMS_ShowerDirectionDiff_canvasMulti;
  TCanvas* Energies_RMS_ShowerMag_canvasMulti;
  TCanvas* Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMulti;
  TCanvas* Energies_RMS_TrueEnergy_canvasMulti;

  std::map<geo::PlaneID,TCanvas*> Energies_RMS_ClusterProjectionMatchedEnergy_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_RMS_ClusterCompletenessEnergy_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_RMS_ClusterPurityEnergy_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_RMS_ClusterCompletnessHits_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_RMS_ClusterPurityHits_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_RMS_ClusterCompPurityEnergy_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_RMS_ClusterCompPurityHits_canvasMulti;

  TCanvas* ShowerDirection_X_canvas;
  TCanvas* ShowerDirection_Y_canvas;
  TCanvas* ShowerDirection_Z_canvas;
  TCanvas* ShowerStart_X_canvas;
  TCanvas* ShowerStart_Y_canvas;
  TCanvas* ShowerStart_Z_canvas;
  TCanvas* ShowerLength_canvas;
  TCanvas* ShowerEnergyDiff_canvas;
  TCanvas* ShowerTotalEnergyDiff_canvas;
  TCanvas* ShowerdEdx_canvas;
  TCanvas* EventSeggy_canvas;
  TCanvas* ShowerCompleteness_canvas;
  TCanvas* ShowerPurity_canvas;
  TCanvas* ShowerEnergy_canvas;
  TCanvas* ShowerHitNum_canvas;
  TCanvas* HitEnergy_canvas;
  TCanvas* ShowerDirectionDiff_canvas;
  TCanvas* ShowerMag_canvas;
  TCanvas* ShowerRecoEnergyVsTrueEnergyinRecoShower_canvas;
  TCanvas* TrueEnergy_canvas;

  std::map<geo::PlaneID,TCanvas*> ClusterProjectionMatchedEnergy_canvas;
  std::map<geo::PlaneID,TCanvas*> ClusterCompletenessEnergy_canvas;
  std::map<geo::PlaneID,TCanvas*> ClusterPurityEnergy_canvas;
  std::map<geo::PlaneID,TCanvas*> ClusterCompletnessHits_canvas;
  std::map<geo::PlaneID,TCanvas*> ClusterPurityHits_canvas;
  std::map<geo::PlaneID,TCanvas*> ClusterCompPurityEnergy_canvas;
  std::map<geo::PlaneID,TCanvas*> ClusterCompPurityHits_canvas;
  

  std::map<float,TCanvas*> Energies_ShowerStart_X_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerStart_Y_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerStart_Z_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerLength_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerEnergyDiff_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerTotalEnergyDiff_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerdEdx_canvasMap;
  std::map<float,TCanvas*> Energies_EventSeggy_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerCompleteness_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerPurity_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerEnergy_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerHitNum_canvasMap;
  std::map<float,TCanvas*> Energies_HitEnergy_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerDirectionDiff_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerMag_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMap;
  std::map<float,TCanvas*> Energies_TrueEnergy_canvasMap;

  std::map<float,std::map<geo::PlaneID,TCanvas*> > Energies_ClusterProjectionMatchedEnergy_canvasMap;
  std::map<float,std::map<geo::PlaneID,TCanvas*> > Energies_ClusterCompletenessEnergy_canvasMap;
  std::map<float,std::map<geo::PlaneID,TCanvas*> > Energies_ClusterPurityEnergy_canvasMap;
  std::map<float,std::map<geo::PlaneID,TCanvas*> > Energies_ClusterCompletnessHits_canvasMap;
  std::map<float,std::map<geo::PlaneID,TCanvas*> > Energies_ClusterPurityHits_canvasMap;
  std::map<float,std::map<geo::PlaneID,TCanvas*> > Energies_ClusterCompPurityEnergy_canvasMap;
  std::map<float,std::map<geo::PlaneID,TCanvas*> > Energies_ClusterCompPurityHits_canvasMap;

  //Service handles   
  art::ServiceHandle<cheat::BackTrackerService> backtracker;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
  art::ServiceHandle<geo::Geometry> geom;

};

ana::ShowerValidation::ShowerValidation(const fhicl::ParameterSet& pset) : EDAnalyzer(pset){

  fGenieGenModuleLabel   = pset.get<std::string>("GenieGenModuleLabel");
  fLArGeantModuleLabel   = pset.get<std::string>("LArGeantModuleLabel"); 
  fHitsModuleLabel       = pset.get<std::string>("HitsModuleLabel");
  fTrackModuleLabel      = pset.get<std::string>("TrackModuleLabel");
  fShowerModuleLabels    = pset.get<std::vector<std::string> >("ShowerModuleLabels");
  fEnergies              = pset.get<std::vector<float> >("Energies");
  fUseBiggestShower      = pset.get<bool>("UseBiggestShower");
  fVerbose               = pset.get<int>("Verbose");
  fEnergyWidth           = pset.get<float>("EnergyWidth");
  
}
  
void ana::ShowerValidation::beginJob() {

  gStyle->SetOptStat(0);

  art::ServiceHandle<art::TFileService> tfs;

  //Create the canvases for different energyies 
  for(unsigned int i=0; i<fEnergies.size(); ++i){

    //Canvas strings
    std::string Energies_ShowerStart_X_canvasMap_string     = "Energies_ShowerStart_X_canvas at Energy " + std::to_string(fEnergies[i]); 
    std::string Energies_ShowerStart_Y_canvasMap_string    = "Energies_ShowerStart_Y_canvas at Energy " + std::to_string(fEnergies[i]);
    std::string Energies_ShowerStart_Z_canvasMap_string    = "Energies_ShowerStart_Z_canvas at Energy " + std::to_string(fEnergies[i]);
    std::string Energies_ShowerLength_canvasMap_string     = "Energies_ShowerLength_canvas at Energy " + std::to_string(fEnergies[i]);
    std::string Energies_ShowerEnergyDiff_canvasMap_string = "Energies_ShowerEnergyDiff_canvas at Energy " + std::to_string(fEnergies[i]);
    std::string Energies_ShowerdEdx_canvasMap_string = "Energies_ShowerdEdx_canvas at Energy " + std::to_string(fEnergies[i]);
    std::string Energies_EventSeggy_canvasMap_string = "Energies_EventSeggy_canvas at Energy " + std::to_string(fEnergies[i]);
    std::string Energies_ShowerCompleteness_canvasMap_string = "Energies_ShowerCompleteness_canvas at Energy " + std::to_string(fEnergies[i]);
    std::string Energies_ShowerPurity_canvasMap_string = "Energies_ShowerPurity_canvas at Energy " + std::to_string(fEnergies[i]);
    std::string Energies_ShowerEnergy_canvasMap_string = "Energies_ShowerEnergy_canvas at Energy " + std::to_string(fEnergies[i]);
    std::string Energies_ShowerHitNum_canvasMap_string = "Energies_ShowerHitNum_canvas at Energy " + std::to_string(fEnergies[i]);
    std::string Energies_HitEnergy_canvasMap_string = "Energies_HitEnergy_canvas at Energy " + std::to_string(fEnergies[i]); 
    std::string Energies_ShowerDirectionDiff_canvasMap_string = "Energies_ShowerDirectionDiff_canvas at Energy " + std::to_string(fEnergies[i]);
    std::string Energies_ShowerMag_canvasMap_string = "Energies_ShowerMag_canvas at Energy " + std::to_string(fEnergies[i]);
    std::string Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMap_string = "Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvas at Energy " + std::to_string(fEnergies[i]);
    std::string Energies_TrueEnergy_canvasMap_string = "Energies_TrueEnergy_canvasMap_canvas at Energy " + std::to_string(fEnergies[i]); 
    
     const char* Energies_ShowerStart_X_canvasMap_name = Energies_ShowerStart_X_canvasMap_string.c_str();
     const char* Energies_ShowerStart_Y_canvasMap_name = Energies_ShowerStart_Y_canvasMap_string.c_str();
     const char* Energies_ShowerStart_Z_canvasMap_name = Energies_ShowerStart_Z_canvasMap_string.c_str();
     const char* Energies_ShowerLength_canvasMap_name = Energies_ShowerLength_canvasMap_string.c_str();
     const char* Energies_ShowerEnergyDiff_canvasMap_name = Energies_ShowerEnergyDiff_canvasMap_string.c_str();
     const char* Energies_ShowerdEdx_canvasMap_name = Energies_ShowerdEdx_canvasMap_string.c_str();
     const char* Energies_EventSeggy_canvasMap_name = Energies_EventSeggy_canvasMap_string.c_str();
     const char* Energies_ShowerCompleteness_canvasMap_name = Energies_ShowerCompleteness_canvasMap_string.c_str();
     const char* Energies_ShowerPurity_canvasMap_name = Energies_ShowerPurity_canvasMap_string.c_str();
     const char* Energies_ShowerEnergy_canvasMap_name = Energies_ShowerEnergy_canvasMap_string.c_str();
     const char* Energies_ShowerHitNum_canvasMap_name = Energies_ShowerHitNum_canvasMap_string.c_str();
     const char* Energies_HitEnergy_canvasMap_name = Energies_HitEnergy_canvasMap_string.c_str();
     const char* Energies_ShowerDirectionDiff_canvasMap_name = Energies_ShowerDirectionDiff_canvasMap_string.c_str();
     const char* Energies_ShowerMag_canvasMap_name = Energies_ShowerMag_canvasMap_string.c_str();
     const char* Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMap_name = Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMap_string.c_str();
     const char* Energies_TrueEnergy_canvasMap_name = Energies_TrueEnergy_canvasMap_string.c_str();
    
     //Create the TCanvas 
     Energies_ShowerStart_X_canvasMap[fEnergies[i]] = tfs->makeAndRegister<TCanvas>(Energies_ShowerStart_X_canvasMap_name,";|True X Position-Reco X Postition| (cm);Enteries");
     Energies_ShowerStart_Y_canvasMap[fEnergies[i]] = tfs->makeAndRegister<TCanvas>(Energies_ShowerStart_Y_canvasMap_name,";|True Y Position-Reco Y Postition| (cm);Enteries");
     Energies_ShowerStart_Z_canvasMap[fEnergies[i]] = tfs->makeAndRegister<TCanvas>(Energies_ShowerStart_Z_canvasMap_name,";|True X Position-Reco X Postition|(cm);Enteries");
     Energies_ShowerLength_canvasMap[fEnergies[i]] = tfs->makeAndRegister<TCanvas>(Energies_ShowerLength_canvasMap_name,";|True Shower Length - Reco Shower Length| (cm); Enteries");
     Energies_ShowerEnergyDiff_canvasMap[fEnergies[i]] = tfs->makeAndRegister<TCanvas>(Energies_ShowerEnergyDiff_canvasMap_name,";(True Energy - Recon Energy)/TrueEnergy;Enteries");
     Energies_ShowerdEdx_canvasMap[fEnergies[i]] = tfs->makeAndRegister<TCanvas>(Energies_ShowerdEdx_canvasMap_name,";dEdx MeV/cm;Enteries");
     Energies_EventSeggy_canvasMap[fEnergies[i]] = tfs->makeAndRegister<TCanvas>(Energies_EventSeggy_canvasMap_name,";Number of Showers;Enteries");
     Energies_ShowerCompleteness_canvasMap[fEnergies[i]] = tfs->makeAndRegister<TCanvas>(Energies_ShowerCompleteness_canvasMap_name,";Completeness;Enteries");
     Energies_ShowerPurity_canvasMap[fEnergies[i]] = tfs->makeAndRegister<TCanvas>(Energies_ShowerPurity_canvasMap_name,";Purity;Enteries;");
     Energies_ShowerEnergy_canvasMap[fEnergies[i]] = tfs->makeAndRegister<TCanvas>(Energies_ShowerEnergy_canvasMap_name,";Energy (MeV);Enteries");
     Energies_ShowerHitNum_canvasMap[fEnergies[i]] = tfs->makeAndRegister<TCanvas>(Energies_ShowerHitNum_canvasMap_name,";Hits;Enteries");
     Energies_HitEnergy_canvasMap[fEnergies[i]] = tfs->makeAndRegister<TCanvas>(Energies_HitEnergy_canvasMap_name,"Energy (MeV)");
     Energies_ShowerDirectionDiff_canvasMap[fEnergies[i]] = tfs->makeAndRegister<TCanvas>(Energies_ShowerDirectionDiff_canvasMap_name,";Angle Between Directions (rad);Enteries");
     Energies_ShowerMag_canvasMap[fEnergies[i]] =tfs->makeAndRegister<TCanvas>(Energies_ShowerMag_canvasMap_name,";|True Start - Reco Start| (cm)","Enteries");
     Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMap[fEnergies[i]]= tfs->makeAndRegister<TCanvas>(Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMap_name,";ShowerRecoEnergy/True Energy in Shower; Enteries");
     Energies_TrueEnergy_canvasMap[fEnergies[i]]= tfs->makeAndRegister<TCanvas>(Energies_TrueEnergy_canvasMap_name,";True Energy;Enteries");


     for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){
       std::string Energies_ClusterProjectionMatchedEnergy_canvasMap_string = "Energies_ClusterProjectionMatchedEnergy_canvas at Energy " + std::to_string(fEnergies[i]) + " For Plane: " + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + std::to_string(plane_id.Cryostat);
       std::string Energies_ClusterCompletenessEnergy_canvasMap_string = "Energies_ClusterCompletenessEnergy_canvas at Energy " + std::to_string(fEnergies[i]) + " For Plane: " + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + std::to_string(plane_id.Cryostat);
       std::string Energies_ClusterPurityEnergy_canvasMap_string = "Energies_ClusterPurityEnergy_canvas at Energy " + std::to_string(fEnergies[i]) + " For Plane: " + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + std::to_string(plane_id.Cryostat);
       std::string Energies_ClusterCompletnessHits_canvasMap_string = "Energies_ClusterCompletnessHits_canvas at Energy " + std::to_string(fEnergies[i]) + " For Plane: " + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + std::to_string(plane_id.Cryostat);
       std::string Energies_ClusterPurityHits_canvasMap_string = "Energies_ClusterPurityHits_canvas at Energy " + std::to_string(fEnergies[i]) + " For Plane: " + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + std::to_string(plane_id.Cryostat);
       std::string Energies_ClusterCompPurityEnergy_canvasMap_string = "Energies_ClusterCompPurityEnergy at Energy " + std::to_string(fEnergies[i]) + " For Plane: " + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + std::to_string(plane_id.Cryostat);
       std::string Energies_ClusterCompPurityHits_canvasMap_string = "Energies_ClusterCompPurityHits_canvas at Energy " + std::to_string(fEnergies[i]) + " For Plane: " + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + std::to_string(plane_id.Cryostat);
       
       const char* Energies_ClusterProjectionMatchedEnergy_canvasMap_name = Energies_ClusterProjectionMatchedEnergy_canvasMap_string.c_str();
       const char* Energies_ClusterCompletenessEnergy_canvasMap_name = Energies_ClusterCompletenessEnergy_canvasMap_string.c_str();
       const char* Energies_ClusterPurityEnergy_canvasMap_name = Energies_ClusterPurityEnergy_canvasMap_string.c_str(); 
       const char* Energies_ClusterCompletnessHits_canvasMap_name = Energies_ClusterCompletnessHits_canvasMap_string.c_str();
       const char* Energies_ClusterPurityHits_canvasMap_name = Energies_ClusterPurityHits_canvasMap_string.c_str(); 
       const char* Energies_ClusterCompPurityEnergy_canvasMap_name = Energies_ClusterCompPurityEnergy_canvasMap_string.c_str();
       const char* Energies_ClusterCompPurityHits_canvasMap_name = Energies_ClusterCompPurityHits_canvasMap_string.c_str();
       
       Energies_ClusterProjectionMatchedEnergy_canvasMap[fEnergies[i]][plane_id]= tfs->makeAndRegister<TCanvas>(Energies_ClusterProjectionMatchedEnergy_canvasMap_name,"Cluster Matched Correctly; Enteries");
       Energies_ClusterCompletenessEnergy_canvasMap[fEnergies[i]][plane_id]= tfs->makeAndRegister<TCanvas>(Energies_ClusterCompletenessEnergy_canvasMap_name,";Completeness;Enteries");
       Energies_ClusterPurityEnergy_canvasMap[fEnergies[i]][plane_id]= tfs->makeAndRegister<TCanvas>(Energies_ClusterPurityEnergy_canvasMap_name,";Purity;Enteries");
       Energies_ClusterCompletnessHits_canvasMap[fEnergies[i]][plane_id]= tfs->makeAndRegister<TCanvas>(Energies_ClusterCompletnessHits_canvasMap_name,";Completeness;Enteries");
       Energies_ClusterPurityHits_canvasMap[fEnergies[i]][plane_id]= tfs->makeAndRegister<TCanvas>(Energies_ClusterPurityHits_canvasMap_name,";Purity;Enteries");
       Energies_ClusterCompPurityEnergy_canvasMap[fEnergies[i]][plane_id]= tfs->makeAndRegister<TCanvas>(Energies_ClusterCompPurityEnergy_canvasMap_name,";Purity * Completeness;Enteries");
       Energies_ClusterCompPurityHits_canvasMap[fEnergies[i]][plane_id]= tfs->makeAndRegister<TCanvas>(Energies_ClusterCompPurityHits_canvasMap_name,";Purity * Completeness;Enteries");
              
     }
  }


  for(unsigned int j=0; j<fShowerModuleLabels.size(); ++j){
    
    std::string ShowerDirection_X_string = "ShowerDirection X " + fShowerModuleLabels[j];
    std::string ShowerDirection_Y_string = "ShowerDirection Y " + fShowerModuleLabels[j];
    std::string ShowerDirection_Z_string = "ShowerDirection Z " + fShowerModuleLabels[j];

    std::string ShowerStartX_string = "Shower Start X " + fShowerModuleLabels[j];
    std::string ShowerStartY_string = "Shower Start Y " + fShowerModuleLabels[j];
    std::string ShowerStartZ_string = "Shower Start Z " + fShowerModuleLabels[j];

    std::string ShowerLength_string          = "Shower Length " + fShowerModuleLabels[j];
    std::string ShowerEnergyDiff_string      =  "Shower Energy Diff " + fShowerModuleLabels[j];
    std::string ShowerdEdX_string            = "Shower dEdX initial track" + fShowerModuleLabels[j];
    std::string EventSeggy_string            = "Seggy " + fShowerModuleLabels[j];
    std::string ShowerCompleteness_string    = "Shower Completeness " + fShowerModuleLabels[j];
    std::string ShowerPurity_string          = "Shower Purity " + fShowerModuleLabels[j];
    std::string ShowerEnergy_string          = "Shower Energy " + fShowerModuleLabels[j];
    std::string HitEnergy_string             = "Shower Hit Energy " +  fShowerModuleLabels[j];
    std::string ShowerHitNum_string          = "Shower Hit Num " + fShowerModuleLabels[j];
    std::string EnergyRes_Mean_string        = "Energy Res Mean " + fShowerModuleLabels[j];
    std::string EnergyRes_RMS_string         = "Energy Res RMS " +  fShowerModuleLabels[j];
    std::string ShowerMag_string             = "Shower Mang Position " + fShowerModuleLabels[j];
    std::string ShowerDirectionDiff_string         = "Shower Direction Difference" +  fShowerModuleLabels[j];
    std::string ShowerRecoEnergyVsTrueEnergyinRecoShower_string = "Reco Energy/ True Energy in Reco Shower" + fShowerModuleLabels[j];
    std::string TrueEnergy_string = "True Energy" + fShowerModuleLabels[j];
 
    std::string ShowerDirection_X_titlestring      = ";X Direction;Enteries";
    std::string ShowerDirection_Y_titlestring      = ";Y Direction;Enteries";
    std::string ShowerDirection_Z_titlestring      = ";Z Direction;Enteries";
    std::string ShowerStartX_titlestring           = ";|True X Position-Reco X Postition(cm)|;Enteries";
    std::string ShowerStartY_titlestring           = ";|True Y Position-Reco Y Postition(cm)|;Enteries";
    std::string ShowerStartZ_titlestring           = ";|True Z Position-Reco Z Postition(cm)|;Enteries";
    std::string ShowerLength_titlestring           = ";|True Shower Length - Reco Shower Length| (cm); Enteries";
    std::string ShowerEnergyDiff_titlestring       = ";(Reco Energy)/TrueEnergy;Enteries";
    std::string ShowerdEdX_titlestring             = ";dEdx (MeV/cm);Enteries";
    std::string EventSeggy_titlestring             = ";Number of Showers;Enteries";
    std::string ShowerCompleteness_titlestring     = ";Completeness;Enteries";
    std::string ShowerPurity_titlestring           = ";Purity;Enteries;";
    std::string ShowerEnergy_titlestring           = ";Energy (MeV);Enteries";
    std::string HitEnergy_titlestring              = ";Energy (MeV);Enteries";
    std::string ShowerHitNum_titlestring           = "; Number of Hits;Enteries";
    std::string ShowerMag_titlestring              = ";|True Shower Start - RecoShower Start| (cm);Enteries";
    std::string ShowerDirectionDiff_titlestring     = ";(True Direction).(Reco Direction) (cm^2); Enteries"; 
    std::string ShowerRecoEnergyVsTrueEnergyinRecoShower_titlestring = "Reco Energy/ True Energy in Reco Shower; Enteries"; 
    std::string TrueEnergy_titlestring = "True Energy (MeV); Enteries"; 

    const char* ShowerDirection_X_titlename = ShowerDirection_X_titlestring.c_str();
    const char* ShowerDirection_Y_titlename = ShowerDirection_Y_titlestring.c_str(); 
    const char* ShowerDirection_Z_titlename = ShowerDirection_Z_titlestring.c_str();  
    const char* ShowerStartX_titlename = ShowerStartX_titlestring.c_str();
    const char* ShowerStartY_titlename = ShowerStartY_titlestring.c_str();
    const char* ShowerStartZ_titlename = ShowerStartZ_titlestring.c_str();
    const char* ShowerLength_titlename = ShowerLength_titlestring.c_str();
    const char* ShowerEnergyDiff_titlename = ShowerEnergyDiff_titlestring.c_str();
    const char* ShowerdEdX_titlename =  ShowerdEdX_titlestring.c_str();
    const char* EventSeggy_titlename = EventSeggy_titlestring.c_str();
    const char* ShowerCompleteness_titlename = ShowerCompleteness_titlestring.c_str();
    const char* ShowerPurity_titlename = ShowerPurity_titlestring.c_str();
    const char* ShowerEnergy_titlename = ShowerEnergy_titlestring.c_str();
    const char* HitEnergy_titlename = HitEnergy_titlestring.c_str();
    const char* ShowerHitNum_titlename = ShowerHitNum_titlestring.c_str();
    const char* ShowerMag_titlename         = ShowerMag_titlestring.c_str();
    const char* ShowerDirectionDiff_titlename   = ShowerDirectionDiff_titlestring.c_str();
    const char* ShowerRecoEnergyVsTrueEnergyinRecoShower_titlename = ShowerRecoEnergyVsTrueEnergyinRecoShower_titlestring.c_str();
    const char* TrueEnergy_titlename = TrueEnergy_titlestring.c_str();

    const char* ShowerDirection_X_name     =  ShowerDirection_X_string.c_str();
    const char* ShowerDirection_Y_name     =  ShowerDirection_Y_string.c_str();
    const char* ShowerDirection_Z_name     =  ShowerDirection_Z_string.c_str();
    const char* ShowerStartX_name          =  ShowerStartX_string.c_str();
    const char* ShowerStartY_name          =  ShowerStartY_string.c_str();
    const char* ShowerStartZ_name          =  ShowerStartZ_string.c_str();
    const char* ShowerLength_name          =  ShowerLength_string.c_str();
    const char* ShowerEnergyDiff_name      =  ShowerEnergyDiff_string.c_str();
    const char* ShowerdEdX_name            = ShowerdEdX_string.c_str();
    const char* EventSeggy_name            = EventSeggy_string.c_str();
    const char* ShowerCompleteness_name    = ShowerCompleteness_string.c_str();
    const char* ShowerPurity_name          = ShowerPurity_string.c_str();
    const char* ShowerEnergy_name          = ShowerEnergy_string.c_str();
    const char* ShowerHitNum_name          = ShowerHitNum_string.c_str();
    const char* HitEnergy_name             = HitEnergy_string.c_str();
    const char* ShowerMag_name             = ShowerMag_string.c_str();
    const char* ShowerDirectionDiff_name   = ShowerDirectionDiff_string.c_str(); 
    const char* ShowerRecoEnergyVsTrueEnergyinRecoShower_name = ShowerRecoEnergyVsTrueEnergyinRecoShower_string.c_str();
    const char* TrueEnergy_name = TrueEnergy_string.c_str();
 
    ShowerDirection_X_HistMap[fShowerModuleLabels[j]]  = new TH1F(ShowerDirection_X_name, ShowerDirection_X_titlename, 100, -1, 1);
    ShowerDirection_Y_HistMap[fShowerModuleLabels[j]]  = new TH1F(ShowerDirection_Y_name, ShowerDirection_Y_titlename, 100, -1, 1);
    ShowerDirection_Z_HistMap[fShowerModuleLabels[j]]  = new TH1F(ShowerDirection_Z_name, ShowerDirection_Z_titlename, 100, -1, 1);

    
    ShowerStart_X_HistMap[fShowerModuleLabels[j]]  = new TH1F(ShowerStartX_name, ShowerStartX_titlename, 400, 0, 20);
    ShowerStart_Y_HistMap[fShowerModuleLabels[j]]  = new TH1F(ShowerStartY_name, ShowerStartY_titlename, 400, 0, 20);
    ShowerStart_Z_HistMap[fShowerModuleLabels[j]]  = new TH1F(ShowerStartZ_name, ShowerStartZ_titlename, 600, 0, 20);
    
    ShowerLength_HistMap[fShowerModuleLabels[j]]          = new TH1F(ShowerLength_name, ShowerLength_titlename, 200, 0, 200);
    ShowerEnergyDiff_HistMap[fShowerModuleLabels[j]]      = new TH1F(ShowerEnergyDiff_name, ShowerEnergyDiff_titlename, 500, -5, 5);
    ShowerdEdx_HistMap[fShowerModuleLabels[j]]            = new TH1F(ShowerdEdX_name, ShowerdEdX_titlename, 500, 0, 20);
    EventSeggy_HistMap[fShowerModuleLabels[j]]            = new TH1F(EventSeggy_name, EventSeggy_titlename, 20, -10, 10);
    ShowerCompleteness_HistMap[fShowerModuleLabels[j]]    = new TH1F(ShowerCompleteness_name, ShowerCompleteness_titlename, 50, 0, 1);
    ShowerPurity_HistMap[fShowerModuleLabels[j]]          = new TH1F(ShowerPurity_name, ShowerPurity_titlename, 100, 0, 1);
    ShowerEnergy_HistMap[fShowerModuleLabels[j]]          = new TH1F(ShowerEnergy_name, ShowerEnergy_titlename, 20, 0, 5000);
    ShowerHitNum_HistMap[fShowerModuleLabels[j]]          = new TH1F(ShowerHitNum_name, ShowerHitNum_titlename, 12000, 0, 12000);
    HitEnergy_HistMap[fShowerModuleLabels[j]]             = new TH1F(HitEnergy_name, HitEnergy_titlename,500,0,5000);

    ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap[fShowerModuleLabels[j]]      = new TH1F(ShowerRecoEnergyVsTrueEnergyinRecoShower_name, ShowerRecoEnergyVsTrueEnergyinRecoShower_titlename, 100,-5,5);
    TrueEnergy_HistMap[fShowerModuleLabels[j]]      = new TH1F(TrueEnergy_name, TrueEnergy_titlename, 500,0,5000);

    ShowerMag_HistMap[fShowerModuleLabels[j]]             = new TH1F(ShowerMag_name, ShowerMag_titlename,20,0,200); 
    ShowerDirectionDiff_HistMap[fShowerModuleLabels[j]]   = new TH1F(ShowerDirectionDiff_name, ShowerDirectionDiff_titlename,100,-1,1);
 

    Energies_Mean_ShowerStart_X_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_Mean_ShowerStart_Y_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_Mean_ShowerStart_Z_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_Mean_ShowerLength_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_Mean_ShowerEnergyDiff_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_Mean_ShowerdEdx_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_Mean_EventSeggy_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_Mean_ShowerCompleteness_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_Mean_ShowerPurity_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_Mean_ShowerEnergy_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_Mean_ShowerHitNum_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_Mean_ShowerTotalEnergyDiff_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_Mean_HitEnergy_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_Mean_ShowerMag_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_Mean_ShowerDirectionDiff_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_Mean_TrueEnergy_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    
    Energies_RMS_ShowerStart_X_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_RMS_ShowerStart_Y_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_RMS_ShowerStart_Z_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_RMS_ShowerLength_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_RMS_ShowerEnergyDiff_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_RMS_ShowerdEdx_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_RMS_EventSeggy_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_RMS_ShowerCompleteness_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_RMS_ShowerPurity_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_RMS_ShowerEnergy_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_RMS_ShowerHitNum_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_RMS_ShowerTotalEnergyDiff_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_RMS_HitEnergy_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_RMS_ShowerMag_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_RMS_ShowerDirectionDiff_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
    Energies_RMS_TrueEnergy_GraphMap[fShowerModuleLabels[j]] = new TGraphErrors(0);
  
    for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){

      std::string ClusterProjectionMatchedEnergy_string = "ClusterProjectionMatchedEnergy " + fShowerModuleLabels[j] + "Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
      std::string ClusterCompletenessEnergy_string      = "ClusterCompletenessEnergy "      + fShowerModuleLabels[j] + "Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
      std::string ClusterPurityEnergy_string            = "ClusterPurityEnergy "            + fShowerModuleLabels[j] + "Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
      std::string ClusterCompletnessHits_string         = "ClusterCompletnessHits "         + fShowerModuleLabels[j] + "Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
      std::string ClusterPurityHits_string              = "ClusterPurityHits "              + fShowerModuleLabels[j] + "Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
      std::string ClusterCompPurityEnergy_string        = "ClusterCompPurityEnergy "        + fShowerModuleLabels[j] + "Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
      std::string ClusterCompPurityHits_string          = "ClusterCompPurityHits "          + fShowerModuleLabels[j] + "Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
      
      std::string ClusterProjectionMatchedEnergy_titlestring = ";Cluster Matched; Enteries";
      std::string ClusterCompletenessEnergy_titlestring      = ";Completeness; Enteries";
      std::string ClusterPurityEnergy_titlestring            = ";Purity; Enteries"; 
      std::string ClusterCompletnessHits_titlestring         = ";Completeness; Enteries";
      std::string ClusterPurityHits_titlestring              = ";Purity; Enteries";
      std::string ClusterCompPurityEnergy_titlestring        = ";Completeness * Purity; Enteries"; 
      std::string ClusterCompPurityHits_titlestring          = ";Completeness * Purity; Enteries";
   
      const char* ClusterProjectionMatchedEnergy_titlename = ClusterProjectionMatchedEnergy_titlestring.c_str();
      const char* ClusterCompletenessEnergy_titlename      = ClusterCompletenessEnergy_titlestring.c_str();
      const char* ClusterPurityEnergy_titlename            = ClusterPurityEnergy_titlestring.c_str();
      const char* ClusterCompletnessHits_titlename         =  ClusterCompletnessHits_titlestring.c_str();  
      const char* ClusterPurityHits_titlename              = ClusterPurityHits_titlestring.c_str();
      const char* ClusterCompPurityEnergy_titlename        = ClusterCompPurityEnergy_titlestring.c_str(); 
      const char* ClusterCompPurityHits_titlename          = ClusterCompPurityHits_titlestring.c_str();
  
      const char* ClusterProjectionMatchedEnergy_name = ClusterProjectionMatchedEnergy_string.c_str();
      const char* ClusterCompletenessEnergy_name      = ClusterCompletenessEnergy_string.c_str();
      const char* ClusterPurityEnergy_name            = ClusterPurityEnergy_string.c_str();
      const char* ClusterCompletnessHits_name         =  ClusterCompletnessHits_string.c_str();  
      const char* ClusterPurityHits_name              = ClusterPurityHits_string.c_str();
      const char* ClusterCompPurityEnergy_name        = ClusterCompPurityEnergy_string.c_str(); 
      const char* ClusterCompPurityHits_name          = ClusterCompPurityHits_string.c_str();
   
      ClusterProjectionMatchedEnergy_HistMap[fShowerModuleLabels[j]][plane_id] = new TH1F(ClusterProjectionMatchedEnergy_name,ClusterProjectionMatchedEnergy_titlename,4,-1,2);
      ClusterCompletenessEnergy_HistMap[fShowerModuleLabels[j]][plane_id]      = new TH1F(ClusterCompletenessEnergy_name,ClusterCompletenessEnergy_titlename,100,0,2);
      ClusterPurityEnergy_HistMap[fShowerModuleLabels[j]][plane_id]            = new TH1F(ClusterPurityEnergy_name,ClusterPurityEnergy_titlename,100,0,2);
      ClusterCompletnessHits_HistMap[fShowerModuleLabels[j]][plane_id]         = new TH1F(ClusterCompletnessHits_name,ClusterCompletnessHits_titlename,100,0,2);
      ClusterPurityHits_HistMap[fShowerModuleLabels[j]][plane_id]              = new TH1F(ClusterPurityHits_name,ClusterPurityHits_titlename,100,0,2);  
      ClusterCompPurityEnergy_HistMap[fShowerModuleLabels[j]][plane_id]        = new TH1F(ClusterCompPurityEnergy_name,ClusterCompPurityEnergy_titlename,100,0,2);
      ClusterCompPurityHits_HistMap[fShowerModuleLabels[j]][plane_id]          = new TH1F(ClusterCompPurityHits_name,ClusterCompPurityHits_titlename,100,0,2); 

      Energies_Mean_ClusterProjectionMatchedEnergy_GraphMap[fShowerModuleLabels[j]][plane_id] = new TGraphErrors(0);
      Energies_Mean_ClusterCompletenessEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]      = new TGraphErrors(0); 
      Energies_Mean_ClusterPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]            = new TGraphErrors(0); 
      Energies_Mean_ClusterCompletnessHits_GraphMap[fShowerModuleLabels[j]][plane_id]         = new TGraphErrors(0);
      Energies_Mean_ClusterPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]              = new TGraphErrors(0); 
      Energies_Mean_ClusterCompPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]        = new TGraphErrors(0);
      Energies_Mean_ClusterCompPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]          = new TGraphErrors(0); 

      Energies_RMS_ClusterProjectionMatchedEnergy_GraphMap[fShowerModuleLabels[j]][plane_id] = new TGraphErrors(0);
      Energies_RMS_ClusterCompletenessEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]      = new TGraphErrors(0); 
      Energies_RMS_ClusterPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]            = new TGraphErrors(0); 
      Energies_RMS_ClusterCompletnessHits_GraphMap[fShowerModuleLabels[j]][plane_id]         = new TGraphErrors(0);
      Energies_RMS_ClusterPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]              = new TGraphErrors(0); 
      Energies_RMS_ClusterCompPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]        = new TGraphErrors(0);
      Energies_RMS_ClusterCompPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]          = new TGraphErrors(0); 
      
    }
                 
    for(unsigned int i=0; i<fEnergies.size(); ++i){
    
      //Create the Histogram maps name and string for the different energies. 
      std::string Energies_ShowerStartX_string = "Shower Start X " + fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)"  + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)";  
      std::string Energies_ShowerStartY_string = "Shower Start Y " + fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
      std::string Energies_ShowerStartZ_string = "Shower Start Z " + fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
      std::string Energies_ShowerLength_string          = "Shower Length " + fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
      std::string Energies_ShowerEnergyDiff_string      =  "Shower Energy Diff " + fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
      std::string Energies_ShowerdEdX_string            = "Shower dEdX initial track" + fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
      std::string Energies_EventSeggy_string            = "Seggy " + fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
      std::string Energies_ShowerCompleteness_string    = "Shower Completeness " + fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
      std::string Energies_ShowerPurity_string          = "Shower Purity " + fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
      std::string Energies_ShowerEnergy_string          = "Shower Energy " + fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
      std::string Energies_HitEnergy_string             = "Shower Hit Energy " +  fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
      std::string Energies_ShowerHitNum_string          = "Shower Hit Num " + fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
      std::string Energies_EnergyRes_Mean_string        = "Energy Res Mean " + fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
      std::string Energies_EnergyRes_RMS_string         = "Energy Res RMS " +  fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
      std::string Energies_ShowerMag_string             = "Shower Mang Position " + fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
      std::string Energies_ShowerDirectionDiff_string         = "Shower Direction Difference" +  fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
      std::string Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_string = "Reco Energy/ True Energy in Reco Shower" + fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
      std::string Energies_TrueEnergy_string = "True Energy (MeV)" + fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
      
      std::string Title = "Simulated Energy " + std::to_string(fEnergies[i]) + " MeV";   
      
      std::string Energies_ShowerStartX_titlestring           = Title +  ";|True X Position-Reco X Postition(cm)|;Enteries";
      std::string Energies_ShowerStartY_titlestring           = Title +  ";|True Y Position-Reco Y Postition(cm)|;Enteries";
      std::string Energies_ShowerStartZ_titlestring           = Title +  ";|True Z Position-Reco Z Postition(cm)|;Enteries";
      std::string Energies_ShowerLength_titlestring           = Title +  ";|True Shower Length - Reco Shower Length| (cm); Enteries";
      std::string Energies_ShowerEnergyDiff_titlestring       = Title +  ";(Reco Energy)/TrueEnergy;Enteries";
      std::string Energies_ShowerdEdX_titlestring             = Title +  ";dEdx (MeV/cm);Enteries";
      std::string Energies_EventSeggy_titlestring             = Title +  ";Number of Showers;Enteries";
      std::string Energies_ShowerCompleteness_titlestring     = Title +  ";Completeness;Enteries";
      std::string Energies_ShowerPurity_titlestring           = Title +  ";Purity;Enteries;";
      std::string Energies_ShowerEnergy_titlestring           = Title +  ";Energy (MeV);Enteries";
      std::string Energies_HitEnergy_titlestring              = Title +  ";Energy (MeV);Enteries";
      std::string Energies_ShowerHitNum_titlestring           = Title +  "; Number of Hits;Enteries";
      std::string Energies_ShowerMag_titlestring              = Title +  ";|True Shower Start - RecoShower Start| (cm);Enteries";
      std::string Energies_ShowerDirectionDiff_titlestring     = Title +  ";(True Direction).(Reco Direction) (cm^2); Enteries"; 
      std::string Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_titlestring = Title +  "Reco Energy/ True Energy in Reco Shower; Enteries"; 
      std::string Energies_TrueEnergy_titlestring = Title +  "True Energy (MeV); Enteries"; 

      const char* Energies_ShowerStartX_titlename = Energies_ShowerStartX_titlestring.c_str();
      const char* Energies_ShowerStartY_titlename = Energies_ShowerStartY_titlestring.c_str();
      const char* Energies_ShowerStartZ_titlename = Energies_ShowerStartZ_titlestring.c_str();
      const char* Energies_ShowerLength_titlename = Energies_ShowerLength_titlestring.c_str();
      const char* Energies_ShowerEnergyDiff_titlename = Energies_ShowerEnergyDiff_titlestring.c_str();
      const char* Energies_ShowerdEdX_titlename =  Energies_ShowerdEdX_titlestring.c_str();
      const char* Energies_EventSeggy_titlename = Energies_EventSeggy_titlestring.c_str();
      const char* Energies_ShowerCompleteness_titlename = Energies_ShowerCompleteness_titlestring.c_str();
      const char* Energies_ShowerPurity_titlename = Energies_ShowerPurity_titlestring.c_str();
      const char* Energies_ShowerEnergy_titlename = Energies_ShowerEnergy_titlestring.c_str();
      const char* Energies_HitEnergy_titlename = Energies_HitEnergy_titlestring.c_str();
      const char* Energies_ShowerHitNum_titlename = Energies_ShowerHitNum_titlestring.c_str();
      const char* Energies_ShowerMag_titlename         = Energies_ShowerMag_titlestring.c_str();
      const char* Energies_ShowerDirectionDiff_titlename   = Energies_ShowerDirectionDiff_titlestring.c_str();
      const char* Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_titlename = Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_titlestring.c_str();
      const char* Energies_TrueEnergy_titlename = Energies_TrueEnergy_titlestring.c_str();
      
      const char* Energies_ShowerStartX_name          = Energies_ShowerStartX_string.c_str();
      const char* Energies_ShowerStartY_name          = Energies_ShowerStartY_string.c_str();
      const char* Energies_ShowerStartZ_name          = Energies_ShowerStartZ_string.c_str();
      const char* Energies_ShowerLength_name          = Energies_ShowerLength_string.c_str();
      const char* Energies_ShowerEnergyDiff_name      = Energies_ShowerEnergyDiff_string.c_str();
      const char* Energies_ShowerdEdX_name            = Energies_ShowerdEdX_string.c_str();
      const char* Energies_EventSeggy_name            = Energies_EventSeggy_string.c_str();
      const char* Energies_ShowerCompleteness_name    = Energies_ShowerCompleteness_string.c_str();
      const char* Energies_ShowerPurity_name          = Energies_ShowerPurity_string.c_str();
      const char* Energies_ShowerEnergy_name          = Energies_ShowerEnergy_string.c_str();
      const char* Energies_ShowerHitNum_name          = Energies_ShowerHitNum_string.c_str();
      const char* Energies_HitEnergy_name             = Energies_HitEnergy_string.c_str();
      const char* Energies_ShowerMag_name             = Energies_ShowerMag_string.c_str();
      const char* Energies_ShowerDirectionDiff_name   = Energies_ShowerDirectionDiff_string.c_str(); 
      const char* Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_name = Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_string.c_str();
      const char* Energies_TrueEnergy_name = Energies_TrueEnergy_string.c_str();
      
      //Create the histograms
    
      Energies_ShowerStart_X_HistMap[fShowerModuleLabels[j]][fEnergies[i]]  = new TH1F(Energies_ShowerStartX_name, Energies_ShowerStartX_titlename, 400, 0, 20);
      Energies_ShowerStart_Y_HistMap[fShowerModuleLabels[j]][fEnergies[i]]  = new TH1F(Energies_ShowerStartY_name, Energies_ShowerStartY_titlename, 400, 0, 20);
      Energies_ShowerStart_Z_HistMap[fShowerModuleLabels[j]][fEnergies[i]]  = new TH1F(Energies_ShowerStartZ_name, Energies_ShowerStartZ_titlename, 600, 0, 20);
      
      Energies_ShowerLength_HistMap[fShowerModuleLabels[j]][fEnergies[i]]          = new TH1F(Energies_ShowerLength_name, Energies_ShowerLength_titlename, 200, 0, 200);
      Energies_ShowerEnergyDiff_HistMap[fShowerModuleLabels[j]][fEnergies[i]]      = new TH1F(Energies_ShowerEnergyDiff_name, Energies_ShowerEnergyDiff_titlename, 500, -5, 5);
      Energies_ShowerdEdx_HistMap[fShowerModuleLabels[j]][fEnergies[i]]            = new TH1F(Energies_ShowerdEdX_name, Energies_ShowerdEdX_titlename, 500, 0, 20);
      Energies_EventSeggy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]            = new TH1F(Energies_EventSeggy_name, Energies_EventSeggy_titlename, 20, -10, 10);
      Energies_ShowerCompleteness_HistMap[fShowerModuleLabels[j]][fEnergies[i]]    = new TH1F(Energies_ShowerCompleteness_name, Energies_ShowerCompleteness_titlename, 50, 0, 1);
      Energies_ShowerPurity_HistMap[fShowerModuleLabels[j]][fEnergies[i]]          = new TH1F(Energies_ShowerPurity_name, Energies_ShowerPurity_titlename, 120, 0, 1.2);
      Energies_ShowerEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]          = new TH1F(Energies_ShowerEnergy_name, Energies_ShowerEnergy_titlename, 20, 0, 5000);
      Energies_ShowerHitNum_HistMap[fShowerModuleLabels[j]][fEnergies[i]]          = new TH1F(Energies_ShowerHitNum_name, Energies_ShowerHitNum_titlename, 12000, 0, 12000);
      Energies_HitEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]             = new TH1F(Energies_HitEnergy_name, Energies_HitEnergy_titlename,100,0,5000);
      
      Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap[fShowerModuleLabels[j]][fEnergies[i]]      = new TH1F(Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_name, Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_titlename, 100,-5,5);
      Energies_TrueEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]      = new TH1F(Energies_TrueEnergy_name, Energies_TrueEnergy_titlename, 500,0,5000);
      
      Energies_ShowerMag_HistMap[fShowerModuleLabels[j]][fEnergies[i]]             = new TH1F(Energies_ShowerMag_name, Energies_ShowerMag_titlename,50,0,200); 
      Energies_ShowerDirectionDiff_HistMap[fShowerModuleLabels[j]][fEnergies[i]]   = new TH1F(Energies_ShowerDirectionDiff_name, Energies_ShowerDirectionDiff_titlename,100,-1,1);


      for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){
	
	std::string Energies_ClusterProjectionMatchedEnergy_string = "ClusterProjectionMatchedEnergy " + fShowerModuleLabels[j] + "Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)   + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)";
	std::string Energies_ClusterCompletenessEnergy_string      = "ClusterCompletenessEnergy "      + fShowerModuleLabels[j] + "Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)   + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
	std::string Energies_ClusterPurityEnergy_string            = "ClusterPurityEnergy "            + fShowerModuleLabels[j] + "Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)   + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
	std::string Energies_ClusterCompletnessHits_string         = "ClusterCompletnessHits "         + fShowerModuleLabels[j] + "Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)   + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
	std::string Energies_ClusterPurityHits_string              = "ClusterPurityHits "              + fShowerModuleLabels[j] + "Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)   + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
	std::string Energies_ClusterCompPurityEnergy_string        = "ClusterCompPurityEnergy "        + fShowerModuleLabels[j] + "Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)   + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;
	std::string Energies_ClusterCompPurityHits_string          = "ClusterCompPurityHits "          + fShowerModuleLabels[j] + "Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)   + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)" ;

	std::string Energies_ClusterProjectionMatchedEnergy_titlestring = Title + "Cluster Matched; Enteries";
	std::string Energies_ClusterCompletenessEnergy_titlestring      = Title + "Completeness; Enteries";
	std::string Energies_ClusterPurityEnergy_titlestring            = Title + "Purity; Enteries"; 
	std::string Energies_ClusterCompletnessHits_titlestring         = Title + "Completeness; Enteries";
	std::string Energies_ClusterPurityHits_titlestring              = Title + "Purity; Enteries";
	std::string Energies_ClusterCompPurityEnergy_titlestring        = Title + "Completeness * Purity; Enteries"; 
	std::string Energies_ClusterCompPurityHits_titlestring          = Title + "Completeness * Purity; Enteries";
	
	const char* Energies_ClusterProjectionMatchedEnergy_titlename = Energies_ClusterProjectionMatchedEnergy_titlestring.c_str();
	const char* Energies_ClusterCompletenessEnergy_titlename      = Energies_ClusterCompletenessEnergy_titlestring.c_str();
	const char* Energies_ClusterPurityEnergy_titlename            = Energies_ClusterPurityEnergy_titlestring.c_str();
	const char* Energies_ClusterCompletnessHits_titlename         = Energies_ClusterCompletnessHits_titlestring.c_str();  
	const char* Energies_ClusterPurityHits_titlename              = Energies_ClusterPurityHits_titlestring.c_str();
	const char* Energies_ClusterCompPurityEnergy_titlename        = Energies_ClusterCompPurityEnergy_titlestring.c_str(); 
	const char* Energies_ClusterCompPurityHits_titlename          = Energies_ClusterCompPurityHits_titlestring.c_str();
	
	const char* Energies_ClusterProjectionMatchedEnergy_name = Energies_ClusterProjectionMatchedEnergy_string.c_str();
	const char* Energies_ClusterCompletenessEnergy_name      = Energies_ClusterCompletenessEnergy_string.c_str();
	const char* Energies_ClusterPurityEnergy_name            = Energies_ClusterPurityEnergy_string.c_str();
	const char* Energies_ClusterCompletnessHits_name         = Energies_ClusterCompletnessHits_string.c_str();  
	const char* Energies_ClusterPurityHits_name              = Energies_ClusterPurityHits_string.c_str();
	const char* Energies_ClusterCompPurityEnergy_name        = Energies_ClusterCompPurityEnergy_string.c_str(); 
	const char* Energies_ClusterCompPurityHits_name          = Energies_ClusterCompPurityHits_string.c_str();
	
	Energies_ClusterProjectionMatchedEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id] = new TH1F(Energies_ClusterProjectionMatchedEnergy_name,Energies_ClusterProjectionMatchedEnergy_titlename,4,-1,2);
	Energies_ClusterCompletenessEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]      = new TH1F(Energies_ClusterCompletenessEnergy_name,Energies_ClusterCompletenessEnergy_titlename,100,0,2);
	Energies_ClusterPurityEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]            = new TH1F(Energies_ClusterPurityEnergy_name,Energies_ClusterPurityEnergy_titlename,100,0,2);
	Energies_ClusterCompletnessHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]         = new TH1F(Energies_ClusterCompletnessHits_name,Energies_ClusterCompletnessHits_titlename,100,0,2);
	Energies_ClusterPurityHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]              = new TH1F(Energies_ClusterPurityHits_name,Energies_ClusterPurityHits_titlename,100,0,2);  
	Energies_ClusterCompPurityEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]        = new TH1F(Energies_ClusterCompPurityEnergy_name,Energies_ClusterCompPurityEnergy_titlename,100,0,2);
	Energies_ClusterCompPurityHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]          = new TH1F(Energies_ClusterCompPurityHits_name,Energies_ClusterCompPurityHits_titlename,100,0,2); 
	
      }
    }
  }//Shower Label Loop

  
 
  Energies_Mean_ShowerStart_X_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_Mean_ShowerStart_X_Multi",";Energy (MeV);Mean |True X Position-Reco X Postition| (cm)");
  Energies_Mean_ShowerStart_Y_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_Mean_ShowerStart_Y_Multi",";Energy (MeV);Mean |True Y Position-Reco Y Postition| (cm)");
  Energies_Mean_ShowerStart_Z_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_Mean_ShowerStart_Z_Multi",";Energy (MeV);Mean |True Z Position-Reco Z Postition| (cm)");
  Energies_Mean_ShowerLength_Multi  = tfs->makeAndRegister<TMultiGraph>("Energy_Mean_ShowerLength_Multi",";Energy (MeV);Mean |True Length - Reco Length |(cm)"); 
  Energies_Mean_ShowerEnergyDiff_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_Mean_ShowerEnergyDiff_Multi",";Energy (MeV); Mean Reco Energy/True Energy");
  Energies_Mean_ShowerdEdx_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_Mean_ShowerdEdx_Multi",";Energy (MeV); Mean dEdx (MeV/cm)");
  Energies_Mean_EventSeggy_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_Mean_EventSeggy_Multi",";Energy (MeV); Mean Number of Clusters");
  Energies_Mean_ShowerCompleteness_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_Mean_ShowerCompleteness_Multi",";Energy (MeV); Mean Completeness");
  Energies_Mean_ShowerPurity_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_Mean_ShowerPurity_Multi",";Energy (MeV); Mean Purity");
  Energies_Mean_ShowerEnergy_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_Mean_ShowerEnergy_Multi",";Energy (MeV); Mean Reco Energy (MeV)");
  Energies_Mean_ShowerHitNum_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_Mean_ShowerHitNum_Multi",";Energy (MeV); Mean Hit Number");
  Energies_Mean_HitEnergy_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_Mean_HitEnergy_Multi", ";Energy (MeV); Mean True Energy in Hits (MeV)");
  Energies_Mean_ShowerMag_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_Mean_ShowerMag_Multi", ";Energy (MeV); Mean |True Start - Reco Start| (cm)");
  Energies_Mean_ShowerDirectionDiff_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_Mean_ShowerDirectionDiff_Multi", ";Energy (Mev); Mean cosine of the angle between true and reco direction ");
  Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_Multi",";Energy (MeV); Mean Shower Reco Energy/ True Energy in Reco Shower");
  Energies_Mean_TrueEnergy_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_Mean_TrueEnergy_Multi",";Energy (MeV); Mean True Energy (MeV)");
  
  Energies_RMS_ShowerStart_X_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_RMS_ShowerStart_X_Multi",";Energy (MeV);RMS |True X Position-Reco X Postition| (cm)");
  Energies_RMS_ShowerStart_Y_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_RMS_ShowerStart_Y_Multi","1Energy (MeV);RMS |True Y Position-Reco Y Postition| (cm)");
  Energies_RMS_ShowerStart_Z_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_RMS_ShowerStart_Z_Multi",";Energy (MeV);RMS |True Z Position-Reco Z Postition| (cm)");
  Energies_RMS_ShowerLength_Multi  = tfs->makeAndRegister<TMultiGraph>("Energy_RMS_ShowerLength_Multi",";Energy (MeV);RMS |True Length - Reco Length |(cm)"); 
  Energies_RMS_ShowerEnergyDiff_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_RMS_ShowerEnergyDiff_Multi",";Energy (MeV); RMS Reco Energy/True Energy");
  Energies_RMS_ShowerdEdx_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_RMS_ShowerdEdx_Multi",";Energy (MeV); RMS dEdx (MeV/cm)");
  Energies_RMS_EventSeggy_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_RMS_EventSeggy_Multi",";Energy (MeV); RMS Number of Clusters");
  Energies_RMS_ShowerCompleteness_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_RMS_ShowerCompleteness_Multi",";Energy (MeV); RMS Completeness");
  Energies_RMS_ShowerPurity_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_RMS_ShowerPurity_Multi",";Energy (MeV); RMS Purity");
  Energies_RMS_ShowerEnergy_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_RMS_ShowerEnergy_Multi",";Energy (MeV); RMS Reco Energy (MeV)");
  Energies_RMS_ShowerHitNum_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_RMS_ShowerHitNum_Multi",";Energy (MeV); RMS Hit Number");
  Energies_RMS_HitEnergy_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_RMS_HitEnergy_Multi", ";Energy (MeV); RMS True Energy in Hits (MeV)");
  Energies_RMS_ShowerMag_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_RMS_ShowerMag_Multi", ";Energy (MeV); RMS |True Start - Reco Start| (cm)");
  Energies_RMS_ShowerDirectionDiff_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_RMS_ShowerDirectionDiff_Multi", ";Energy (Mev); RMS cosine of the angle between true and reco direction");
  Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_Multi",";Energy (MeV); RMS Shower Reco Energy/ True Energy in Reco Shower");
  Energies_RMS_TrueEnergy_Multi  = tfs->makeAndRegister<TMultiGraph>("Energies_RMS_TrueEnergy_Multi",";Energy (MeV); RMS True Energy");

  Energies_Mean_ShowerStart_X_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_Mean_ShowerStart_X_canvasMulti",";|True X Position-Reco X Postition| (cm);Enteries");
  Energies_Mean_ShowerStart_Y_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_Mean_ShowerStart_Y_canvasMulti",";|True Y Position-Reco Y Postition| (cm);Enteries");
  Energies_Mean_ShowerStart_Z_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_Mean_ShowerStart_Z_canvasMulti",";|True X Position-Reco X Postition|(cm);Enteries");
  Energies_Mean_ShowerLength_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_Mean_ShowerLength_canvasMulti",";|True Shower Length - Reco Shower Length| (cm); Enteries");
  Energies_Mean_ShowerEnergyDiff_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_Mean_ShowerEnergyDiff_canvasMulti",";(True Energy - Recon Energy)/TrueEnergy;Enteries");
  Energies_Mean_ShowerdEdx_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_Mean_ShowerdEdx_canvasMulti",";dEdx MeV/cm;Enteries");
  Energies_Mean_EventSeggy_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_Mean_EventSeggy_canvasMulti",";Number of Showers;Enteries");
  Energies_Mean_ShowerCompleteness_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_Mean_ShowerCompleteness_canvasMulti",";Completeness;Enteries");
  Energies_Mean_ShowerPurity_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_Mean_ShowerPurity_canvasMulti",";Purity;Enteries;");
  Energies_Mean_ShowerEnergy_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_Mean_ShowerEnergy_canvasMulti",";Energy (MeV);Enteries");
  Energies_Mean_ShowerHitNum_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_Mean_ShowerHitNum_canvasMulti",";Hits;Enteries");
  Energies_Mean_HitEnergy_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_Mean_HitEnergy_canvasMulti","Energy (MeV)");
  Energies_Mean_ShowerDirectionDiff_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_Mean_ShowerDirectionDiff_canvasMulti",";Angle Between Directions (rad);Enteries");
  Energies_Mean_ShowerMag_canvasMulti =tfs->makeAndRegister<TCanvas>("Energies_Mean_ShowerMag_canvasMulti",";|True Start - Reco Start| (cm)","Enteries");
  Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMulti= tfs->makeAndRegister<TCanvas>("Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMulti","");
  Energies_Mean_TrueEnergy_canvasMulti= tfs->makeAndRegister<TCanvas>("Energies_Mean_TrueEnergy_canvasMulti","");
 
  Energies_RMS_ShowerStart_X_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_RMS_ShowerStart_X_canvasMulti",";|True X Position-Reco X Postition| (cm);Enteries");
  Energies_RMS_ShowerStart_Y_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_RMS_ShowerStart_Y_canvasMulti",";|True Y Position-Reco Y Postition| (cm);Enteries");
  Energies_RMS_ShowerStart_Z_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_RMS_ShowerStart_Z_canvasMulti",";|True X Position-Reco X Postition|(cm);Enteries");
  Energies_RMS_ShowerLength_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_RMS_ShowerLength_canvasMulti",";|True Shower Length - Reco Shower Length| (cm); Enteries");
  Energies_RMS_ShowerEnergyDiff_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_RMS_ShowerEnergyDiff_canvasMulti",";(True Energy - Recon Energy)/TrueEnergy;Enteries");
  Energies_RMS_ShowerdEdx_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_RMS_ShowerdEdx_canvasMulti",";dEdx MeV/cm;Enteries");
  Energies_RMS_EventSeggy_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_RMS_EventSeggy_canvasMulti",";Number of Showers;Enteries");
  Energies_RMS_ShowerCompleteness_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_RMS_ShowerCompleteness_canvasMulti",";Completeness;Enteries");
  Energies_RMS_ShowerPurity_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_RMS_ShowerPurity_canvasMulti",";Purity;Enteries;");
  Energies_RMS_ShowerEnergy_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_RMS_ShowerEnergy_canvasMulti",";Energy (MeV);Enteries");
  Energies_RMS_ShowerHitNum_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_RMS_ShowerHitNum_canvasMulti",";Hits;Enteries");
  Energies_RMS_HitEnergy_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_RMS_HitEnergy_canvasMulti","Energy (MeV)");
  Energies_RMS_ShowerDirectionDiff_canvasMulti = tfs->makeAndRegister<TCanvas>("Energies_RMS_ShowerDirectionDiff_canvasMulti",";Angle Between Directions (rad);Enteries");
  Energies_RMS_ShowerMag_canvasMulti =tfs->makeAndRegister<TCanvas>("Energies_RMS_ShowerMag_canvasMulti",";|True Start - Reco Start| (cm)","Enteries");
  Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMulti= tfs->makeAndRegister<TCanvas>("Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMulti","");
  Energies_RMS_TrueEnergy_canvasMulti= tfs->makeAndRegister<TCanvas>("Energies_RMS_TrueEnergy_canvasMulti","");

  ShowerDirection_X_canvas = tfs->makeAndRegister<TCanvas>("ShowerDirection_X_canvas","Shower Direction X");
  ShowerDirection_Y_canvas = tfs->makeAndRegister<TCanvas>("ShowerDirection_Y_canvas","Shower Direction Y");
  ShowerDirection_Z_canvas = tfs->makeAndRegister<TCanvas>("ShowerDirection_Z_canvas","Shower Direction Z");
  ShowerStart_X_canvas = tfs->makeAndRegister<TCanvas>("ShowerStart_X_canvas",";|True X Position-Reco X Postition| (cm);Enteries");
  ShowerStart_Y_canvas = tfs->makeAndRegister<TCanvas>("ShowerStart_Y_canvas",";|True Y Position-Reco Y Postition| (cm);Enteries");
  ShowerStart_Z_canvas = tfs->makeAndRegister<TCanvas>("ShowerStart_Z_canvas",";|True X Position-Reco X Postition| (cm);Enteries");
  ShowerLength_canvas = tfs->makeAndRegister<TCanvas>("ShowerLength_canvas",";|True Shower Length - Reco Shower Length| (cm); Enteries");
  ShowerEnergyDiff_canvas = tfs->makeAndRegister<TCanvas>("ShowerEnergyDiff_canvas",";(True Energy - Recon Energy)/TrueEnergy;Enteries");
  ShowerdEdx_canvas = tfs->makeAndRegister<TCanvas>("ShowerdEdx_canvas",";dEdx MeV/cm;Enteries");
  EventSeggy_canvas = tfs->makeAndRegister<TCanvas>("EventSeggy_canvas",";Number of Showers;Enteries");
  ShowerCompleteness_canvas = tfs->makeAndRegister<TCanvas>("ShowerCompleteness_canvas",";Completeness;Enteries");
  ShowerPurity_canvas = tfs->makeAndRegister<TCanvas>("ShowerPurity_canvas",";Purity;Enteries;");
  ShowerEnergy_canvas = tfs->makeAndRegister<TCanvas>("ShowerEnergy_canvas",";Energy (MeV);Enteries");
  ShowerHitNum_canvas = tfs->makeAndRegister<TCanvas>("ShowerHitNum_canvas",";Hits;Enteries");
  HitEnergy_canvas = tfs->makeAndRegister<TCanvas>("HitEnergy_canvas","Energy (MeV)");
  ShowerDirectionDiff_canvas = tfs->makeAndRegister<TCanvas>("ShowerDirectionDiff_canvas",";Angle Between Directions (rad);Enteries");
  ShowerMag_canvas =tfs->makeAndRegister<TCanvas>("ShowerMag_canvas",";|True Start - Reco Start| (cm)","Enteries");
  ShowerRecoEnergyVsTrueEnergyinRecoShower_canvas= tfs->makeAndRegister<TCanvas>("ShowerRecoEnergyVsTrueEnergyinRecoShower","");
  TrueEnergy_canvas= tfs->makeAndRegister<TCanvas>("TrueEnergy","");

  for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){

    std::string TFileEnergies_ClusterProjectionMatchedEnergy_stringMean = "MeanClusterProjectionMatchedEnergy Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    std::string TFileEnergies_ClusterCompletenessEnergy_stringMean       = "MeanClusterCompletenessEnergy Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    std::string TFileEnergies_ClusterPurityEnergy_stringMean             = "MeanClusterPurityEnergy  Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    std::string TFileEnergies_ClusterCompletnessHits_stringMean          = "MeanClusterCompletnessHits Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    std::string TFileEnergies_ClusterPurityHits_stringMean               = "MeanClusterPurityHits Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    std::string TFileEnergies_ClusterCompPurityEnergy_stringMean         = "MeanClusterCompPurityEnergy Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    std::string TFileEnergies_ClusterCompPurityHits_stringMean           = "MeanClusterCompPurityHits Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;

    std::string TFileEnergies_ClusterProjectionMatchedEnergy_stringRMS = "RMSClusterProjectionMatchedEnergy Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    std::string TFileEnergies_ClusterCompletenessEnergy_stringRMS       = "RMSClusterCompletenessEnergy Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    std::string TFileEnergies_ClusterPurityEnergy_stringRMS             = "RMSClusterPurityEnergy Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    std::string TFileEnergies_ClusterCompletnessHits_stringRMS          = "RMSClusterCompletnessHits Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    std::string TFileEnergies_ClusterPurityHits_stringRMS               = "RMSClusterPurityHits Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    std::string TFileEnergies_ClusterCompPurityEnergy_stringRMS         = "RMSClusterCompPurityEnergy Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    std::string TFileEnergies_ClusterCompPurityHits_stringRMS           = "RMSClusterCompPurityHits Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    
    std::string TFile_ClusterProjectionMatchedEnergy_string = "ClusterProjectionMatchedEnergy Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    std::string TFile_ClusterCompletenessEnergy_string      = "ClusterCompletenessEnergy Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    std::string TFile_ClusterPurityEnergy_string            = "ClusterPurityEnergy Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    std::string TFile_ClusterCompletnessHits_string         = "ClusterCompletnessHits Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    std::string TFile_ClusterPurityHits_string              = "ClusterPurityHits Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    std::string TFile_ClusterCompPurityEnergy_string        = "ClusterCompPurityEnergy Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;
    std::string TFile_ClusterCompPurityHits_string          = "ClusterCompPurityHits Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat)  ;

    std::string TFileEnergies_ClusterProjectionMatchedEnergy_titlestring = ";Energy (MeV);Cluster Matched";
    std::string TFileEnergies_ClusterCompletenessEnergy_titlestring      = ";Energy (MeV);Completeness";
    std::string TFileEnergies_ClusterPurityEnergy_titlestring            = ";Energy (MeV);Purity"; 
    std::string TFileEnergies_ClusterCompletnessHits_titlestring         = ";Energy (MeV);Completeness";
    std::string TFileEnergies_ClusterPurityHits_titlestring              = ";Energy (MeV);Purity";
    std::string TFileEnergies_ClusterCompPurityEnergy_titlestring        = ";Energy (MeV);Completeness * Purity"; 
    std::string TFileEnergies_ClusterCompPurityHits_titlestring          = ";Energy (MeV);Completeness * Purity;";

    const char* TFileEnergies_ClusterProjectionMatchedEnergy_titlename = TFileEnergies_ClusterProjectionMatchedEnergy_titlestring.c_str();
    const char* TFileEnergies_ClusterCompletenessEnergy_titlename      = TFileEnergies_ClusterCompletenessEnergy_titlestring.c_str();
    const char* TFileEnergies_ClusterPurityEnergy_titlename            = TFileEnergies_ClusterPurityEnergy_titlestring.c_str();
    const char* TFileEnergies_ClusterCompletnessHits_titlename         = TFileEnergies_ClusterCompletnessHits_titlestring.c_str();  
    const char* TFileEnergies_ClusterPurityHits_titlename              = TFileEnergies_ClusterPurityHits_titlestring.c_str();
    const char* TFileEnergies_ClusterCompPurityEnergy_titlename        = TFileEnergies_ClusterCompPurityEnergy_titlestring.c_str(); 
    const char* TFileEnergies_ClusterCompPurityHits_titlename          = TFileEnergies_ClusterCompPurityHits_titlestring.c_str();

    const char* TFile_ClusterProjectionMatchedEnergy_name = TFile_ClusterProjectionMatchedEnergy_string.c_str();
    const char* TFile_ClusterCompletenessEnergy_name      = TFile_ClusterCompletenessEnergy_string.c_str();
    const char* TFile_ClusterPurityEnergy_name            = TFile_ClusterPurityEnergy_string.c_str();
    const char* TFile_ClusterCompletnessHits_name         = TFile_ClusterCompletnessHits_string.c_str();  
    const char* TFile_ClusterPurityHits_name              = TFile_ClusterPurityHits_string.c_str();
    const char* TFile_ClusterCompPurityEnergy_name        = TFile_ClusterCompPurityEnergy_string.c_str(); 
    const char* TFile_ClusterCompPurityHits_name          = TFile_ClusterCompPurityHits_string.c_str();

    const char* TFileEnergies_ClusterProjectionMatchedEnergy_nameMean = TFileEnergies_ClusterProjectionMatchedEnergy_stringMean.c_str();
    const char* TFileEnergies_ClusterCompletenessEnergy_nameMean      = TFileEnergies_ClusterCompletenessEnergy_stringMean.c_str();
    const char* TFileEnergies_ClusterPurityEnergy_nameMean            = TFileEnergies_ClusterPurityEnergy_stringMean.c_str();
    const char* TFileEnergies_ClusterCompletnessHits_nameMean         = TFileEnergies_ClusterCompletnessHits_stringMean.c_str();  
    const char* TFileEnergies_ClusterPurityHits_nameMean              = TFileEnergies_ClusterPurityHits_stringMean.c_str();
    const char* TFileEnergies_ClusterCompPurityEnergy_nameMean        = TFileEnergies_ClusterCompPurityEnergy_stringMean.c_str(); 
    const char* TFileEnergies_ClusterCompPurityHits_nameMean          = TFileEnergies_ClusterCompPurityHits_stringMean.c_str();

    const char* TFileEnergies_ClusterProjectionMatchedEnergy_nameRMS = TFileEnergies_ClusterProjectionMatchedEnergy_stringRMS.c_str();
    const char* TFileEnergies_ClusterCompletenessEnergy_nameRMS      = TFileEnergies_ClusterCompletenessEnergy_stringRMS.c_str();
    const char* TFileEnergies_ClusterPurityEnergy_nameRMS            = TFileEnergies_ClusterPurityEnergy_stringRMS.c_str();
    const char* TFileEnergies_ClusterCompletnessHits_nameRMS         = TFileEnergies_ClusterCompletnessHits_stringRMS.c_str();  
    const char* TFileEnergies_ClusterPurityHits_nameRMS              = TFileEnergies_ClusterPurityHits_stringRMS.c_str();
    const char* TFileEnergies_ClusterCompPurityEnergy_nameRMS        = TFileEnergies_ClusterCompPurityEnergy_stringRMS.c_str(); 
    const char* TFileEnergies_ClusterCompPurityHits_nameRMS          = TFileEnergies_ClusterCompPurityHits_stringRMS.c_str();
    
    
    Energies_Mean_ClusterProjectionMatchedEnergy_Multi[plane_id]  = tfs->makeAndRegister<TMultiGraph>(TFileEnergies_ClusterProjectionMatchedEnergy_nameMean,TFileEnergies_ClusterProjectionMatchedEnergy_titlename);
    Energies_Mean_ClusterCompletenessEnergy_Multi[plane_id]  = tfs->makeAndRegister<TMultiGraph>(TFileEnergies_ClusterCompletenessEnergy_nameMean,TFileEnergies_ClusterCompletenessEnergy_titlename);
    Energies_Mean_ClusterPurityEnergy_Multi[plane_id]  = tfs->makeAndRegister<TMultiGraph>(TFileEnergies_ClusterPurityEnergy_nameMean,TFileEnergies_ClusterPurityEnergy_titlename);
    Energies_Mean_ClusterCompletnessHits_Multi[plane_id]  = tfs->makeAndRegister<TMultiGraph>(TFileEnergies_ClusterCompletnessHits_nameMean,TFileEnergies_ClusterCompletnessHits_titlename);
    Energies_Mean_ClusterPurityHits_Multi[plane_id]  = tfs->makeAndRegister<TMultiGraph>(TFileEnergies_ClusterPurityHits_nameMean,TFileEnergies_ClusterPurityHits_titlename);
    Energies_Mean_ClusterCompPurityEnergy_Multi[plane_id]  = tfs->makeAndRegister<TMultiGraph>(TFileEnergies_ClusterCompPurityEnergy_nameMean,TFileEnergies_ClusterCompPurityEnergy_titlename);
    Energies_Mean_ClusterCompPurityHits_Multi[plane_id]  = tfs->makeAndRegister<TMultiGraph>(TFileEnergies_ClusterCompPurityHits_nameMean,TFileEnergies_ClusterCompPurityHits_titlename); 

    
    Energies_RMS_ClusterProjectionMatchedEnergy_Multi[plane_id]  = tfs->makeAndRegister<TMultiGraph>(TFileEnergies_ClusterProjectionMatchedEnergy_nameRMS,TFileEnergies_ClusterProjectionMatchedEnergy_titlename);
    Energies_RMS_ClusterCompletenessEnergy_Multi[plane_id]  = tfs->makeAndRegister<TMultiGraph>(TFileEnergies_ClusterCompletenessEnergy_nameRMS,TFileEnergies_ClusterCompletenessEnergy_titlename);
    Energies_RMS_ClusterPurityEnergy_Multi[plane_id]  = tfs->makeAndRegister<TMultiGraph>(TFileEnergies_ClusterPurityEnergy_nameRMS,TFileEnergies_ClusterPurityEnergy_titlename);
    Energies_RMS_ClusterCompletnessHits_Multi[plane_id]  = tfs->makeAndRegister<TMultiGraph>(TFileEnergies_ClusterCompletnessHits_nameRMS,TFileEnergies_ClusterCompletnessHits_titlename);
    Energies_RMS_ClusterPurityHits_Multi[plane_id]  = tfs->makeAndRegister<TMultiGraph>(TFileEnergies_ClusterPurityHits_nameRMS,TFileEnergies_ClusterPurityHits_titlename);
    Energies_RMS_ClusterCompPurityEnergy_Multi[plane_id]  = tfs->makeAndRegister<TMultiGraph>(TFileEnergies_ClusterCompPurityEnergy_nameRMS,TFileEnergies_ClusterCompPurityEnergy_titlename);
    Energies_RMS_ClusterCompPurityHits_Multi[plane_id]  = tfs->makeAndRegister<TMultiGraph>(TFileEnergies_ClusterCompPurityHits_nameRMS,TFileEnergies_ClusterCompPurityHits_titlename); 


    Energies_Mean_ClusterProjectionMatchedEnergy_canvasMulti[plane_id]  = tfs->makeAndRegister<TCanvas>(TFileEnergies_ClusterProjectionMatchedEnergy_nameMean,TFileEnergies_ClusterProjectionMatchedEnergy_titlename);
    Energies_Mean_ClusterCompletenessEnergy_canvasMulti[plane_id]  = tfs->makeAndRegister<TCanvas>(TFileEnergies_ClusterCompletenessEnergy_nameMean,TFileEnergies_ClusterCompletenessEnergy_titlename);
    Energies_Mean_ClusterPurityEnergy_canvasMulti[plane_id]  = tfs->makeAndRegister<TCanvas>(TFileEnergies_ClusterPurityEnergy_nameMean,TFileEnergies_ClusterPurityEnergy_titlename);
    Energies_Mean_ClusterCompletnessHits_canvasMulti[plane_id]  = tfs->makeAndRegister<TCanvas>(TFileEnergies_ClusterCompletnessHits_nameMean,TFileEnergies_ClusterCompletnessHits_titlename);
    Energies_Mean_ClusterPurityHits_canvasMulti[plane_id]  = tfs->makeAndRegister<TCanvas>(TFileEnergies_ClusterPurityHits_nameMean,TFileEnergies_ClusterPurityHits_titlename);
    Energies_Mean_ClusterCompPurityEnergy_canvasMulti[plane_id]  = tfs->makeAndRegister<TCanvas>(TFileEnergies_ClusterCompPurityEnergy_nameMean,TFileEnergies_ClusterCompPurityEnergy_titlename);
    Energies_Mean_ClusterCompPurityHits_canvasMulti[plane_id]  = tfs->makeAndRegister<TCanvas>(TFileEnergies_ClusterCompPurityHits_nameMean,TFileEnergies_ClusterCompPurityHits_titlename); 

    
    Energies_RMS_ClusterProjectionMatchedEnergy_canvasMulti[plane_id]  = tfs->makeAndRegister<TCanvas>(TFileEnergies_ClusterProjectionMatchedEnergy_nameRMS,TFileEnergies_ClusterProjectionMatchedEnergy_titlename);
    Energies_RMS_ClusterCompletenessEnergy_canvasMulti[plane_id]  = tfs->makeAndRegister<TCanvas>(TFileEnergies_ClusterCompletenessEnergy_nameRMS,TFileEnergies_ClusterCompletenessEnergy_titlename);
    Energies_RMS_ClusterPurityEnergy_canvasMulti[plane_id]  = tfs->makeAndRegister<TCanvas>(TFileEnergies_ClusterPurityEnergy_nameRMS,TFileEnergies_ClusterPurityEnergy_titlename);
    Energies_RMS_ClusterCompletnessHits_canvasMulti[plane_id]  = tfs->makeAndRegister<TCanvas>(TFileEnergies_ClusterCompletnessHits_nameRMS,TFileEnergies_ClusterCompletnessHits_titlename);
    Energies_RMS_ClusterPurityHits_canvasMulti[plane_id]  = tfs->makeAndRegister<TCanvas>(TFileEnergies_ClusterPurityHits_nameRMS,TFileEnergies_ClusterPurityHits_titlename);
    Energies_RMS_ClusterCompPurityEnergy_canvasMulti[plane_id]  = tfs->makeAndRegister<TCanvas>(TFileEnergies_ClusterCompPurityEnergy_nameRMS,TFileEnergies_ClusterCompPurityEnergy_titlename);
    Energies_RMS_ClusterCompPurityHits_canvasMulti[plane_id]  = tfs->makeAndRegister<TCanvas>(TFileEnergies_ClusterCompPurityHits_nameRMS,TFileEnergies_ClusterCompPurityHits_titlename); 

    ClusterProjectionMatchedEnergy_canvas[plane_id]= tfs->makeAndRegister<TCanvas>(TFile_ClusterProjectionMatchedEnergy_name,"Cluster Matched Correctly; Enteries");
    ClusterCompletenessEnergy_canvas[plane_id]= tfs->makeAndRegister<TCanvas>(TFile_ClusterCompletenessEnergy_name,";Completeness;Enteries");
    ClusterPurityEnergy_canvas[plane_id]= tfs->makeAndRegister<TCanvas>(TFile_ClusterPurityEnergy_name,";Purity;Enteries");
    ClusterCompletnessHits_canvas[plane_id]= tfs->makeAndRegister<TCanvas>(TFile_ClusterCompletnessHits_name,";Completeness;Enteries");
    ClusterPurityHits_canvas[plane_id]= tfs->makeAndRegister<TCanvas>(TFile_ClusterPurityHits_name,";Purity;Enteries");
    ClusterCompPurityEnergy_canvas[plane_id]= tfs->makeAndRegister<TCanvas>(TFile_ClusterCompPurityEnergy_name,";Purity * Completeness;Enteries");
    ClusterCompPurityHits_canvas[plane_id]= tfs->makeAndRegister<TCanvas>(TFile_ClusterCompPurityHits_name,";Purity * Completeness;Enteries");
  }
}


void ana::ShowerValidation::analyze(const art::Event& evt) {

  //Getting  MC truth information
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if(evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    {art::fill_ptr_vector(mclist, mctruthListHandle);}
  
  //Getting the SimWire Information
  //Get the SimChannels so that we can find the IDEs deposited on them.                                                                          
  art::Handle<std::vector<sim::SimChannel> > simChannelHandle;
  std::vector<art::Ptr<sim::SimChannel> > simchannels;
  if(evt.getByLabel(fLArGeantModuleLabel,simChannelHandle))
    {art::fill_ptr_vector(simchannels, simChannelHandle);}

  //Getting the Hit Information 
  art::Handle<std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hits;
  if(evt.getByLabel(fHitsModuleLabel,hitListHandle))
    {art::fill_ptr_vector(hits, hitListHandle);}
  std::cout << "hits.size(): " << hits.size() << std::endl;

  //I think that doing getManyByType kind of initalises the handles giving every particle product id. Doing this allows us to find handles for the individal hits later.

  //Get all the hits
  std::vector<art::Handle<std::vector<recob::Hit> > > hitHandles;
  evt.getManyByType(hitHandles);  

  //Get all the clusters
  std::vector<art::Handle<std::vector<recob::Cluster> > > clusterHandles;
  evt.getManyByType(clusterHandles);

  //Get all the pfparticles
  std::vector<art::Handle<std::vector<recob::PFParticle> > > pfpHandles;
  evt.getManyByType(pfpHandles);

  //Get the energy deposited in the hits.                                                                                                                                          
  float EnergyinHits = RecoUtils::TotalEnergyDepinHits(hits,2);

  //###############################################
  //### Get the Truth information for the event ### 
  //###############################################

  //List the particles in the event                                                                                                                                                
  const sim::ParticleList& particles = particleInventory->ParticleList();

  //Loop over the particles                                                                                                                                                        
  std::map<int,const simb::MCParticle*> trueParticles;
  std::map<int,const simb::MCParticle*> trueInitialParticles;
  std::map<int,float> trueParticleEnergy;
  std::map<int,std::vector<int> > ShowersMothers; //Mothers are the key Daughters are in the vector. 
  int num_of_showers_viaEcut = 0;
  int num_of_showers_viaDensitycut = 0;
  float simenergy=-99999999;

  std::map<int,std::map<geo::PlaneID,int> > MCWires_Track_Map = RecoUtils::NumberofMCWiresHitMap(simchannels);

  //Make a map of Track id and pdgcode
  for (sim::ParticleList::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt){
    const simb::MCParticle *particle = particleIt->second;
    trueParticleEnergy[particle->TrackId()] = 0;
    trueParticles[particle->TrackId()] = particle;
  }

  for (sim::ParticleList::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt){
    const simb::MCParticle *particle = particleIt->second;
    
    //Particles with mother 0 are the initial particles (neutrino events this is the particles generated after the interaction. Keep note of these. 
    if(particle->Mother() == 0){
      
      simenergy = particle->E(); 
      trueInitialParticles[particle->TrackId()] = particle;
      
      //Check to see if the particle is contained and find the trajectory of the particle                                                                                                                                            
      //Get the number of Traj points to loop over       
      unsigned int TrajPoints = particle->NumberTrajectoryPoints();
      
      //Get the startpoistion so we can get the initial tpc.
      const TLorentzVector StartPositionTrajP = particle->Position(0);
      double start_vtx[3] = {StartPositionTrajP.X() ,StartPositionTrajP.Y(), StartPositionTrajP.Z()};
      geo::TPCID init_idtpc = geom->FindTPCAtPosition(start_vtx);
      
      //Loop over the trajectory points (they are in order). Loop to find the start point. 
      for(unsigned int TrajPoints_it=0; TrajPoints_it<TrajPoints; ++TrajPoints_it){
	
	//Find the vertex of the vector                                                            
	const TLorentzVector PositionTrajP = particle->Position(TrajPoints_it);
	double vtx[3] = {PositionTrajP.X() ,PositionTrajP.Y(), PositionTrajP.Z()};
	
	//Find if the vertex is in the TPC. If so make it the start point.                 
	geo::TPCID idtpc = geom->FindTPCAtPosition(vtx);
	
	if(idtpc != init_idtpc ){std::cout <<"Particle outside the TPC" << std::endl; continue;}
      }
    }

    if(fVerbose > 1){std::cout << "True Particle with track ID: " << particle->TrackId() << " Has code of: " << particle->PdgCode() << " and Energy of: " << particle->E() << " With Mother: " << particle->Mother() << std::endl;}


    //Find The electron or photon mother of the shower.   
    const simb::MCParticle * particle_temp = particle;
    int ShowerMotherID = particle_temp->TrackId();



    if(particle_temp->Mother() != 0){
      if(trueParticles.find(ShowerMotherID) == trueParticles.end() || trueParticles.find(particle_temp->Mother()) == trueParticles.end()){continue;}

    
      while(particle_temp->Mother() != 0 && (TMath::Abs(trueParticles[particle_temp->Mother()]->PdgCode()) == 11 || trueParticles[particle_temp->Mother()]->PdgCode() == 22)){

	ShowerMotherID = particle_temp->Mother();

	particle_temp =  trueParticles[particle_temp->Mother()];

	if(trueParticles.find(particle_temp->Mother()) == trueParticles.end()){break;}
      }

    }

    if(ShowersMothers.find(ShowerMotherID) == ShowersMothers.end() && (TMath::Abs(trueParticles[ShowerMotherID]->PdgCode()) == 11 || trueParticles[ShowerMotherID]->PdgCode() == 22)){ShowersMothers[ShowerMotherID].push_back(ShowerMotherID);}
  }

  std::map<int,int> Daughters_used;

  //Find the all the daughter particles associated with the mothers
  for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end(); ++showermother){

    int Gen_Num = 0;
    int Daughter_id = showermother->first;
    while(Gen_Num!=0 || Daughters_used[showermother->first] != trueParticles[showermother->first]->NumberDaughters()){

      if(Daughters_used[Daughter_id] == trueParticles[Daughter_id]->NumberDaughters() || (trueParticles[Daughter_id]->PdgCode() != 22 && TMath::Abs(trueParticles[Daughter_id]->PdgCode()) != 11 )){
	if(trueParticles.find(Daughter_id) == trueParticles.end()){continue;}
	(showermother->second).push_back(Daughter_id);
	Daughter_id = trueParticles[Daughter_id]->Mother();
	--Gen_Num;
	++Daughters_used[Daughter_id];
	continue;
      }

      //Sometimes the duaghter is not in the particle list and so breaks the code. I presume its because its too small in energy to propergate by the simulation.      
      if(trueParticles.find(trueParticles[Daughter_id]->Daughter(Daughters_used[Daughter_id])) != trueParticles.end()){
      Daughter_id = trueParticles[Daughter_id]->Daughter(Daughters_used[Daughter_id]);
      ++Gen_Num;
      }
      else {
	++Daughters_used[Daughter_id];
      }
    }
  } 


  if(fVerbose > 1){
    for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end(); ++showermother){
      std::cout << " A Mother has track id: " << showermother->first << std::endl;
      for(std::vector<int>::iterator daughterID=(showermother->second).begin(); daughterID!=(showermother->second).end();++daughterID){
	std::cout << " Has a daughter of id: " << *daughterID << std::endl;
      }
    }
  }
   
  //Time to cut the true showers and make sure they are a shower.
  for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end();){
  
    //I've read that pair production starts to dominate at around ~100 MeV so to find how many showers we expect loop over the mother particle. Pi0=143.97 MeV min gammas = 71.985 MeV which is greater than that from electrons at ~100MeV so pi0 should always shower? So cut on anything below 100MeV in energy.
      
    //It ain't a shower I'm interested in if it didn't start with a pi0 or electron...probably.
    const simb::MCParticle *motherparticle = trueParticles[showermother->first];
    int pdgcode = motherparticle->PdgCode();
    if(pdgcode == 11 || pdgcode == 22){
      
      if(motherparticle->E() > 0.1){
	++num_of_showers_viaEcut;

      
	//using the RecoUtil function calculate the number of hits that see a charge deposition from the track.
	std::map<geo::PlaneID,int> Hit_num_map = RecoUtils::NumberofHitsThatContainEnergyDepositedByTracks(showermother->second, hits);
	
	//Calculaute the number of wires hit.  
	std::map<geo::PlaneID,int> Wire_num_map = ShowerUtils::NumberofWiresHitByShower(showermother->second, hits); 
   
	int low_density=0;

	//Compare hit density on the collection plane;
	for(std::map<geo::PlaneID,int>::iterator Hitnum_iter=Hit_num_map.begin(); Hitnum_iter!=Hit_num_map.end(); ++Hitnum_iter){
	  if(Wire_num_map[Hitnum_iter->first] == 0){continue;} 
	  double Hit_num = (Hitnum_iter->second);
	  double Wire_num = Wire_num_map[Hitnum_iter->first];
	  double Hit_Density = Hit_num/Wire_num;
	  if(Hit_Density > 1){
	    ++num_of_showers_viaDensitycut; 
	    ++low_density;
	    break;
	  }
	}
      
      	//If we don't have a density bigger than one in at least one plane then it aint a shower. This could be due to hit reco. 
	if(low_density == 0){
	  if(fVerbose > 0){std::cout << "Mother removed with id: " << showermother->first << " becuase the density is too low in the hit reconstruction" << std::endl;}
      	  showermother = ShowersMothers.erase(showermother);
	  continue;
	}
	++showermother;
      }
      else {
	if(fVerbose > 0){std::cout << "Mother removed with id: " << showermother->first << " becuase the true energy is too low" << std::endl;}
	showermother = ShowersMothers.erase(showermother);
      }
    }
    else {
      if(fVerbose > 0){std::cout << "Mother removed with id: " << showermother->first << " becuase it is a electron or photon" << std::endl;}
      showermother = ShowersMothers.erase(showermother);
    }
  }

  //Get the MC Energy deposited for each MC track. 
  std::map<int,float> MCTrack_Energy_map = RecoUtils::TrueEnergyDepositedFromMCTracks(simchannels);

  //Get the number of hits associated wit each Track. This is done for every hits handle in the event.
  std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > > MCTrack_hit_map;
  for(auto& handle : hitHandles){

    if(!handle.isValid()){
      mf::LogError("ShowerValidation") << "Bad hit handle from the all the hit handles" << std::endl;
      continue;
    }

    //RawHitFinder is a bit silly at the moment so ignore it for the time being - REMOVE 
    if(handle.provenance()->moduleLabel() == "fasthit"){continue;}
   
    //Getting the Hit Information 
    art::Handle<std::vector<recob::Hit> > hitHandle;
    std::vector<art::Ptr<recob::Hit> > hits_fromhandle;
    if(evt.getByLabel(handle.provenance()->moduleLabel(),hitHandle))
      {art::fill_ptr_vector(hits_fromhandle, hitHandle);}

    
    //Get a map of the number of hits per plane each track has deposited. 
    MCTrack_hit_map[handle.id()] = RecoUtils::NumberofPlaneHitsPerTrack(hits_fromhandle);
  }
  
  //##############################################
  //Get the reconstructed info and Match with MC #
  //##############################################
  
  for(unsigned int shwrlab_it=0; shwrlab_it<fShowerModuleLabels.size(); ++shwrlab_it){ 
    
    if(ShowersMothers.size() == 0){
      if(fVerbose > 0){std::cout << "No Mothers to Match to the showers" << std::endl;}
      continue;
    }

    //Getting the Shower Information
    art::Handle<std::vector<recob::Shower> > showerListHandle;
    std::vector<art::Ptr<recob::Shower> > showers;
    std::string fShowerModuleLabel = fShowerModuleLabels[shwrlab_it];
    if(evt.getByLabel(fShowerModuleLabel,showerListHandle))
      {art::fill_ptr_vector(showers,showerListHandle);}

    if(showers.size() == 0){
      if(fVerbose){std::cout << "No Shower in the Event" << std::endl;}
      continue;
    }

    //Association between Showers and 2d Hits
    art::FindManyP<recob::Hit>  fmh(showerListHandle, evt, fShowerModuleLabel);

    //Association between Showers and clusters   
    art::FindManyP<recob::Cluster> fmch(showerListHandle, evt, fShowerModuleLabel);
      
    //Association between Showers and pfParticle                                                                                            
    art::FindManyP<recob::PFParticle> fmpf(showerListHandle, evt, fShowerModuleLabel);

    HitEnergy_HistMap[fShowerModuleLabel]->Fill(EnergyinHits/1000*simenergy);

    std::vector< art::Ptr<recob::Hit> > showerhits; //hits in the shower    
    unsigned int max_hitnum=0;
    unsigned int biggest_shower_iter = 9999;

    for(unsigned int shower_iter = 0; shower_iter < showers.size(); ++shower_iter){
      //Get the shower  
      art::Ptr<recob::Shower>& shower = showers.at(shower_iter);
      
      //Get the hits vector from the shower
      showerhits = fmh.at(shower.key());
      
      if(showerhits.size() > max_hitnum){
	max_hitnum = showerhits.size();
	biggest_shower_iter = shower_iter;
      }
    }

    //Loop over hits associated with the shower add up the IDEs energy for each of the "track ID" and find the purity and compare other properites.
    
    //Loop over the showers in the event
    for(unsigned int shower_iter = 0; shower_iter < showers.size(); ++shower_iter){

      if(fUseBiggestShower == true){
	if(shower_iter != biggest_shower_iter){continue;}
      }

      //Get the shower
      art::Ptr<recob::Shower>& shower = showers.at(shower_iter);

        
      //#########################
      //### Shower Validation ###
      //#########################
      
      //Get the hits vector from the shower                                                          
      showerhits = fmh.at(shower.key());

      if(showerhits.size() == 0) {continue;}
      //      if(showerhits.size() < 200) {continue;}
      
      //Function from RecoUtils, finds the most probable track ID associated with the set of hits from there true energy depositons. The pair returns the energy as well. 
      std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain(ShowersMothers,showerhits,2); 
      int ShowerTrackID = ShowerTrackInfo.first;
      double TrueEnergyDepWithinShower_FromTrueShower = ShowerTrackInfo.second;

      double TrueEnergyDep_FromShower = 0;
      //Calculate the true Energy deposited By Shower 
      for(std::vector<int>::iterator daughterID=ShowersMothers[ShowerTrackID].begin(); daughterID!=ShowersMothers[ShowerTrackID].end(); ++daughterID){
	TrueEnergyDep_FromShower += MCTrack_Energy_map[*daughterID];
      }

      //Energy deposited within the set of Hits associated to the shower.
      double TrueEnergyDep_WithinRecoShower = 0; 
      
      //Loop over the hits and find the IDEs 
      for(std::vector< art::Ptr<recob::Hit> >::iterator hitIt=showerhits.begin(); hitIt!=showerhits.end(); ++hitIt){
	
	//Get the plane ID 
	geo::WireID wireid = (*hitIt)->WireID();
	int PlaneID = wireid.Plane;
	if(PlaneID != 2){continue;}
	
	//Split the Hit into its IDE for each track it associates with.                               
	std::vector<sim::TrackIDE> trackIDEs = backtracker->HitToTrackIDEs((*hitIt));
	for (unsigned int idIt = 0; idIt < trackIDEs.size(); ++idIt){ 
	  //Find the true total energy deposited in a set of hits. 
	  TrueEnergyDep_WithinRecoShower += trackIDEs.at(idIt).energy;
	}
      }//Hit Loop
      
      double completeness =  TrueEnergyDepWithinShower_FromTrueShower/TrueEnergyDep_FromShower;
      double purity       =  TrueEnergyDepWithinShower_FromTrueShower/TrueEnergyDep_WithinRecoShower;
      
      //Add to the respective hitograms.
      ShowerCompleteness_HistMap[fShowerModuleLabel]->Fill(completeness);
      ShowerPurity_HistMap[fShowerModuleLabel]->Fill(purity);
      
      //Find the MCParticle this shower associates to
      const simb::MCParticle* MCShowerParticle = trueParticles.at(ShowerTrackID);
      
      //Find the Energy of the particle: 
      //double Energy = MCShowerParticle->E();
      
      //Get the number of Traj points to loop over                                                    
      unsigned int TrajPoints = MCShowerParticle->NumberTrajectoryPoints();
      
      //Find the start and end points of the initial particle in order to compare the track length. 
      const TLorentzVector PositionTrajStart =  MCShowerParticle->Position(0);
      const TLorentzVector PositionTrajEnd   =  MCShowerParticle->Position(TrajPoints-1);
      
      //The three vecotor for track length is the shower direction 
      //TVector3  TrueShowerDirection = (PositionTrajEnd - PositionTrajStart).Vect();
      TVector3  TrueShowerDirection(MCShowerParticle->Px(), MCShowerParticle->Py(),MCShowerParticle->Pz());

      //Initial track lentgh of the shower.
      double TrueTrackLength = TrueShowerDirection.Mag();
      
      //Get the information for the shower  
      const int ShowerBest_Plane                       = shower->best_plane(); 
      const TVector3 & ShowerDirection                 = shower->Direction();
      const TVector3 & ShowerStart                     = shower->ShowerStart();//cm
      const double   & ShowerTrackLength               = shower->Length();//cm        
      const  std::vector< double > & ShowerEnergyPlanes = shower->Energy();//MeV
      const std::vector< double > & ShowerdEdX_vec     = shower->dEdx();//MeV/cm  

      std::vector<double> ShowerEnergyPlanes_remove(3);
      ShowerEnergyPlanes_remove[0] = ShowerEnergyPlanes[0];
      ShowerEnergyPlanes_remove[1] = ShowerEnergyPlanes[1];
      ShowerEnergyPlanes_remove[2] = (ShowerEnergyPlanes[2] - 0.00155171)*0.00155171/4.39964 + 4.39964;

      //Get the Errror in the position 
      //      double Position_Error = TMath::Sqrt(TMath::Power(TrueShowerDirection.X()-ShowerDirection.X(),2) + TMath::Power(TrueShowerDirection.Y()-ShowerDirection.Y(),2) + TMath::Power(TrueShowerDirection.Z()-ShowerDirection.Z(),2));

      double Start_diff =  TMath::Sqrt(TMath::Power(PositionTrajStart.X()-ShowerStart.X(),2) + TMath::Power(PositionTrajStart.Y()-ShowerStart.Y(),2) + TMath::Power(PositionTrajStart.Z()-ShowerStart.Z(),2));

      
      //Get the Fraction off the true value and fill the histograms.
      ShowerDirection_X_HistMap[fShowerModuleLabel]->Fill((TrueShowerDirection.Y()*ShowerDirection.Y() + TrueShowerDirection.Z()*ShowerDirection.Z())/(TMath::Sqrt((TrueShowerDirection.Y()*TrueShowerDirection.Y() + TrueShowerDirection.Z()*TrueShowerDirection.Z()))*TMath::Sqrt((ShowerDirection.Y()*ShowerDirection.Y() + ShowerDirection.Z()*ShowerDirection.Z()))));
      ShowerDirection_Y_HistMap[fShowerModuleLabel]->Fill((TrueShowerDirection.X()*ShowerDirection.X() + TrueShowerDirection.Z()*ShowerDirection.Z())/(TMath::Sqrt((TrueShowerDirection.X()*TrueShowerDirection.X() + TrueShowerDirection.Z()*TrueShowerDirection.Z()))*TMath::Sqrt((ShowerDirection.X()*ShowerDirection.X() + ShowerDirection.Z()*ShowerDirection.Z()))));
      ShowerDirection_Z_HistMap[fShowerModuleLabel]->Fill((TrueShowerDirection.Y()*ShowerDirection.Y() + TrueShowerDirection.X()*ShowerDirection.X())/(TMath::Sqrt((TrueShowerDirection.Y()*TrueShowerDirection.Y() + TrueShowerDirection.X()*TrueShowerDirection.X()))*TMath::Sqrt((ShowerDirection.Y()*ShowerDirection.Y() + ShowerDirection.X()*ShowerDirection.X()))));
      ShowerStart_X_HistMap[fShowerModuleLabel]->Fill(TMath::Abs(PositionTrajStart.X()-ShowerStart.X()));
      ShowerStart_Y_HistMap[fShowerModuleLabel]->Fill(TMath::Abs(PositionTrajStart.Y()-ShowerStart.Y()));
      ShowerStart_Z_HistMap[fShowerModuleLabel]->Fill(TMath::Abs(PositionTrajStart.Z()-ShowerStart.Z()));
      ShowerMag_HistMap[fShowerModuleLabel]->Fill(Start_diff);
      ShowerDirectionDiff_HistMap[fShowerModuleLabel]->Fill(TrueShowerDirection.Dot(ShowerDirection)/(TrueShowerDirection.Mag()*ShowerDirection.Mag()));
      ShowerLength_HistMap[fShowerModuleLabel]->Fill(TMath::Abs(TrueTrackLength-ShowerTrackLength));
      ShowerEnergyDiff_HistMap[fShowerModuleLabel]->Fill(ShowerEnergyPlanes_remove[ShowerBest_Plane]/TrueEnergyDep_FromShower);
      ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap[fShowerModuleLabel]->Fill(ShowerEnergyPlanes_remove[ShowerBest_Plane]/TrueEnergyDepWithinShower_FromTrueShower);
      ShowerEnergy_HistMap[fShowerModuleLabel]->Fill(ShowerEnergyPlanes_remove[ShowerBest_Plane]);
      ShowerHitNum_HistMap[fShowerModuleLabel]->Fill(showerhits.size());
      ShowerdEdx_HistMap[fShowerModuleLabel]->Fill((ShowerdEdX_vec[ShowerBest_Plane]));
   
      //Fill the Energy dependent Histograms. 
      if(fEnergies.size() != 0){
	for(unsigned int i=0; i<fEnergies.size(); ++i){
	  if(TMath::Abs(simenergy*1000 - fEnergies[i])< fEnergyWidth){

	    Energies_ShowerStart_X_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(TMath::Abs(PositionTrajStart.X()-ShowerStart.X()));
	    Energies_ShowerStart_Y_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(TMath::Abs(PositionTrajStart.Y()-ShowerStart.Y()));
	    Energies_ShowerStart_Z_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(TMath::Abs(PositionTrajStart.Z()-ShowerStart.Z()));
	    Energies_ShowerMag_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(Start_diff );
	    Energies_ShowerDirectionDiff_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(TrueShowerDirection.Dot(ShowerDirection)/(TrueShowerDirection.Mag()*ShowerDirection.Mag()));
	    Energies_ShowerEnergy_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(ShowerEnergyPlanes_remove[ShowerBest_Plane]);
	    Energies_ShowerLength_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(TMath::Abs(TrueTrackLength-ShowerTrackLength));
	    Energies_ShowerEnergyDiff_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(ShowerEnergyPlanes_remove[ShowerBest_Plane]/TrueEnergyDep_FromShower);
	    Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(ShowerEnergyPlanes_remove[ShowerBest_Plane]/3/TrueEnergyDepWithinShower_FromTrueShower);
	    Energies_ShowerHitNum_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(showerhits.size());
	    Energies_ShowerCompleteness_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(completeness);
	    Energies_ShowerPurity_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(purity);
	    Energies_ShowerdEdx_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill((ShowerdEdX_vec[ShowerBest_Plane]));
	  }
	}
      }

      if(fVerbose > 0){
	std::cout << "#################################################" << std::endl;
	std::cout << "Shower Label: " <<  fShowerModuleLabel << std::endl;
	std::cout << "True Start: " << PositionTrajStart.X() << " Shower Start: " << ShowerStart.X() << std::endl;
	std::cout << "X Poisition: " <<  ShowerStart.X() << "Y Position " << ShowerStart.Y() << " Z Poistion: " << ShowerDirection.Z() << std::endl;
	std::cout << "ShowerBest_Plane: " << ShowerBest_Plane << std::endl;
	std::cout << "TrueTrackLength: " << TrueTrackLength << " ShowerTrackLength: " << ShowerTrackLength << std::endl;
	std::cout << "Number of True Showers That pass the E cut: " <<  num_of_showers_viaEcut << std::endl;
	std::cout << "Number of True Showers That pass the Density cut: " << num_of_showers_viaDensitycut << std::endl;
	std::cout << "Best Plane Shower energy " << ShowerEnergyPlanes_remove[ShowerBest_Plane] << std::endl;
	std::cout << "TrueEnergyDepWithinShower_FromTrueShower : " <<  TrueEnergyDepWithinShower_FromTrueShower  << std::endl;
	std::cout << "TrueEnergyDep_FromShower: " << TrueEnergyDep_FromShower << std::endl;
	std::cout << "TrueEnergy deposited by hits in shower: " <<  TrueEnergyDep_WithinRecoShower << std::endl;
	std::cout << "Purity: " << purity << " completeness: " << completeness << std::endl;
	std::cout << "Hit Size: " << showerhits.size() << std::endl;
	std::cout << "Energy in Hits: " << EnergyinHits << std::endl;
	std::cout << "Energy Simulated: " << simenergy << std::endl;
	std::cout << "#################################################" <<std::endl;
      }

      //##########################
      //### Cluster Validation ###
      //##########################

      //Get the clusters associated to the shower.
      if(fmch.isValid()){
	art::Handle<std::vector<recob::Cluster > > clusterHandle;
	evt.get(fmch.at(shower.key()).front().id(),clusterHandle);
	if(clusterHandle.isValid()){
	  std::vector<art::Ptr<recob::Cluster> > showerclusters = fmch.at(shower.key());
	  ana::ShowerValidation::ClusterValidation(showerclusters,evt,clusterHandle,ShowersMothers,MCTrack_Energy_map,MCTrack_hit_map,ShowerTrackID,simenergy,fShowerModuleLabel);
	}
	else{
	  mf::LogError("ShowerValidation") << "Cluster handle is stale. No clustering validation done" << std::endl;
	}
      }
      else if(fmpf.isValid()){

	//Find the Clusters associated to PF particle.
	art::Handle<std::vector<recob::PFParticle> > pfpHandle;  
	evt.get(fmpf.at(shower.key()).front().id(),pfpHandle);
	if(pfpHandle.isValid()){
	  art::FindManyP<recob::Cluster> fmcpf(pfpHandle, evt, pfpHandle.provenance()->moduleLabel());
	  if(fmcpf.isValid()){ 
	    art::Handle<std::vector<recob::Cluster > > clusterHandle;
	    evt.get(fmcpf.at(0).front().id(),clusterHandle);
	    if(clusterHandle.isValid()){
	      std::vector< art::Ptr<recob::Cluster> > showerclusters = fmcpf.at(shower.key());
	      ana::ShowerValidation::ClusterValidation(showerclusters,evt,clusterHandle,ShowersMothers,MCTrack_Energy_map,MCTrack_hit_map,ShowerTrackID,simenergy,fShowerModuleLabel);
	    }
	    else{
	      mf::LogError("ShowerValidation") << "Cluster handle is stale. No clustering validation done" << std::endl;
	    }
	  }
	  else{
	    mf::LogError("ShowerValidation") << "No Assosoication between pf particles and clusters was found for shower made from a pf particle. No clustering validation was done." << std::endl;
	  }
	}
	else{
	  mf::LogError("ShowerValidation") << "pf particle handle is stale" << std::endl;
	}
      }
      else{
	mf::LogError("ShowerValidation") << "No cluster or pandora association" << std::endl;
      }
    


    }//Shower Loop 
    
    //Whats the segementyness of the event.
    EventSeggy_HistMap[fShowerModuleLabel]->Fill(showers.size()/num_of_showers_viaDensitycut);
    TrueEnergy_HistMap[fShowerModuleLabel]->Fill(simenergy*1000);

    if(fEnergies.size() != 0){
      for(unsigned int i=0; i<fEnergies.size(); ++i){
	if(TMath::Abs(simenergy*1000 - fEnergies[i])< fEnergyWidth){
	  Energies_HitEnergy_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(EnergyinHits/simenergy*1000);
	  Energies_EventSeggy_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(showers.size()/num_of_showers_viaDensitycut);
	  Energies_TrueEnergy_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(simenergy*1000);
	}
      }
    }
  }
  return; 
}

void ana::ShowerValidation::ClusterValidation(std::vector< art::Ptr<recob::Cluster> >& clusters, const art::Event& evt, art::Handle<std::vector<recob::Cluster> >& clusterHandle, std::map<int,std::vector<int> >& ShowerMotherTrackIDs, std::map<int,float>& MCTrack_Energy_map, std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > >& MCTrack_hit_map, int& TrueShowerID, float& simenergy, std::string& fShowerModuleLabel){

  //Get the associated hits 
  art::FindManyP<recob::Hit> fmhc(clusterHandle, evt, clusterHandle.provenance()->moduleLabel());

  //Holder for cluster hits
  std::vector< art::Ptr<recob::Hit> > clusterhits; 

  //Get the Hits Handle used for this cluster type 
  art::Handle<std::vector<recob::Hit > > hitHandle;
  evt.get(fmhc.at(0).front().id(),hitHandle);

  if(!hitHandle.isValid()){   
    mf::LogError("ShowerValidation") << "Hits handle is stale. No clustering validation done" << std::endl;
    return; 
  }

  //Get the hits vector from the shower
  for(auto const& cluster : clusters){

    clusterhits = fmhc.at(cluster.key());
    
    //Function from RecoUtils, finds the most probable track ID associated with the set of hits from there true energy depositons. The pair returns the energy in the hits as well.
    std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain(ShowerMotherTrackIDs,clusterhits, cluster->Plane().Plane);
    float TotalTrueEnergy = 0;
    float signalhits      = 0;
    float totalhits       = 0;
    for(std::vector<int>::iterator daughterID=ShowerMotherTrackIDs[ShowerTrackInfo.first].begin(); daughterID!=ShowerMotherTrackIDs[ShowerTrackInfo.first].end(); ++daughterID){
      
      //Calculate the true Energy deposited By Shower  
      TotalTrueEnergy += MCTrack_Energy_map[*daughterID];

      //Count how many hits are from the true shower.
      signalhits += RecoUtils::NumberofHitsFromTrack(*daughterID, clusterhits);

      //Count how many hits are missed in the plane from the true track id.
      totalhits += MCTrack_hit_map[hitHandle.id()][*daughterID][cluster->Plane()];
    }

    int projection_match = -999;

    //Have we matched the 2D cluster to the correct shower correctly. In terms of Energy depositions: 
    if(ShowerTrackInfo.first == TrueShowerID){projection_match = 1;}
    else{projection_match = 0;}

    //Calculate the purity and completeness metrics
    float completeness_hits  = signalhits/totalhits;
    float purity_hits        = signalhits/clusterhits.size();

    float completness_energy = ShowerTrackInfo.second/TotalTrueEnergy; 
    float purity_energy      = ShowerTrackInfo.second/RecoUtils::TotalEnergyDepinHits(clusterhits,cluster->Plane().Plane);

    //Add to the histograms    
    ClusterProjectionMatchedEnergy_HistMap[fShowerModuleLabel][cluster->Plane()]->Fill(projection_match);
    ClusterCompletenessEnergy_HistMap[fShowerModuleLabel][cluster->Plane()]     ->Fill(completness_energy);
    ClusterPurityEnergy_HistMap[fShowerModuleLabel][cluster->Plane()]           ->Fill(purity_energy);
    ClusterCompletnessHits_HistMap[fShowerModuleLabel][cluster->Plane()]        ->Fill(completeness_hits);
    ClusterPurityHits_HistMap[fShowerModuleLabel][cluster->Plane()]             ->Fill(purity_hits);
    ClusterCompPurityEnergy_HistMap[fShowerModuleLabel][cluster->Plane()]       ->Fill(completness_energy*purity_energy);
    ClusterCompPurityHits_HistMap[fShowerModuleLabel][cluster->Plane()]         ->Fill(completeness_hits*purity_hits);
    
    //Add to the Energy dependent histograms
    if(fEnergies.size() != 0){
      for(unsigned int i=0; i<fEnergies.size(); ++i){
        if(TMath::Abs(simenergy*1000 - fEnergies[i])< fEnergyWidth){
	  Energies_ClusterProjectionMatchedEnergy_HistMap[fShowerModuleLabel][fEnergies[i]][cluster->Plane()]->Fill(projection_match);
	  Energies_ClusterCompletenessEnergy_HistMap[fShowerModuleLabel][fEnergies[i]][cluster->Plane()]     ->Fill(completness_energy);
	  Energies_ClusterPurityEnergy_HistMap[fShowerModuleLabel][fEnergies[i]][cluster->Plane()]           ->Fill(purity_energy);
	  Energies_ClusterCompletnessHits_HistMap[fShowerModuleLabel][fEnergies[i]][cluster->Plane()]        ->Fill(completeness_hits);
	  Energies_ClusterPurityHits_HistMap[fShowerModuleLabel][fEnergies[i]][cluster->Plane()]             ->Fill(purity_hits);
	  Energies_ClusterCompPurityEnergy_HistMap[fShowerModuleLabel][fEnergies[i]][cluster->Plane()]       ->Fill(completness_energy*purity_energy);
	  Energies_ClusterCompPurityHits_HistMap[fShowerModuleLabel][fEnergies[i]][cluster->Plane()]         ->Fill(completeness_hits*purity_hits);
	}
      }
    }

  }//Cluster Loop
    return;
}

void ana::ShowerValidation::endJob() {
  
  std::vector<int> colours; 
    for(unsigned int j=0; j<fShowerModuleLabels.size(); ++j){
      colours.push_back(j+2);
    }
 
  //Create the legend                                   
  TLegend *leg = new TLegend(0.65, 0.1, 0.9, 0.35);
  //leg->SetFillColor(0);
  //leg->SetHeader("Shower Modules");
  //leg->SetBorderSize(1);
  //leg->SetTextSize(0.025);
  
  int fillstyle = 3000;

  for(unsigned int j=0; j<fShowerModuleLabels.size(); ++j){

    ++fillstyle; 

    const char* name_legend = fShowerModuleLabels[j].c_str();
    leg->AddEntry(ShowerDirection_X_HistMap[fShowerModuleLabels[j]], name_legend);

    ShowerDirection_X_canvas->cd();
    ShowerDirection_X_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    ShowerDirection_X_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    ShowerDirection_X_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    ShowerDirection_X_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    ShowerDirection_X_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();

    ShowerDirection_Y_canvas->cd();
    ShowerDirection_Y_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    ShowerDirection_Y_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    ShowerDirection_Y_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    ShowerDirection_Y_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    ShowerDirection_Y_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();

    ShowerDirection_Z_canvas->cd();
    ShowerDirection_Z_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    ShowerDirection_Z_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    ShowerDirection_Z_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    ShowerDirection_Z_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    ShowerDirection_Z_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();

    ShowerStart_X_canvas->cd();
    ShowerStart_X_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    ShowerStart_X_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    ShowerStart_X_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    ShowerStart_X_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    ShowerStart_X_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();    
  
    ShowerStart_Y_canvas->cd();
    ShowerStart_Y_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    ShowerStart_Y_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    ShowerStart_Y_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    ShowerStart_Y_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    ShowerStart_Y_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();

    ShowerStart_Z_canvas->cd();
    ShowerStart_Z_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    ShowerStart_Z_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    ShowerStart_Z_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    ShowerStart_Z_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    ShowerStart_Z_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();

    ShowerLength_canvas->cd();
    ShowerLength_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    ShowerLength_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    ShowerLength_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    ShowerLength_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    ShowerLength_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();

    ShowerEnergyDiff_canvas->cd();
    ShowerEnergyDiff_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    ShowerEnergyDiff_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    ShowerEnergyDiff_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    ShowerEnergyDiff_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    ShowerEnergyDiff_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();

    ShowerdEdx_canvas->cd();
    ShowerdEdx_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    ShowerdEdx_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    ShowerdEdx_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    ShowerdEdx_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    ShowerdEdx_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();

    EventSeggy_canvas->cd();
    EventSeggy_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    EventSeggy_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    EventSeggy_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    EventSeggy_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    EventSeggy_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();

    ShowerCompleteness_canvas->cd();
    ShowerCompleteness_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    ShowerCompleteness_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    ShowerCompleteness_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    ShowerCompleteness_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    ShowerCompleteness_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();

    ShowerPurity_canvas->cd();
    ShowerPurity_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    ShowerPurity_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    ShowerPurity_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    ShowerPurity_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    ShowerPurity_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();

    ShowerEnergy_canvas->cd();
    ShowerEnergy_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    ShowerEnergy_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    ShowerEnergy_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    ShowerEnergy_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    ShowerEnergy_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();

    ShowerHitNum_canvas->cd();
    ShowerHitNum_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    ShowerHitNum_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    ShowerHitNum_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    ShowerHitNum_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    ShowerHitNum_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();

    HitEnergy_canvas->cd();
    HitEnergy_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    HitEnergy_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    HitEnergy_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    HitEnergy_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    HitEnergy_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();

    ShowerMag_canvas->cd();
    ShowerMag_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    ShowerMag_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    ShowerMag_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    ShowerMag_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    ShowerMag_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();
  

    ShowerDirectionDiff_canvas->cd();
    ShowerDirectionDiff_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    ShowerDirectionDiff_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    ShowerDirectionDiff_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    ShowerDirectionDiff_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    ShowerDirectionDiff_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();

    ShowerRecoEnergyVsTrueEnergyinRecoShower_canvas->cd();
    ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();

    TrueEnergy_canvas->cd();
    TrueEnergy_HistMap[fShowerModuleLabels[j]]->SetDirectory(0);
    TrueEnergy_HistMap[fShowerModuleLabels[j]]->SetFillColor(colours[j]);
    TrueEnergy_HistMap[fShowerModuleLabels[j]]->SetFillStyle(fillstyle);
    TrueEnergy_HistMap[fShowerModuleLabels[j]]->SetLineColor(1);
    TrueEnergy_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    leg->Draw();

    Energies_Mean_ShowerStart_X_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_Mean_ShowerStart_Y_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_Mean_ShowerStart_Z_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_Mean_ShowerLength_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_Mean_ShowerEnergyDiff_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_Mean_ShowerdEdx_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_Mean_EventSeggy_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_Mean_ShowerCompleteness_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_Mean_ShowerPurity_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_Mean_ShowerEnergy_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_Mean_ShowerHitNum_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_Mean_ShowerTotalEnergyDiff_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_Mean_HitEnergy_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_Mean_ShowerMag_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_Mean_ShowerDirectionDiff_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_Mean_TrueEnergy_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);

    Energies_RMS_ShowerStart_X_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_RMS_ShowerStart_Y_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_RMS_ShowerStart_Z_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_RMS_ShowerLength_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_RMS_ShowerEnergyDiff_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_RMS_ShowerdEdx_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_RMS_EventSeggy_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_RMS_ShowerCompleteness_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_RMS_ShowerPurity_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_RMS_ShowerEnergy_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_RMS_ShowerHitNum_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_RMS_ShowerTotalEnergyDiff_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_RMS_HitEnergy_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_RMS_ShowerMag_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_RMS_ShowerDirectionDiff_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_RMS_TrueEnergy_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);

    Energies_Mean_ShowerStart_X_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_Mean_ShowerStart_Y_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_Mean_ShowerStart_Z_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_Mean_ShowerLength_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_Mean_ShowerEnergyDiff_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_Mean_ShowerdEdx_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_Mean_EventSeggy_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_Mean_ShowerCompleteness_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_Mean_ShowerPurity_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_Mean_ShowerEnergy_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_Mean_ShowerHitNum_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_Mean_ShowerTotalEnergyDiff_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_Mean_HitEnergy_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_Mean_ShowerMag_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_Mean_ShowerDirectionDiff_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_Mean_TrueEnergy_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    
    Energies_RMS_ShowerStart_X_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_RMS_ShowerStart_Y_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_RMS_ShowerStart_Z_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_RMS_ShowerLength_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_RMS_ShowerEnergyDiff_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_RMS_ShowerdEdx_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_RMS_EventSeggy_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_RMS_ShowerCompleteness_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_RMS_ShowerPurity_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_RMS_ShowerEnergy_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_RMS_ShowerHitNum_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_RMS_ShowerTotalEnergyDiff_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_RMS_HitEnergy_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_RMS_ShowerMag_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_RMS_ShowerDirectionDiff_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);
    Energies_RMS_TrueEnergy_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);

    for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){

    	ClusterProjectionMatchedEnergy_canvas[plane_id]->cd();
    	ClusterProjectionMatchedEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->SetDirectory(0);
    	ClusterProjectionMatchedEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->SetFillColor(colours[j]);
    	ClusterProjectionMatchedEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->SetFillStyle(fillstyle);
    	ClusterProjectionMatchedEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->SetLineColor(1);
    	ClusterProjectionMatchedEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->Draw("SAME");
    	leg->Draw();

    	ClusterCompletenessEnergy_canvas[plane_id]->cd();
    	ClusterCompletenessEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->SetDirectory(0);
    	ClusterCompletenessEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->SetFillColor(colours[j]);
    	ClusterCompletenessEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->SetFillStyle(fillstyle);
    	ClusterCompletenessEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->SetLineColor(1);
    	ClusterCompletenessEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->Draw("SAME");
    	leg->Draw();

    	ClusterPurityEnergy_canvas[plane_id]->cd();
    	ClusterPurityEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->SetDirectory(0);
    	ClusterPurityEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->SetFillColor(colours[j]);
    	ClusterPurityEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->SetFillStyle(fillstyle);
    	ClusterPurityEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->SetLineColor(1);
    	ClusterPurityEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->Draw("SAME");
    	leg->Draw();

    	ClusterCompletnessHits_canvas[plane_id]->cd();
    	ClusterCompletnessHits_HistMap[fShowerModuleLabels[j]][plane_id]->SetDirectory(0);
    	ClusterCompletnessHits_HistMap[fShowerModuleLabels[j]][plane_id]->SetFillColor(colours[j]);
    	ClusterCompletnessHits_HistMap[fShowerModuleLabels[j]][plane_id]->SetFillStyle(fillstyle);
    	ClusterCompletnessHits_HistMap[fShowerModuleLabels[j]][plane_id]->SetLineColor(1);
    	ClusterCompletnessHits_HistMap[fShowerModuleLabels[j]][plane_id]->Draw("SAME");
    	leg->Draw();

    	ClusterPurityHits_canvas[plane_id]->cd();
    	ClusterPurityHits_HistMap[fShowerModuleLabels[j]][plane_id]->SetDirectory(0);
    	ClusterPurityHits_HistMap[fShowerModuleLabels[j]][plane_id]->SetFillColor(colours[j]);
    	ClusterPurityHits_HistMap[fShowerModuleLabels[j]][plane_id]->SetFillStyle(fillstyle);
    	ClusterPurityHits_HistMap[fShowerModuleLabels[j]][plane_id]->SetLineColor(1);
    	ClusterPurityHits_HistMap[fShowerModuleLabels[j]][plane_id]->Draw("SAME");
    	leg->Draw();

    	ClusterCompPurityEnergy_canvas[plane_id]->cd();
    	ClusterCompPurityEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->SetDirectory(0);
    	ClusterCompPurityEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->SetFillColor(colours[j]);
    	ClusterCompPurityEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->SetFillStyle(fillstyle);
    	ClusterCompPurityEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->SetLineColor(1);
    	ClusterCompPurityEnergy_HistMap[fShowerModuleLabels[j]][plane_id]->Draw("SAME");
    	leg->Draw();

    	ClusterCompPurityHits_canvas[plane_id]->cd();
    	ClusterCompPurityHits_HistMap[fShowerModuleLabels[j]][plane_id]->SetDirectory(0);
    	ClusterCompPurityHits_HistMap[fShowerModuleLabels[j]][plane_id]->SetFillColor(colours[j]);
    	ClusterCompPurityHits_HistMap[fShowerModuleLabels[j]][plane_id]->SetFillStyle(fillstyle);
    	ClusterCompPurityHits_HistMap[fShowerModuleLabels[j]][plane_id]->SetLineColor(1);
    	ClusterCompPurityHits_HistMap[fShowerModuleLabels[j]][plane_id]->Draw("SAME");
    	leg->Draw();

    	Energies_Mean_ClusterProjectionMatchedEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerColor(colours[j]);
    	Energies_Mean_ClusterCompletenessEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerColor(colours[j]);
    	Energies_Mean_ClusterPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerColor(colours[j]);
    	Energies_Mean_ClusterCompletnessHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerColor(colours[j]);
    	Energies_Mean_ClusterPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerColor(colours[j]);
    	Energies_Mean_ClusterCompPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerColor(colours[j]);
    	Energies_Mean_ClusterCompPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerColor(colours[j]);
    	Energies_Mean_ClusterProjectionMatchedEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerStyle(8);
    	Energies_Mean_ClusterCompletenessEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerStyle(8);
    	Energies_Mean_ClusterPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerStyle(8);
    	Energies_Mean_ClusterCompletnessHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerStyle(8);
    	Energies_Mean_ClusterPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerStyle(8);
    	Energies_Mean_ClusterCompPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerStyle(8);
    	Energies_Mean_ClusterCompPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerStyle(8);

    	Energies_RMS_ClusterProjectionMatchedEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerColor(colours[j]);
    	Energies_RMS_ClusterCompletenessEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerColor(colours[j]);
    	Energies_RMS_ClusterPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerColor(colours[j]);
    	Energies_RMS_ClusterCompletnessHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerColor(colours[j]);
    	Energies_RMS_ClusterPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerColor(colours[j]);
    	Energies_RMS_ClusterCompPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerColor(colours[j]);
    	Energies_RMS_ClusterCompPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerColor(colours[j]);
    	Energies_RMS_ClusterProjectionMatchedEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerStyle(8);
    	Energies_RMS_ClusterCompletenessEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerStyle(8);
    	Energies_RMS_ClusterPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerStyle(8);
    	Energies_RMS_ClusterCompletnessHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerStyle(8);
    	Energies_RMS_ClusterPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerStyle(8);
    	Energies_RMS_ClusterCompPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerStyle(8);
    	Energies_RMS_ClusterCompPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerStyle(8);
	
    }

    
    for(unsigned int i=0; i<fEnergies.size(); ++i){

      float Enteries = Energies_TrueEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetEntries();

      if(Enteries == 0){continue;}
      
      Energies_TrueEnergy_canvasMap[fEnergies[i]]->cd();
      Energies_TrueEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetDirectory(0);
      Energies_TrueEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillColor(colours[j]);
      Energies_TrueEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillStyle(fillstyle);
      Energies_TrueEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetLineColor(1);
      Energies_TrueEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->Draw("SAME");
      leg->Draw();
      float Energies_TrueEnergy_Mean = Energies_TrueEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();
      float Energies_TrueEnergy_RMS = Energies_TrueEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
      float Error_on_True_Energy = Energies_TrueEnergy_RMS/TMath::Sqrt(Enteries);
      Energies_Mean_TrueEnergy_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_Mean_TrueEnergy_GraphMap[fShowerModuleLabels[j]]->GetN(),fEnergies[i],Energies_TrueEnergy_Mean);
      Energies_RMS_TrueEnergy_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_RMS_TrueEnergy_GraphMap[fShowerModuleLabels[j]]->GetN(),fEnergies[i],Energies_TrueEnergy_RMS);

      Energies_ShowerStart_X_canvasMap[fEnergies[i]]->cd();
      Energies_ShowerStart_X_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetDirectory(0);
      Energies_ShowerStart_X_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillColor(colours[j]);
      Energies_ShowerStart_X_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillStyle(fillstyle);
      Energies_ShowerStart_X_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetLineColor(1);
      Energies_ShowerStart_X_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->Draw("SAME");
      leg->Draw();
      float Energies_ShowerStart_X_Mean = Energies_ShowerStart_X_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();
      float Energies_ShowerStart_X_RMS = Energies_ShowerStart_X_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
      Energies_Mean_ShowerStart_X_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_Mean_ShowerStart_X_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerStart_X_Mean);
      Energies_RMS_ShowerStart_X_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_RMS_ShowerStart_X_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerStart_X_RMS);
      Energies_Mean_ShowerStart_X_GraphMap[fShowerModuleLabels[j]]->SetPointError(Energies_Mean_ShowerStart_X_GraphMap[fShowerModuleLabels[j]]->GetN()-1,Error_on_True_Energy,Energies_ShowerStart_X_RMS/TMath::Sqrt(Enteries));
          
      Energies_ShowerStart_Y_canvasMap[fEnergies[i]]->cd();
      Energies_ShowerStart_Y_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetDirectory(0);
      Energies_ShowerStart_Y_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillColor(colours[j]);
      Energies_ShowerStart_Y_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillStyle(fillstyle);
      Energies_ShowerStart_Y_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetLineColor(1);
      Energies_ShowerStart_Y_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->Draw("SAME");
      leg->Draw();
      float Energies_ShowerStart_Y_Mean = Energies_ShowerStart_Y_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();
      float Energies_ShowerStart_Y_RMS = Energies_ShowerStart_Y_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
      Energies_Mean_ShowerStart_Y_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_Mean_ShowerStart_Y_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerStart_Y_Mean);
      Energies_RMS_ShowerStart_Y_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_RMS_ShowerStart_Y_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerStart_Y_RMS);
      Energies_Mean_ShowerStart_Y_GraphMap[fShowerModuleLabels[j]]->SetPointError(Energies_Mean_ShowerStart_Y_GraphMap[fShowerModuleLabels[j]]->GetN()-1,Error_on_True_Energy,Energies_ShowerStart_Y_RMS/TMath::Sqrt(Enteries));

      Energies_ShowerStart_Z_canvasMap[fEnergies[i]]->cd();
      Energies_ShowerStart_Z_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetDirectory(0);
      Energies_ShowerStart_Z_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillColor(colours[j]);
      Energies_ShowerStart_Z_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillStyle(fillstyle);
      Energies_ShowerStart_Z_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetLineColor(1);
      Energies_ShowerStart_Z_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->Draw("SAME");
      leg->Draw();
      float Energies_ShowerStart_Z_Mean = Energies_ShowerStart_Z_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();
      float Energies_ShowerStart_Z_RMS = Energies_ShowerStart_Z_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
      Energies_Mean_ShowerStart_Z_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_Mean_ShowerStart_Z_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerStart_Z_Mean);
      Energies_RMS_ShowerStart_Z_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_RMS_ShowerStart_Z_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerStart_Z_RMS);
      Energies_Mean_ShowerStart_Z_GraphMap[fShowerModuleLabels[j]]->SetPointError(Energies_Mean_ShowerStart_Z_GraphMap[fShowerModuleLabels[j]]->GetN()-1,Error_on_True_Energy,Energies_ShowerStart_Z_RMS/TMath::Sqrt(Enteries));

      

      Energies_ShowerLength_canvasMap[fEnergies[i]]->cd();
      Energies_ShowerLength_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetDirectory(0);
      Energies_ShowerLength_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillColor(colours[j]);
      Energies_ShowerLength_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillStyle(fillstyle);
      Energies_ShowerLength_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetLineColor(1);
      Energies_ShowerLength_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->Draw("SAME");
      leg->Draw();
      float Energies_ShowerLength_Mean = Energies_ShowerLength_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();
      float Energies_ShowerLength_RMS = Energies_ShowerLength_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
      Energies_Mean_ShowerLength_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_Mean_ShowerStart_Z_GraphMap[fShowerModuleLabels[j]]->GetN(), fEnergies[i],Energies_ShowerLength_Mean);
      Energies_RMS_ShowerLength_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_RMS_ShowerStart_Z_GraphMap[fShowerModuleLabels[j]]->GetN(), fEnergies[i],Energies_ShowerLength_RMS);
      Energies_Mean_ShowerLength_GraphMap[fShowerModuleLabels[j]]->SetPointError(Energies_Mean_ShowerLength_GraphMap[fShowerModuleLabels[j]]->GetN()-1,Error_on_True_Energy,Energies_ShowerLength_RMS/TMath::Sqrt(Enteries));

    
      Energies_ShowerEnergyDiff_canvasMap[fEnergies[i]]->cd();
      Energies_ShowerEnergyDiff_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetDirectory(0);
      Energies_ShowerEnergyDiff_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillColor(colours[j]);
      Energies_ShowerEnergyDiff_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillStyle(fillstyle);
      Energies_ShowerEnergyDiff_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetLineColor(1);
      Energies_ShowerEnergyDiff_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->Draw("SAME");
      leg->Draw();
      float Energies_ShowerEnergyDiff_Mean = Energies_ShowerEnergyDiff_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();
      float Energies_ShowerEnergyDiff_RMS = Energies_ShowerEnergyDiff_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
      Energies_Mean_ShowerEnergyDiff_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_Mean_ShowerEnergyDiff_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerEnergyDiff_Mean);
      Energies_RMS_ShowerEnergyDiff_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_RMS_ShowerEnergyDiff_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerEnergyDiff_RMS);
      Energies_Mean_ShowerEnergyDiff_GraphMap[fShowerModuleLabels[j]]->SetPointError(Energies_Mean_ShowerEnergyDiff_GraphMap[fShowerModuleLabels[j]]->GetN()-1,Error_on_True_Energy,Energies_ShowerEnergyDiff_RMS/TMath::Sqrt(Enteries));

      Energies_ShowerdEdx_canvasMap[fEnergies[i]]->cd();
      Energies_ShowerdEdx_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetDirectory(0);
      Energies_ShowerdEdx_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillColor(colours[j]);
      Energies_ShowerdEdx_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillStyle(fillstyle);
      Energies_ShowerdEdx_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetLineColor(1);
      Energies_ShowerdEdx_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->Draw("SAME");
      leg->Draw();
      float Energies_ShowerdEdx_Mean = Energies_ShowerdEdx_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();
      float Energies_ShowerdEdx_RMS = Energies_ShowerdEdx_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
      Energies_Mean_ShowerdEdx_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_Mean_ShowerdEdx_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerdEdx_Mean);
      Energies_RMS_ShowerdEdx_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_RMS_ShowerdEdx_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerdEdx_RMS);
      Energies_Mean_ShowerdEdx_GraphMap[fShowerModuleLabels[j]]->SetPointError(Energies_Mean_ShowerdEdx_GraphMap[fShowerModuleLabels[j]]->GetN()-1,Error_on_True_Energy,Energies_ShowerdEdx_RMS/TMath::Sqrt(Enteries));
	
      Energies_EventSeggy_canvasMap[fEnergies[i]]->cd();
      Energies_EventSeggy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetDirectory(0);
      Energies_EventSeggy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillColor(colours[j]);
      Energies_EventSeggy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillStyle(fillstyle);
      Energies_EventSeggy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetLineColor(1);
      Energies_EventSeggy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->Draw("SAME");
      leg->Draw();
      float Energies_EventSeggy_Mean = Energies_EventSeggy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();
      float Energies_EventSeggy_RMS = Energies_EventSeggy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
      Energies_Mean_EventSeggy_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_Mean_EventSeggy_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_EventSeggy_Mean);
      Energies_RMS_EventSeggy_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_RMS_EventSeggy_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_EventSeggy_RMS);
      Energies_Mean_EventSeggy_GraphMap[fShowerModuleLabels[j]]->SetPointError(Energies_Mean_EventSeggy_GraphMap[fShowerModuleLabels[j]]->GetN()-1,Error_on_True_Energy,Energies_EventSeggy_RMS/TMath::Sqrt(Enteries));

      Energies_ShowerCompleteness_canvasMap[fEnergies[i]]->cd();
      Energies_ShowerCompleteness_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetDirectory(0);
      Energies_ShowerCompleteness_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillColor(colours[j]);
      Energies_ShowerCompleteness_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillStyle(fillstyle);
      Energies_ShowerCompleteness_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetLineColor(1);
      Energies_ShowerCompleteness_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->Draw("SAME");
      leg->Draw();
      float Energies_ShowerCompleteness_Mean =  Energies_ShowerCompleteness_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();
      float Energies_ShowerCompleteness_RMS = Energies_ShowerCompleteness_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
      Energies_Mean_ShowerCompleteness_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_Mean_ShowerCompleteness_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerCompleteness_Mean);
      Energies_RMS_ShowerCompleteness_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_RMS_ShowerCompleteness_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerCompleteness_RMS);
      Energies_Mean_ShowerCompleteness_GraphMap[fShowerModuleLabels[j]]->SetPointError(Energies_Mean_ShowerCompleteness_GraphMap[fShowerModuleLabels[j]]->GetN()-1,Error_on_True_Energy,Energies_ShowerCompleteness_RMS/TMath::Sqrt(Enteries));


      Energies_ShowerPurity_canvasMap[fEnergies[i]]->cd();
      Energies_ShowerPurity_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetDirectory(0);
      Energies_ShowerPurity_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillColor(colours[j]);
      Energies_ShowerPurity_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillStyle(fillstyle);
      Energies_ShowerPurity_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetLineColor(1);
      Energies_ShowerPurity_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->Draw("SAME");
      leg->Draw();
      float Energies_ShowerPurity_Mean = Energies_ShowerPurity_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();
      float Energies_ShowerPurity_RMS = Energies_ShowerPurity_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
      Energies_Mean_ShowerPurity_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_Mean_ShowerPurity_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerPurity_Mean);
      Energies_RMS_ShowerPurity_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_RMS_ShowerPurity_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerPurity_RMS);
      Energies_Mean_ShowerPurity_GraphMap[fShowerModuleLabels[j]]->SetPointError(Energies_Mean_ShowerPurity_GraphMap[fShowerModuleLabels[j]]->GetN()-1,Error_on_True_Energy,Energies_ShowerPurity_RMS/TMath::Sqrt(Enteries));


      Energies_ShowerEnergy_canvasMap[fEnergies[i]]->cd();
      Energies_ShowerEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetDirectory(0);
      Energies_ShowerEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillColor(colours[j]);
      Energies_ShowerEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillStyle(fillstyle);
      Energies_ShowerEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetLineColor(1);
      Energies_ShowerEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->Draw("SAME");
      leg->Draw();
      float Energies_ShowerEnergy_Mean = Energies_ShowerEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();
      float Energies_ShowerEnergy_RMS = Energies_ShowerEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
      Energies_Mean_ShowerEnergy_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_Mean_ShowerEnergy_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerEnergy_Mean);
      Energies_RMS_ShowerEnergy_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_RMS_ShowerEnergy_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerEnergy_RMS);
      Energies_Mean_ShowerEnergy_GraphMap[fShowerModuleLabels[j]]->SetPointError(Energies_Mean_ShowerEnergy_GraphMap[fShowerModuleLabels[j]]->GetN()-1,Error_on_True_Energy,Energies_ShowerEnergy_RMS/TMath::Sqrt(Enteries));


      Energies_ShowerHitNum_canvasMap[fEnergies[i]]->cd();
      Energies_ShowerHitNum_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetDirectory(0);
      Energies_ShowerHitNum_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillColor(colours[j]);
      Energies_ShowerHitNum_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillStyle(fillstyle);
      Energies_ShowerHitNum_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetLineColor(1);
      Energies_ShowerHitNum_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->Draw("SAME");
      leg->Draw();
      float Energies_ShowerHitNum_Mean =  Energies_ShowerHitNum_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();
      float Energies_ShowerHitNum_RMS =  Energies_ShowerHitNum_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
      Energies_Mean_ShowerHitNum_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_Mean_ShowerHitNum_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerHitNum_Mean);
      Energies_RMS_ShowerHitNum_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_RMS_ShowerHitNum_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerHitNum_RMS);

      Energies_HitEnergy_canvasMap[fEnergies[i]]->cd();
      Energies_HitEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetDirectory(0);
      Energies_HitEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillColor(colours[j]);
      Energies_HitEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillStyle(fillstyle);
      Energies_HitEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetLineColor(1);
      Energies_HitEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->Draw("SAME");
      leg->Draw();
      float Energies_HitEnergy_Mean =  Energies_HitEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();
      float Energies_HitEnergy_RMS =  Energies_HitEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
      Energies_Mean_HitEnergy_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_Mean_HitEnergy_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_HitEnergy_Mean);
      Energies_RMS_HitEnergy_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_RMS_HitEnergy_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_HitEnergy_RMS);
      Energies_Mean_HitEnergy_GraphMap[fShowerModuleLabels[j]]->SetPointError(Energies_Mean_HitEnergy_GraphMap[fShowerModuleLabels[j]]->GetN()-1,Error_on_True_Energy,Energies_HitEnergy_RMS/TMath::Sqrt(Enteries));


      Energies_ShowerMag_canvasMap[fEnergies[i]]->cd();
      Energies_ShowerMag_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetDirectory(0);
      Energies_ShowerMag_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillColor(colours[j]);
      Energies_ShowerMag_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillStyle(fillstyle);
      Energies_ShowerMag_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetLineColor(1);
      Energies_ShowerMag_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->Draw("SAME");
      leg->Draw();
      float Energies_ShowerMag_Mean =  Energies_ShowerMag_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();
      float Energies_ShowerMag_RMS =  Energies_ShowerMag_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
      Energies_Mean_ShowerMag_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_Mean_ShowerMag_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerMag_Mean);
      Energies_RMS_ShowerMag_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_RMS_ShowerMag_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerMag_RMS);
      Energies_Mean_ShowerMag_GraphMap[fShowerModuleLabels[j]]->SetPointError(Energies_Mean_ShowerMag_GraphMap[fShowerModuleLabels[j]]->GetN()-1,Error_on_True_Energy,Energies_ShowerMag_RMS/TMath::Sqrt(Enteries));

      Energies_ShowerDirectionDiff_canvasMap[fEnergies[i]]->cd();
      Energies_ShowerDirectionDiff_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetDirectory(0);
      Energies_ShowerDirectionDiff_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillColor(colours[j]);
      Energies_ShowerDirectionDiff_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillStyle(fillstyle);
      Energies_ShowerDirectionDiff_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetLineColor(1);
      Energies_ShowerDirectionDiff_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->Draw("SAME");
      leg->Draw();
      float Energies_ShowerDirectionDiff_Mean =  Energies_ShowerDirectionDiff_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();
      float Energies_ShowerDirectionDiff_RMS =  Energies_ShowerDirectionDiff_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
      Energies_Mean_ShowerDirectionDiff_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_Mean_ShowerDirectionDiff_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerDirectionDiff_Mean);
      Energies_RMS_ShowerDirectionDiff_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_RMS_ShowerDirectionDiff_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerDirectionDiff_RMS);
      Energies_Mean_ShowerDirectionDiff_GraphMap[fShowerModuleLabels[j]]->SetPointError(Energies_Mean_ShowerDirectionDiff_GraphMap[fShowerModuleLabels[j]]->GetN()-1,Error_on_True_Energy,Energies_ShowerDirectionDiff_RMS/TMath::Sqrt(Enteries));
       
      Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMap[fEnergies[i]]->cd();
      Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetDirectory(0);
      Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillColor(colours[j]);
      Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetFillStyle(fillstyle);
      Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetLineColor(1);
      Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->Draw("SAME");
      leg->Draw();
      float Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_Mean =  Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();
      float Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_RMS =  Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
      Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_Mean);
      Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_RMS);
      Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap[fShowerModuleLabels[j]]->SetPointError(Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap[fShowerModuleLabels[j]]->GetN()-1,Error_on_True_Energy,Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_RMS/TMath::Sqrt(Enteries));
     
      for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){
	
      	Energies_ClusterProjectionMatchedEnergy_canvasMap[fEnergies[i]][plane_id]->cd();
      	Energies_ClusterProjectionMatchedEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetDirectory(0);
      	Energies_ClusterProjectionMatchedEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetFillColor(colours[j]);
      	Energies_ClusterProjectionMatchedEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetFillStyle(fillstyle);
      	Energies_ClusterProjectionMatchedEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetLineColor(1);
      	Energies_ClusterProjectionMatchedEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->Draw("SAME");
      	leg->Draw();
      	float Energies_ClusterProjectionMatchedEnergy_Mean =  Energies_ClusterProjectionMatchedEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->GetMean();
      	float Energies_ClusterProjectionMatchedEnergy_RMS =  Energies_ClusterProjectionMatchedEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->GetRMS();
      	Energies_Mean_ClusterProjectionMatchedEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPoint(Energies_Mean_ClusterProjectionMatchedEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN(),Energies_TrueEnergy_Mean,Energies_ClusterProjectionMatchedEnergy_Mean);
      	Energies_RMS_ClusterProjectionMatchedEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPoint(Energies_RMS_ClusterProjectionMatchedEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN(),Energies_TrueEnergy_Mean,Energies_ClusterProjectionMatchedEnergy_RMS);
      	Energies_Mean_ClusterProjectionMatchedEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPointError(Energies_Mean_ClusterProjectionMatchedEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN()-1,Error_on_True_Energy,Energies_ClusterProjectionMatchedEnergy_RMS/TMath::Sqrt(Enteries));

      	Energies_ClusterCompletenessEnergy_canvasMap[fEnergies[i]][plane_id]->cd();
      	Energies_ClusterCompletenessEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetDirectory(0);
      	Energies_ClusterCompletenessEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetFillColor(colours[j]);
      	Energies_ClusterCompletenessEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetFillStyle(fillstyle);
      	Energies_ClusterCompletenessEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetLineColor(1);
      	Energies_ClusterCompletenessEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->Draw("SAME");
      	leg->Draw();
      	float Energies_ClusterCompletenessEnergy_Mean =  Energies_ClusterCompletenessEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->GetMean();
      	float Energies_ClusterCompletenessEnergy_RMS =  Energies_ClusterCompletenessEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->GetRMS();
      	Energies_Mean_ClusterCompletenessEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPoint(Energies_Mean_ClusterCompletenessEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN(),Energies_TrueEnergy_Mean,Energies_ClusterCompletenessEnergy_Mean);
      	Energies_RMS_ClusterCompletenessEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPoint(Energies_RMS_ClusterCompletenessEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN(),Energies_TrueEnergy_Mean,Energies_ClusterCompletenessEnergy_RMS);
      	Energies_Mean_ClusterCompletenessEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPointError(Energies_Mean_ClusterCompletenessEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN()-1,Error_on_True_Energy,Energies_ClusterCompletenessEnergy_RMS/TMath::Sqrt(Enteries));

      	Energies_ClusterPurityEnergy_canvasMap[fEnergies[i]][plane_id]->cd();
      	Energies_ClusterPurityEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetDirectory(0);
      	Energies_ClusterPurityEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetFillColor(colours[j]);
      	Energies_ClusterPurityEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetFillStyle(fillstyle);
      	Energies_ClusterPurityEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetLineColor(1);
      	Energies_ClusterPurityEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->Draw("SAME");
      	leg->Draw();
      	float Energies_ClusterPurityEnergy_Mean =  Energies_ClusterPurityEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->GetMean();
      	float Energies_ClusterPurityEnergy_RMS =  Energies_ClusterPurityEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->GetRMS();
      	Energies_Mean_ClusterPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPoint(Energies_Mean_ClusterPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN(),Energies_TrueEnergy_Mean,Energies_ClusterPurityEnergy_Mean);
      	Energies_RMS_ClusterPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPoint(Energies_RMS_ClusterPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN(),Energies_TrueEnergy_Mean,Energies_ClusterPurityEnergy_RMS);
      	Energies_Mean_ClusterPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPointError(Energies_Mean_ClusterPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN()-1,Error_on_True_Energy,Energies_ClusterPurityEnergy_RMS/TMath::Sqrt(Enteries));

      	Energies_ClusterCompletnessHits_canvasMap[fEnergies[i]][plane_id]->cd();
      	Energies_ClusterCompletnessHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetDirectory(0);
      	Energies_ClusterCompletnessHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetFillColor(colours[j]);
      	Energies_ClusterCompletnessHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetFillStyle(fillstyle);
      	Energies_ClusterCompletnessHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetLineColor(1);
      	Energies_ClusterCompletnessHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->Draw("SAME");
      	leg->Draw();
      	float Energies_ClusterCompletnessHits_Mean =  Energies_ClusterCompletnessHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->GetMean();
      	float Energies_ClusterCompletnessHits_RMS =  Energies_ClusterCompletnessHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->GetRMS();
      	Energies_Mean_ClusterCompletnessHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPoint(Energies_Mean_ClusterCompletnessHits_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN(),Energies_TrueEnergy_Mean,Energies_ClusterCompletnessHits_Mean);
      	Energies_RMS_ClusterCompletnessHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPoint(Energies_RMS_ClusterCompletnessHits_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN(),Energies_TrueEnergy_Mean,Energies_ClusterCompletnessHits_RMS);
      	Energies_Mean_ClusterCompletnessHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPointError(Energies_Mean_ClusterCompletnessHits_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN()-1,Error_on_True_Energy,Energies_ClusterCompletnessHits_RMS/TMath::Sqrt(Enteries));

      	Energies_ClusterPurityHits_canvasMap[fEnergies[i]][plane_id]->cd();
      	Energies_ClusterPurityHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetDirectory(0);
      	Energies_ClusterPurityHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetFillColor(colours[j]);
      	Energies_ClusterPurityHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetFillStyle(fillstyle);
      	Energies_ClusterPurityHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetLineColor(1);
      	Energies_ClusterPurityHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->Draw("SAME");
      	leg->Draw();
      	float Energies_ClusterPurityHits_Mean =  Energies_ClusterPurityHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->GetMean();
      	float Energies_ClusterPurityHits_RMS =  Energies_ClusterPurityHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->GetRMS();
      	Energies_Mean_ClusterPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPoint(Energies_Mean_ClusterPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN(),Energies_TrueEnergy_Mean,Energies_ClusterPurityHits_Mean);
      	Energies_RMS_ClusterPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPoint(Energies_RMS_ClusterPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN(),Energies_TrueEnergy_Mean,Energies_ClusterPurityHits_RMS);
      	Energies_Mean_ClusterPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPointError(Energies_Mean_ClusterPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN()-1,Error_on_True_Energy,Energies_ClusterPurityHits_RMS/TMath::Sqrt(Enteries));

      	Energies_ClusterCompPurityEnergy_canvasMap[fEnergies[i]][plane_id]->cd();
      	Energies_ClusterCompPurityEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetDirectory(0);
      	Energies_ClusterCompPurityEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetFillColor(colours[j]);
      	Energies_ClusterCompPurityEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetFillStyle(fillstyle);
      	Energies_ClusterCompPurityEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetLineColor(1);
      	Energies_ClusterCompPurityEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->Draw("SAME");
      	leg->Draw();
      	float Energies_ClusterCompPurityEnergy_Mean =  Energies_ClusterCompPurityEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->GetMean();
      	float Energies_ClusterCompPurityEnergy_RMS =  Energies_ClusterCompPurityEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->GetRMS();
      	Energies_Mean_ClusterCompPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPoint(Energies_Mean_ClusterCompPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN(),Energies_TrueEnergy_Mean,Energies_ClusterCompPurityEnergy_Mean);
      	Energies_RMS_ClusterCompPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPoint(Energies_RMS_ClusterCompPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN(),Energies_TrueEnergy_Mean,Energies_ClusterCompPurityEnergy_RMS);
      	Energies_Mean_ClusterCompPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPointError(Energies_Mean_ClusterCompPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN()-1,Error_on_True_Energy,Energies_ClusterCompPurityEnergy_RMS/TMath::Sqrt(Enteries));

      	Energies_ClusterCompPurityHits_canvasMap[fEnergies[i]][plane_id]->cd();
      	Energies_ClusterCompPurityHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetDirectory(0);
      	Energies_ClusterCompPurityHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetFillColor(colours[j]);
      	Energies_ClusterCompPurityHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetFillStyle(fillstyle);
      	Energies_ClusterCompPurityHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetLineColor(1);
      	Energies_ClusterCompPurityHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->Draw("SAME");
      	leg->Draw();
      	float Energies_ClusterCompPurityHits_Mean =  Energies_ClusterCompPurityHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->GetMean();
      	float Energies_ClusterCompPurityHits_RMS =  Energies_ClusterCompPurityHits_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->GetRMS();
      	Energies_Mean_ClusterCompPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPoint(Energies_Mean_ClusterCompPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN(),Energies_TrueEnergy_Mean,Energies_ClusterCompPurityHits_Mean);
      	Energies_RMS_ClusterCompPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPoint(Energies_RMS_ClusterCompPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN(),Energies_TrueEnergy_Mean,Energies_ClusterCompPurityHits_RMS);
      	Energies_Mean_ClusterCompPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPointError(Energies_Mean_ClusterCompPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN()-1,Error_on_True_Energy,Energies_ClusterCompPurityHits_RMS/TMath::Sqrt(Enteries));
      }
    }
 
    Energies_Mean_ShowerStart_X_Multi->Add(Energies_Mean_ShowerStart_X_GraphMap[fShowerModuleLabels[j]]);
    Energies_RMS_ShowerStart_X_Multi->Add(Energies_RMS_ShowerStart_X_GraphMap[fShowerModuleLabels[j]]);
    Energies_Mean_ShowerStart_Y_Multi->Add(Energies_Mean_ShowerStart_Y_GraphMap[fShowerModuleLabels[j]]);
    Energies_RMS_ShowerStart_Y_Multi->Add(Energies_RMS_ShowerStart_Y_GraphMap[fShowerModuleLabels[j]]);
    Energies_Mean_ShowerStart_Z_Multi->Add(Energies_Mean_ShowerStart_Z_GraphMap[fShowerModuleLabels[j]]);
    Energies_RMS_ShowerStart_Z_Multi->Add(Energies_RMS_ShowerStart_Z_GraphMap[fShowerModuleLabels[j]]);
    Energies_Mean_ShowerLength_Multi->Add(Energies_Mean_ShowerLength_GraphMap[fShowerModuleLabels[j]]);
    Energies_RMS_ShowerLength_Multi->Add(Energies_RMS_ShowerLength_GraphMap[fShowerModuleLabels[j]]);
    Energies_Mean_ShowerEnergyDiff_Multi->Add(Energies_Mean_ShowerEnergyDiff_GraphMap[fShowerModuleLabels[j]]);
    Energies_RMS_ShowerEnergyDiff_Multi->Add(Energies_RMS_ShowerEnergyDiff_GraphMap[fShowerModuleLabels[j]]);
    Energies_Mean_ShowerdEdx_Multi->Add(Energies_Mean_ShowerdEdx_GraphMap[fShowerModuleLabels[j]]);
    Energies_RMS_ShowerdEdx_Multi->Add(Energies_RMS_ShowerdEdx_GraphMap[fShowerModuleLabels[j]]);
    Energies_Mean_EventSeggy_Multi->Add(Energies_Mean_EventSeggy_GraphMap[fShowerModuleLabels[j]]);
    Energies_RMS_EventSeggy_Multi->Add(Energies_RMS_EventSeggy_GraphMap[fShowerModuleLabels[j]]);
    Energies_Mean_ShowerCompleteness_Multi->Add(Energies_Mean_ShowerCompleteness_GraphMap[fShowerModuleLabels[j]]);
    Energies_RMS_ShowerCompleteness_Multi->Add(Energies_RMS_ShowerCompleteness_GraphMap[fShowerModuleLabels[j]]);
    Energies_Mean_ShowerPurity_Multi->Add(Energies_Mean_ShowerPurity_GraphMap[fShowerModuleLabels[j]]);
    Energies_RMS_ShowerPurity_Multi->Add(Energies_RMS_ShowerPurity_GraphMap[fShowerModuleLabels[j]]);
    Energies_Mean_ShowerEnergy_Multi->Add(Energies_Mean_ShowerEnergy_GraphMap[fShowerModuleLabels[j]]);
    Energies_RMS_ShowerEnergy_Multi->Add(Energies_RMS_ShowerEnergy_GraphMap[fShowerModuleLabels[j]]);
    Energies_Mean_ShowerHitNum_Multi->Add(Energies_Mean_ShowerHitNum_GraphMap[fShowerModuleLabels[j]]);
    Energies_RMS_ShowerHitNum_Multi->Add(Energies_RMS_ShowerHitNum_GraphMap[fShowerModuleLabels[j]]);
    Energies_Mean_HitEnergy_Multi->Add(Energies_Mean_HitEnergy_GraphMap[fShowerModuleLabels[j]]);
    Energies_RMS_HitEnergy_Multi->Add(Energies_RMS_HitEnergy_GraphMap[fShowerModuleLabels[j]]);
    Energies_Mean_ShowerMag_Multi->Add(Energies_Mean_ShowerMag_GraphMap[fShowerModuleLabels[j]]);
    Energies_RMS_ShowerMag_Multi->Add(Energies_RMS_ShowerMag_GraphMap[fShowerModuleLabels[j]]);
    Energies_Mean_ShowerDirectionDiff_Multi->Add(Energies_Mean_ShowerDirectionDiff_GraphMap[fShowerModuleLabels[j]]);
    Energies_RMS_ShowerDirectionDiff_Multi->Add(Energies_RMS_ShowerDirectionDiff_GraphMap[fShowerModuleLabels[j]]);
    Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_Multi->Add(Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap[fShowerModuleLabels[j]]);
    Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_Multi->Add(Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap[fShowerModuleLabels[j]]);
    Energies_Mean_TrueEnergy_Multi->Add(Energies_Mean_TrueEnergy_GraphMap[fShowerModuleLabels[j]]);
    Energies_RMS_TrueEnergy_Multi->Add(Energies_RMS_TrueEnergy_GraphMap[fShowerModuleLabels[j]]);
 
    for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){
      Energies_Mean_ClusterProjectionMatchedEnergy_Multi[plane_id]->Add(Energies_Mean_ClusterProjectionMatchedEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]);
      Energies_Mean_ClusterCompletenessEnergy_Multi[plane_id]->Add(Energies_Mean_ClusterCompletenessEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]);
      Energies_Mean_ClusterPurityEnergy_Multi[plane_id]->Add(Energies_Mean_ClusterPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]);
      Energies_Mean_ClusterCompletnessHits_Multi[plane_id]->Add(Energies_Mean_ClusterCompletnessHits_GraphMap[fShowerModuleLabels[j]][plane_id]);
      Energies_Mean_ClusterPurityHits_Multi[plane_id]->Add(Energies_Mean_ClusterPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]);
      Energies_Mean_ClusterCompPurityEnergy_Multi[plane_id]->Add(Energies_Mean_ClusterCompPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]);
      Energies_Mean_ClusterCompPurityHits_Multi[plane_id]->Add(Energies_Mean_ClusterCompPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]);
      Energies_RMS_ClusterProjectionMatchedEnergy_Multi[plane_id]->Add(Energies_RMS_ClusterProjectionMatchedEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]);
      Energies_RMS_ClusterCompletenessEnergy_Multi[plane_id]->Add(Energies_RMS_ClusterCompletenessEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]);
      Energies_RMS_ClusterPurityEnergy_Multi[plane_id]->Add(Energies_RMS_ClusterPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]);
      Energies_RMS_ClusterCompletnessHits_Multi[plane_id]->Add(Energies_RMS_ClusterCompletnessHits_GraphMap[fShowerModuleLabels[j]][plane_id]);
      Energies_RMS_ClusterPurityHits_Multi[plane_id]->Add(Energies_RMS_ClusterPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]);
      Energies_RMS_ClusterCompPurityEnergy_Multi[plane_id]->Add(Energies_RMS_ClusterCompPurityEnergy_GraphMap[fShowerModuleLabels[j]][plane_id]);
      Energies_RMS_ClusterCompPurityHits_Multi[plane_id]->Add(Energies_RMS_ClusterCompPurityHits_GraphMap[fShowerModuleLabels[j]][plane_id]);
    }
  }

  Energies_Mean_ShowerStart_X_canvasMulti->cd();
  Energies_Mean_ShowerStart_X_Multi->Draw("AP");
  leg->Draw();
  Energies_RMS_ShowerStart_X_canvasMulti->cd();
  Energies_RMS_ShowerStart_X_Multi->Draw("AP");
  leg->Draw();

  Energies_Mean_ShowerStart_Y_canvasMulti->cd();
  Energies_Mean_ShowerStart_Y_Multi->Draw("AP");
  leg->Draw();
  Energies_RMS_ShowerStart_Y_canvasMulti->cd();
  Energies_RMS_ShowerStart_Y_Multi->Draw("AP");
  leg->Draw();

  Energies_Mean_ShowerStart_Z_canvasMulti->cd();
  Energies_Mean_ShowerStart_Z_Multi->Draw("AP");
  leg->Draw();
  Energies_RMS_ShowerStart_Z_canvasMulti->cd();
  Energies_RMS_ShowerStart_Z_Multi->Draw("AP");
  leg->Draw();

  Energies_Mean_ShowerLength_canvasMulti->cd();
  Energies_Mean_ShowerLength_Multi->Draw("AP");
  leg->Draw();
  Energies_RMS_ShowerLength_canvasMulti->cd();
  Energies_RMS_ShowerLength_Multi->Draw("AP");
  leg->Draw();

  Energies_Mean_ShowerEnergyDiff_canvasMulti->cd();
  Energies_Mean_ShowerEnergyDiff_Multi->Draw("AP");
  leg->Draw();
  Energies_RMS_ShowerEnergyDiff_canvasMulti->cd();
  Energies_RMS_ShowerEnergyDiff_Multi->Draw("AP");
  leg->Draw();

  Energies_Mean_ShowerdEdx_canvasMulti->cd();
  Energies_Mean_ShowerdEdx_Multi->Draw("AP");
  leg->Draw();
  Energies_RMS_ShowerdEdx_canvasMulti->cd();
  Energies_RMS_ShowerdEdx_Multi->Draw("AP");
  leg->Draw();

  Energies_Mean_EventSeggy_canvasMulti->cd();
  Energies_Mean_EventSeggy_Multi->Draw("AP");
  leg->Draw();
  Energies_RMS_EventSeggy_canvasMulti->cd();
  Energies_RMS_EventSeggy_Multi->Draw("AP");
  leg->Draw();

  Energies_Mean_ShowerCompleteness_canvasMulti->cd();
  Energies_Mean_ShowerCompleteness_Multi->Draw("AP");
  leg->Draw();
  Energies_RMS_ShowerCompleteness_canvasMulti->cd();
  Energies_RMS_ShowerCompleteness_Multi->Draw("AP");
  leg->Draw();

  Energies_Mean_ShowerPurity_canvasMulti->cd();
  Energies_Mean_ShowerPurity_Multi->Draw("AP");
  leg->Draw();
  Energies_RMS_ShowerPurity_canvasMulti->cd();
  Energies_RMS_ShowerPurity_Multi->Draw("AP");
  leg->Draw();

  Energies_Mean_ShowerEnergy_canvasMulti->cd();
  Energies_Mean_ShowerEnergy_Multi->Draw("AP");
  leg->Draw();
  Energies_RMS_ShowerEnergy_canvasMulti->cd();
  Energies_RMS_ShowerEnergy_Multi->Draw("AP");
  leg->Draw();

  Energies_Mean_ShowerHitNum_canvasMulti->cd();
  Energies_Mean_ShowerHitNum_Multi->Draw("AP");
  leg->Draw();
  Energies_RMS_ShowerHitNum_canvasMulti->cd();
  Energies_RMS_ShowerHitNum_Multi->Draw("AP");
  leg->Draw();

  Energies_Mean_HitEnergy_canvasMulti->cd();
  Energies_Mean_HitEnergy_Multi->Draw("AP");
  leg->Draw();
  Energies_RMS_HitEnergy_canvasMulti->cd();
  Energies_RMS_HitEnergy_Multi->Draw("AP");
  leg->Draw();

  Energies_Mean_ShowerMag_canvasMulti->cd();
  Energies_Mean_ShowerMag_Multi->Draw("AP");
  leg->Draw();
  Energies_RMS_ShowerMag_canvasMulti->cd();
  Energies_RMS_ShowerMag_Multi->Draw("AP");
  leg->Draw();

  Energies_Mean_ShowerDirectionDiff_canvasMulti->cd();
  Energies_Mean_ShowerDirectionDiff_Multi->Draw("AP");
  leg->Draw();
  Energies_RMS_ShowerDirectionDiff_canvasMulti->cd();
  Energies_RMS_ShowerDirectionDiff_Multi->Draw("AP");
  leg->Draw();

  Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMulti->cd();
  Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_Multi->Draw("AP");
  leg->Draw();
  Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMulti->cd();
  Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_Multi->Draw("AP");
  leg->Draw();

  Energies_Mean_TrueEnergy_canvasMulti->cd();
  Energies_Mean_TrueEnergy_Multi->Draw("AP");
  leg->Draw();
  Energies_RMS_TrueEnergy_canvasMulti->cd();
  Energies_RMS_TrueEnergy_Multi->Draw("AP");
  leg->Draw();
  
  for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){

    Energies_Mean_ClusterProjectionMatchedEnergy_canvasMulti[plane_id]->cd();
    Energies_Mean_ClusterProjectionMatchedEnergy_Multi[plane_id]->Draw("AP");
    leg->Draw();
    Energies_RMS_ClusterProjectionMatchedEnergy_canvasMulti[plane_id]->cd();
    Energies_RMS_ClusterProjectionMatchedEnergy_Multi[plane_id]->Draw("AP");
    leg->Draw();
    
    Energies_Mean_ClusterCompletenessEnergy_canvasMulti[plane_id]->cd();
    Energies_Mean_ClusterCompletenessEnergy_Multi[plane_id]->Draw("AP");
    leg->Draw();
    Energies_RMS_ClusterCompletenessEnergy_canvasMulti[plane_id]->cd();
    Energies_RMS_ClusterCompletenessEnergy_Multi[plane_id]->Draw("AP");
    leg->Draw();
    
    Energies_Mean_ClusterPurityEnergy_canvasMulti[plane_id]->cd();
    Energies_Mean_ClusterPurityEnergy_Multi[plane_id]->Draw("AP");
    leg->Draw();
    Energies_RMS_ClusterPurityEnergy_canvasMulti[plane_id]->cd();
    Energies_RMS_ClusterPurityEnergy_Multi[plane_id]->Draw("AP");
    leg->Draw();
    
    Energies_Mean_ClusterCompletnessHits_canvasMulti[plane_id]->cd();
    Energies_Mean_ClusterCompletnessHits_Multi[plane_id]->Draw("AP");
    leg->Draw();
    Energies_RMS_ClusterCompletnessHits_canvasMulti[plane_id]->cd();
    Energies_RMS_ClusterCompletnessHits_Multi[plane_id]->Draw("AP");
    leg->Draw();
    
    Energies_Mean_ClusterPurityHits_canvasMulti[plane_id]->cd();
    Energies_Mean_ClusterPurityHits_Multi[plane_id]->Draw("AP");
    leg->Draw();
    Energies_RMS_ClusterPurityHits_canvasMulti[plane_id]->cd();
    Energies_RMS_ClusterPurityHits_Multi[plane_id]->Draw("AP");
    leg->Draw();
    
    Energies_Mean_ClusterCompPurityEnergy_canvasMulti[plane_id]->cd();
    Energies_Mean_ClusterCompPurityEnergy_Multi[plane_id]->Draw("AP");
    leg->Draw();
    Energies_RMS_ClusterCompPurityEnergy_canvasMulti[plane_id]->cd();
    Energies_RMS_ClusterCompPurityEnergy_Multi[plane_id]->Draw("AP");
    leg->Draw();
    
    Energies_Mean_ClusterCompPurityHits_canvasMulti[plane_id]->cd();
    Energies_Mean_ClusterCompPurityHits_Multi[plane_id]->Draw("AP");
    leg->Draw();
    Energies_RMS_ClusterCompPurityHits_canvasMulti[plane_id]->cd();
    Energies_RMS_ClusterCompPurityHits_Multi[plane_id]->Draw("AP");
    leg->Draw();
  }
}


DEFINE_ART_MODULE(ana::ShowerValidation)
