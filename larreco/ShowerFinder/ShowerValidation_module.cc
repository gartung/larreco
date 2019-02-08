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
//              13.11.2018  Hit Validation, 2D Histograms and refactoring of the code                     //
//              22.01.2019  TTre output added                                                             //
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
#include "TH1F.h"
#include "TH2F.h"
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
  
  void InitaliseGraphs(std::string Name, std::string TitleName,std::map<std::string,TH1F*>& Name_HistMap, 
		       std::map<std::string,std::map<float,TH1F*> >& Energies_name_HistMap, 
		       std::map<std::string,TGraphErrors*>& Energies_Mean_Name_GraphMap,
		       std::map<std::string,TGraphErrors*>& Energies_RMS_Name_GraphMap,
		       TMultiGraph*& Energies_Mean_Name_Multi,TMultiGraph*& Energies_RMS_Name_Multi,
		       TCanvas*& Energies_Name_canvasMulti,TCanvas*& Energies_RMS_Name_canvasMulti,
		       std::map<float,TCanvas*>& Energies_Mean_Name_canvasMap, TCanvas*& Name_canvas,
		       std::map<std::string,TH2F*>& Name_2dHistMap,
		       std::map<std::string,TCanvas*>& Name_2dCanvasMap,
		       int x_numbins, float x_start, float x_end,
		       int y_numbins, float y_start, float y_end,
		       std::map<std::string, std::vector<float> >& MetricVector
		       );
  
  void InitaliseGraphs(std::string Name, std::string TitleName,
		       std::map<std::string,std::map<geo::PlaneID,TH1F*> >& Name_HistMap,
		       std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > >& Energies_Name_HistMap,
		       std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> >& Energies_Mean_Name_GraphMap,
		       std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> >& Energies_RMS_Name_GraphMap,
		       std::map<geo::PlaneID,TMultiGraph*>& Energies_Mean_Name_Multi,
		       std::map<geo::PlaneID,TMultiGraph*>& Energies_RMS_Name_Multi,
		       std::map<geo::PlaneID,TCanvas*>& Energies_Mean_Name_canvasMulti,
		       std::map<geo::PlaneID,TCanvas*>& Energies_RMS_Name_canvasMulti,
		       std::map<float,std::map<geo::PlaneID,TCanvas*> >& Energies_Name_canvasMap,
		       std::map<geo::PlaneID,TCanvas*>& Name_canvas,
		       std::map<std::string,std::map<geo::PlaneID,TH2F*> >& Name_2dHistMap,
		       std::map<std::string,std::map<geo::PlaneID,TCanvas*> >& Name_2dCanvasMap,
		       int x_numbins, float x_start, float x_end,
		       int y_numbins, float y_start, float y_end,
		       std::map<std::string,std::vector<std::vector<float> > >& MetricVector
		      );
  
  void InitaliseHitGraphs(std::string Name, std::string TitleName,
			  std::map<std::string,std::map<geo::PlaneID,TH1F*> >& Name_HistMap,
			  std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > >& Energies_Name_HistMap,
			  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> >& Energies_Mean_Name_GraphMap,
			  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> >& Energies_RMS_Name_GraphMap,
			  std::map<geo::PlaneID,TMultiGraph*>& Energies_Mean_Name_Multi,
			  std::map<geo::PlaneID,TMultiGraph*>& Energies_RMS_Name_Multi,
			  std::map<geo::PlaneID,TCanvas*>& Energies_Mean_Name_canvasMulti,
			  std::map<geo::PlaneID,TCanvas*>& Energies_RMS_Name_canvasMulti,
			  std::map<float,std::map<geo::PlaneID,TCanvas*> >& Energies_Name_canvasMap,
			  std::map<geo::PlaneID,TCanvas*>& Name_canvas,
			  std::map<std::string,std::map<geo::PlaneID,TH2F*> >& Name_2dHistMap,
			  std::map<std::string,std::map<geo::PlaneID,TCanvas*> >& Name_2dCanvasMap,
			  int x_numbins, float x_start, float x_end,
			  int y_numbins, float y_start, float y_end,
			  std::map<std::string,std::vector<std::vector<float> > >& MetricVector
			  );


  void DrawGraphs(std::map<std::string,TH1F*>& Name_HistMap, 
		  std::map<std::string,std::map<float,TH1F*> >& Energies_name_HistMap, 
		  std::map<std::string,TGraphErrors*>& Energies_Mean_Name_GraphMap,
		  std::map<std::string,TGraphErrors*>& Energies_RMS_Name_GraphMap,
		  TMultiGraph*& Energies_Mean_Name_Multi,TMultiGraph*& Energies_RMS_Name_Multi,
		  TCanvas*& Energies_Name_canvasMulti,TCanvas*& Energies_RMS_Name_canvasMulti,
		  std::map<float,TCanvas*>& Energies_Mean_Name_canvasMap, TCanvas*& Name_canvas,
		  std::map<std::string,TH2F*>& Name_2dHistMap,
		  std::map<std::string,TCanvas*>& Name_2dCanvasMap
		  );
  
  void DrawGraphs(std::map<std::string,std::map<geo::PlaneID,TH1F*> >& Name_HistMap,
		  std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > >& Energies_Name_HistMap,
		  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> >& Energies_Mean_Name_GraphMap,
		  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> >& Energies_RMS_Name_GraphMap,
		  std::map<geo::PlaneID,TMultiGraph*>& Energies_Mean_Name_Multi,
		  std::map<geo::PlaneID,TMultiGraph*>& Energies_RMS_Name_Multi,
		  std::map<geo::PlaneID,TCanvas*>& Energies_Mean_Name_canvasMulti,
		  std::map<geo::PlaneID,TCanvas*>& Energies_RMS_Name_canvasMulti,
		  std::map<float,std::map<geo::PlaneID,TCanvas*> >& Energies_Mean_Name_canvasMap,
		  std::map<geo::PlaneID,TCanvas*>& Name_canvas,
		  std::map<std::string,std::map<geo::PlaneID,TH2F*> >& Name_2dHistMap,
		  std::map<std::string,std::map<geo::PlaneID,TCanvas*> >& Name_2dCanvasMap
		  );

  void DrawHitGraphs(std::map<std::string,std::map<geo::PlaneID,TH1F*> >& Name_HistMap,
		     std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > >& Energies_Name_HistMap,
		     std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> >& Energies_Mean_Name_GraphMap,
		     std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> >& Energies_RMS_Name_GraphMap,
		     std::map<geo::PlaneID,TMultiGraph*>& Energies_Mean_Name_Multi,
		     std::map<geo::PlaneID,TMultiGraph*>& Energies_RMS_Name_Multi,
		     std::map<geo::PlaneID,TCanvas*>& Energies_Mean_Name_canvasMulti,
		     std::map<geo::PlaneID,TCanvas*>& Energies_RMS_Name_canvasMulti,
		     std::map<float,std::map<geo::PlaneID,TCanvas*> >& Energies_Mean_Name_canvasMap,
		     std::map<geo::PlaneID,TCanvas*>& Name_canvas,
		     std::map<std::string,std::map<geo::PlaneID,TH2F*> >& Name_2dHistMap,
		     std::map<std::string,std::map<geo::PlaneID,TCanvas*> >& Name_2dCanvasMap
		     );
    


  void ClusterValidation(std::vector<art::Ptr<recob::Cluster> >& clusters, 
			 const art::Event& evt, art::Handle<std::vector<recob::Cluster> >& clusterHandle, 
			 std::map<int,std::vector<int> >& ShowerMotherTrackIDs,
			 std::map<int,float>& MCTrack_Energy_map, 
			 std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > > & MCTrack_hit_map, 
			 int& TrueShowerID,
			 float& simenergy, 
			 std::string & fShowerModuleLabel
			 );

private:

  //fcl parameters 
  std::string fGenieGenModuleLabel;
  std::string fLArGeantModuleLabel;
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;

  bool  fUseBiggestShower;
  bool  fDrawCanvases;
  bool  fFillOnlyClosestShower;
  int   fVerbose;
  int   fMinHitSize;
  float fEnergyWidth;
  float fSimEnergyCut;
  float fDensityCut;
  float fMaxSimEnergy;

  std::vector<float>       fEnergies;
  std::vector<std::string> fShowerModuleLabels;
  std::vector<std::string> fHitModuleLabels;

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
  std::map<std::string,TH1F*> ShowerEnergyCompleteness_HistMap;
  std::map<std::string,TH1F*> ShowerEnergyPurity_HistMap;
  std::map<std::string,TH1F*> ShowerHitsCompleteness_HistMap;
  std::map<std::string,TH1F*> ShowerHitsPurity_HistMap;
  std::map<std::string,TH1F*> ShowerEnergy_HistMap;
  std::map<std::string,TH1F*> ShowerHitNum_HistMap;
  std::map<std::string,TH1F*> ShowerTotalEnergyDiff_HistMap;
  std::map<std::string,TH1F*> ShowerMag_HistMap;
  std::map<std::string,TH1F*> ShowerDirectionDiff_HistMap;
  std::map<std::string,TH1F*> ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap;
  std::map<std::string,TH1F*> ShowerTrueEnergy_HistMap;
  std::map<std::string,TH1F*> TrueEnergy_HistMap;
  std::map<std::string,TH1F*> TrueHitNum_HistMap;
  std::map<std::string,TH1F*> ShowerBestPlane_HistMap;
  std::map<std::string,TH1F*> GeoProjectionMatched_HistMap;

  std::map<std::string,std::map<geo::PlaneID,TH1F*> > ClusterProjectionMatchedEnergy_HistMap;
  std::map<std::string,std::map<geo::PlaneID,TH1F*> > ClusterCompletenessEnergy_HistMap;
  std::map<std::string,std::map<geo::PlaneID,TH1F*> > ClusterPurityEnergy_HistMap;
  std::map<std::string,std::map<geo::PlaneID,TH1F*> > ClusterCompletenessHits_HistMap;
  std::map<std::string,std::map<geo::PlaneID,TH1F*> > ClusterPurityHits_HistMap;
  std::map<std::string,std::map<geo::PlaneID,TH1F*> > ClusterCompPurityEnergy_HistMap;
  std::map<std::string,std::map<geo::PlaneID,TH1F*> > ClusterCompPurityHits_HistMap;
  
  std::map<std::string,std::map<geo::PlaneID,TH1F*> > HitCompletenessEnergy_HistMap;

  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerDirection_X_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerDirection_Y_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerDirection_Z_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerStart_X_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerStart_Y_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerStart_Z_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerLength_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerEnergyDiff_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerdEdx_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_EventSeggy_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerEnergyCompleteness_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerEnergyPurity_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerHitsCompleteness_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerHitsPurity_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerEnergy_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerHitNum_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerTotalEnergyDiff_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerMag_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerDirectionDiff_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerTrueEnergy_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_TrueEnergy_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_TrueHitNum_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_ShowerBestPlane_HistMap;
  std::map<std::string,std::map<float,TH1F*> > Energies_GeoProjectionMatched_HistMap;

  std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > > Energies_ClusterProjectionMatchedEnergy_HistMap;
  std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > > Energies_ClusterCompletenessEnergy_HistMap;
  std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > > Energies_ClusterPurityEnergy_HistMap;
  std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > > Energies_ClusterCompletenessHits_HistMap;
  std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > > Energies_ClusterPurityHits_HistMap;
  std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > > Energies_ClusterCompPurityEnergy_HistMap;
  std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > > Energies_ClusterCompPurityHits_HistMap;
  
  std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > > Energies_HitCompletenessEnergy_HistMap;

  
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerDirection_X_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerDirection_Y_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerDirection_Z_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerStart_X_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerStart_Y_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerStart_Z_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerLength_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerEnergyDiff_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerdEdx_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_EventSeggy_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerEnergyCompleteness_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerEnergyPurity_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerHitsCompleteness_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerHitsPurity_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerEnergy_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerHitNum_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerTotalEnergyDiff_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerMag_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerDirectionDiff_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerTrueEnergy_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_TrueEnergy_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_TrueHitNum_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_ShowerBestPlane_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_Mean_GeoProjectionMatched_GraphMap;

  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_Mean_ClusterProjectionMatchedEnergy_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_Mean_ClusterCompletenessEnergy_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_Mean_ClusterPurityEnergy_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_Mean_ClusterCompletenessHits_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_Mean_ClusterPurityHits_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_Mean_ClusterCompPurityEnergy_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_Mean_ClusterCompPurityHits_GraphMap;

  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_Mean_HitCompletenessEnergy_GraphMap;

  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerDirection_X_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerDirection_Y_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerDirection_Z_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerStart_X_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerStart_Y_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerStart_Z_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerLength_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerEnergyDiff_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerdEdx_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_EventSeggy_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerEnergyCompleteness_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerEnergyPurity_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerHitsCompleteness_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerHitsPurity_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerEnergy_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerHitNum_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerTotalEnergyDiff_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerMag_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerDirectionDiff_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerTrueEnergy_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_TrueEnergy_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_TrueHitNum_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_ShowerBestPlane_GraphMap;
  std::map<std::string,TGraphErrors*> Energies_RMS_GeoProjectionMatched_GraphMap;

  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_RMS_ClusterProjectionMatchedEnergy_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_RMS_ClusterCompletenessEnergy_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_RMS_ClusterPurityEnergy_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_RMS_ClusterCompletenessHits_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_RMS_ClusterPurityHits_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_RMS_ClusterCompPurityEnergy_GraphMap;
  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_RMS_ClusterCompPurityHits_GraphMap;

  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> > Energies_RMS_HitCompletenessEnergy_GraphMap;


  TMultiGraph* Energies_Mean_ShowerDirection_X_Multi;
  TMultiGraph* Energies_Mean_ShowerDirection_Y_Multi;
  TMultiGraph* Energies_Mean_ShowerDirection_Z_Multi;
  TMultiGraph* Energies_Mean_ShowerStart_X_Multi;
  TMultiGraph* Energies_Mean_ShowerStart_Y_Multi;
  TMultiGraph* Energies_Mean_ShowerStart_Z_Multi;
  TMultiGraph* Energies_Mean_ShowerLength_Multi;
  TMultiGraph* Energies_Mean_ShowerEnergyDiff_Multi;
  TMultiGraph* Energies_Mean_ShowerdEdx_Multi;
  TMultiGraph* Energies_Mean_EventSeggy_Multi;
  TMultiGraph* Energies_Mean_ShowerEnergyCompleteness_Multi;
  TMultiGraph* Energies_Mean_ShowerEnergyPurity_Multi;
  TMultiGraph* Energies_Mean_ShowerHitsCompleteness_Multi;
  TMultiGraph* Energies_Mean_ShowerHitsPurity_Multi;
  TMultiGraph* Energies_Mean_ShowerEnergy_Multi;
  TMultiGraph* Energies_Mean_ShowerHitNum_Multi;
  TMultiGraph* Energies_Mean_ShowerTotalEnergyDiff_Multi;
  TMultiGraph* Energies_Mean_ShowerMag_Multi;
  TMultiGraph* Energies_Mean_ShowerDirectionDiff_Multi;
  TMultiGraph* Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_Multi;
  TMultiGraph* Energies_Mean_ShowerTrueEnergy_Multi;
  TMultiGraph* Energies_Mean_TrueEnergy_Multi;
  TMultiGraph* Energies_Mean_TrueHitNum_Multi;
  TMultiGraph* Energies_Mean_ShowerBestPlane_Multi;
  TMultiGraph* Energies_Mean_GeoProjectionMatched_Multi;

  std::map<geo::PlaneID,TMultiGraph*> Energies_Mean_ClusterProjectionMatchedEnergy_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_Mean_ClusterCompletenessEnergy_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_Mean_ClusterPurityEnergy_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_Mean_ClusterCompletenessHits_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_Mean_ClusterPurityHits_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_Mean_ClusterCompPurityEnergy_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_Mean_ClusterCompPurityHits_Multi;
  
  std::map<geo::PlaneID,TMultiGraph*> Energies_Mean_HitCompletenessEnergy_Multi;

  TMultiGraph* Energies_RMS_ShowerDirection_X_Multi;
  TMultiGraph* Energies_RMS_ShowerDirection_Y_Multi;
  TMultiGraph* Energies_RMS_ShowerDirection_Z_Multi;
  TMultiGraph* Energies_RMS_ShowerStart_X_Multi;
  TMultiGraph* Energies_RMS_ShowerStart_Y_Multi;
  TMultiGraph* Energies_RMS_ShowerStart_Z_Multi;
  TMultiGraph* Energies_RMS_ShowerLength_Multi;
  TMultiGraph* Energies_RMS_ShowerEnergyDiff_Multi;
  TMultiGraph* Energies_RMS_ShowerdEdx_Multi;
  TMultiGraph* Energies_RMS_EventSeggy_Multi;
  TMultiGraph* Energies_RMS_ShowerEnergyCompleteness_Multi;
  TMultiGraph* Energies_RMS_ShowerEnergyPurity_Multi;
  TMultiGraph* Energies_RMS_ShowerHitsCompleteness_Multi;
  TMultiGraph* Energies_RMS_ShowerHitsPurity_Multi;
  TMultiGraph* Energies_RMS_ShowerEnergy_Multi;
  TMultiGraph* Energies_RMS_ShowerHitNum_Multi;
  TMultiGraph* Energies_RMS_ShowerTotalEnergyDiff_Multi;
  TMultiGraph* Energies_RMS_ShowerMag_Multi;
  TMultiGraph* Energies_RMS_ShowerDirectionDiff_Multi;
  TMultiGraph* Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_Multi;
  TMultiGraph* Energies_RMS_ShowerTrueEnergy_Multi;
  TMultiGraph* Energies_RMS_TrueEnergy_Multi;
  TMultiGraph* Energies_RMS_TrueHitNum_Multi;
  TMultiGraph* Energies_RMS_ShowerBestPlane_Multi;
  TMultiGraph* Energies_RMS_GeoProjectionMatched_Multi;

  std::map<geo::PlaneID,TMultiGraph*> Energies_RMS_ClusterProjectionMatchedEnergy_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_RMS_ClusterCompletenessEnergy_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_RMS_ClusterPurityEnergy_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_RMS_ClusterCompletenessHits_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_RMS_ClusterPurityHits_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_RMS_ClusterCompPurityEnergy_Multi;
  std::map<geo::PlaneID,TMultiGraph*> Energies_RMS_ClusterCompPurityHits_Multi;

  std::map<geo::PlaneID,TMultiGraph*> Energies_RMS_HitCompletenessEnergy_Multi;

  TCanvas* Energies_Mean_ShowerDirection_X_canvasMulti;
  TCanvas* Energies_Mean_ShowerDirection_Y_canvasMulti;
  TCanvas* Energies_Mean_ShowerDirection_Z_canvasMulti;
  TCanvas* Energies_Mean_ShowerStart_X_canvasMulti;
  TCanvas* Energies_Mean_ShowerStart_Y_canvasMulti;
  TCanvas* Energies_Mean_ShowerStart_Z_canvasMulti;
  TCanvas* Energies_Mean_ShowerLength_canvasMulti;
  TCanvas* Energies_Mean_ShowerEnergyDiff_canvasMulti;
  TCanvas* Energies_Mean_ShowerTotalEnergyDiff_canvasMulti;
  TCanvas* Energies_Mean_ShowerdEdx_canvasMulti;
  TCanvas* Energies_Mean_EventSeggy_canvasMulti;
  TCanvas* Energies_Mean_ShowerEnergyCompleteness_canvasMulti;
  TCanvas* Energies_Mean_ShowerEnergyPurity_canvasMulti;
  TCanvas* Energies_Mean_ShowerHitsCompleteness_canvasMulti;
  TCanvas* Energies_Mean_ShowerHitsPurity_canvasMulti;
  TCanvas* Energies_Mean_ShowerEnergy_canvasMulti;
  TCanvas* Energies_Mean_ShowerHitNum_canvasMulti;
  TCanvas* Energies_Mean_ShowerDirectionDiff_canvasMulti;
  TCanvas* Energies_Mean_ShowerMag_canvasMulti;
  TCanvas* Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMulti;
  TCanvas* Energies_Mean_ShowerTrueEnergy_canvasMulti;
  TCanvas* Energies_Mean_TrueEnergy_canvasMulti;
  TCanvas* Energies_Mean_TrueHitNum_canvasMulti;
  TCanvas* Energies_Mean_ShowerBestPlane_canvasMulti;
  TCanvas* Energies_Mean_GeoProjectionMatched_canvasMulti;

  std::map<geo::PlaneID,TCanvas*> Energies_Mean_ClusterProjectionMatchedEnergy_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_Mean_ClusterCompletenessEnergy_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_Mean_ClusterPurityEnergy_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_Mean_ClusterCompletenessHits_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_Mean_ClusterPurityHits_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_Mean_ClusterCompPurityEnergy_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_Mean_ClusterCompPurityHits_canvasMulti;

  std::map<geo::PlaneID,TCanvas*> Energies_Mean_HitCompletenessEnergy_canvasMulti;

  TCanvas* Energies_RMS_ShowerDirection_X_canvasMulti;
  TCanvas* Energies_RMS_ShowerDirection_Y_canvasMulti;
  TCanvas* Energies_RMS_ShowerDirection_Z_canvasMulti;
  TCanvas* Energies_RMS_ShowerStart_X_canvasMulti;
  TCanvas* Energies_RMS_ShowerStart_Y_canvasMulti;
  TCanvas* Energies_RMS_ShowerStart_Z_canvasMulti;
  TCanvas* Energies_RMS_ShowerLength_canvasMulti;
  TCanvas* Energies_RMS_ShowerEnergyDiff_canvasMulti;
  TCanvas* Energies_RMS_ShowerTotalEnergyDiff_canvasMulti;
  TCanvas* Energies_RMS_ShowerdEdx_canvasMulti;
  TCanvas* Energies_RMS_EventSeggy_canvasMulti;
  TCanvas* Energies_RMS_ShowerEnergyCompleteness_canvasMulti;
  TCanvas* Energies_RMS_ShowerEnergyPurity_canvasMulti;
  TCanvas* Energies_RMS_ShowerHitsCompleteness_canvasMulti;
  TCanvas* Energies_RMS_ShowerHitsPurity_canvasMulti;
  TCanvas* Energies_RMS_ShowerEnergy_canvasMulti;
  TCanvas* Energies_RMS_ShowerHitNum_canvasMulti;
  TCanvas* Energies_RMS_ShowerDirectionDiff_canvasMulti;
  TCanvas* Energies_RMS_ShowerMag_canvasMulti;
  TCanvas* Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMulti;
  TCanvas* Energies_RMS_ShowerTrueEnergy_canvasMulti;
  TCanvas* Energies_RMS_TrueEnergy_canvasMulti;
  TCanvas* Energies_RMS_TrueHitNum_canvasMulti;
  TCanvas* Energies_RMS_ShowerBestPlane_canvasMulti;
  TCanvas* Energies_RMS_GeoProjectionMatched_canvasMulti;

  std::map<geo::PlaneID,TCanvas*> Energies_RMS_ClusterProjectionMatchedEnergy_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_RMS_ClusterCompletenessEnergy_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_RMS_ClusterPurityEnergy_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_RMS_ClusterCompletenessHits_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_RMS_ClusterPurityHits_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_RMS_ClusterCompPurityEnergy_canvasMulti;
  std::map<geo::PlaneID,TCanvas*> Energies_RMS_ClusterCompPurityHits_canvasMulti;

  std::map<geo::PlaneID,TCanvas*> Energies_RMS_HitCompletenessEnergy_canvasMulti;


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
  TCanvas* ShowerEnergyCompleteness_canvas;
  TCanvas* ShowerEnergyPurity_canvas;
  TCanvas* ShowerHitsCompleteness_canvas;
  TCanvas* ShowerHitsPurity_canvas;
  TCanvas* ShowerEnergy_canvas;
  TCanvas* ShowerHitNum_canvas;
  TCanvas* ShowerDirectionDiff_canvas;
  TCanvas* ShowerMag_canvas;
  TCanvas* ShowerRecoEnergyVsTrueEnergyinRecoShower_canvas;
  TCanvas* ShowerTrueEnergy_canvas;
  TCanvas* TrueEnergy_canvas;
  TCanvas* TrueHitNum_canvas;
  TCanvas* ShowerBestPlane_canvas;
  TCanvas* GeoProjectionMatched_canvas;

  std::map<geo::PlaneID,TCanvas*> ClusterProjectionMatchedEnergy_canvas;
  std::map<geo::PlaneID,TCanvas*> ClusterCompletenessEnergy_canvas;
  std::map<geo::PlaneID,TCanvas*> ClusterPurityEnergy_canvas;
  std::map<geo::PlaneID,TCanvas*> ClusterCompletenessHits_canvas;
  std::map<geo::PlaneID,TCanvas*> ClusterPurityHits_canvas;
  std::map<geo::PlaneID,TCanvas*> ClusterCompPurityEnergy_canvas;
  std::map<geo::PlaneID,TCanvas*> ClusterCompPurityHits_canvas;
  
  std::map<geo::PlaneID,TCanvas*> HitCompletenessEnergy_canvas;

  std::map<float,TCanvas*> Energies_ShowerDirection_X_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerDirection_Y_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerDirection_Z_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerStart_X_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerStart_Y_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerStart_Z_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerLength_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerEnergyDiff_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerTotalEnergyDiff_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerdEdx_canvasMap;
  std::map<float,TCanvas*> Energies_EventSeggy_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerEnergyCompleteness_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerHitsPurity_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerHitsCompleteness_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerEnergyPurity_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerEnergy_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerHitNum_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerDirectionDiff_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerMag_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerTrueEnergy_canvasMap;
  std::map<float,TCanvas*> Energies_TrueEnergy_canvasMap;
  std::map<float,TCanvas*> Energies_TrueHitNum_canvasMap;
  std::map<float,TCanvas*> Energies_ShowerBestPlane_canvasMap;
  std::map<float,TCanvas*> Energies_GeoProjectionMatched_canvasMap;

  std::map<float,std::map<geo::PlaneID,TCanvas*> > Energies_ClusterProjectionMatchedEnergy_canvasMap;
  std::map<float,std::map<geo::PlaneID,TCanvas*> > Energies_ClusterCompletenessEnergy_canvasMap;
  std::map<float,std::map<geo::PlaneID,TCanvas*> > Energies_ClusterPurityEnergy_canvasMap;
  std::map<float,std::map<geo::PlaneID,TCanvas*> > Energies_ClusterCompletenessHits_canvasMap;
  std::map<float,std::map<geo::PlaneID,TCanvas*> > Energies_ClusterPurityHits_canvasMap;
  std::map<float,std::map<geo::PlaneID,TCanvas*> > Energies_ClusterCompPurityEnergy_canvasMap;
  std::map<float,std::map<geo::PlaneID,TCanvas*> > Energies_ClusterCompPurityHits_canvasMap;

  std::map<float,std::map<geo::PlaneID,TCanvas*> > Energies_HitCompletenessEnergy_canvasMap;

  
  std::map<std::string,TH2F*> ShowerDirection_X_2dHistMap;
  std::map<std::string,TH2F*> ShowerDirection_Y_2dHistMap;
  std::map<std::string,TH2F*> ShowerDirection_Z_2dHistMap;
  std::map<std::string,TH2F*> ShowerStart_X_2dHistMap;
  std::map<std::string,TH2F*> ShowerStart_Y_2dHistMap;
  std::map<std::string,TH2F*> ShowerStart_Z_2dHistMap;
  std::map<std::string,TH2F*> ShowerLength_2dHistMap;
  std::map<std::string,TH2F*> ShowerEnergyDiff_2dHistMap;
  std::map<std::string,TH2F*> ShowerdEdx_2dHistMap;
  std::map<std::string,TH2F*> EventSeggy_2dHistMap;
  std::map<std::string,TH2F*> ShowerEnergyCompleteness_2dHistMap;
  std::map<std::string,TH2F*> ShowerEnergyPurity_2dHistMap;
  std::map<std::string,TH2F*> ShowerHitsCompleteness_2dHistMap;
  std::map<std::string,TH2F*> ShowerHitsPurity_2dHistMap;
  std::map<std::string,TH2F*> ShowerEnergy_2dHistMap;
  std::map<std::string,TH2F*> ShowerHitNum_2dHistMap;
  std::map<std::string,TH2F*> ShowerTotalEnergyDiff_2dHistMap;
  std::map<std::string,TH2F*> ShowerMag_2dHistMap;
  std::map<std::string,TH2F*> ShowerDirectionDiff_2dHistMap;
  std::map<std::string,TH2F*> ShowerRecoEnergyVsTrueEnergyinRecoShower_2dHistMap;
  std::map<std::string,TH2F*> ShowerTrueEnergy_2dHistMap;
  std::map<std::string,TH2F*> TrueEnergy_2dHistMap;
  std::map<std::string,TH2F*> TrueHitNum_2dHistMap;
  std::map<std::string,TH2F*> ShowerBestPlane_2dHistMap;
  std::map<std::string,TH2F*> GeoProjectionMatched_2dHistMap;

  std::map<std::string,std::map<geo::PlaneID,TH2F*> > ClusterProjectionMatchedEnergy_2dHistMap;
  std::map<std::string,std::map<geo::PlaneID,TH2F*> > ClusterCompletenessEnergy_2dHistMap;
  std::map<std::string,std::map<geo::PlaneID,TH2F*> > ClusterPurityEnergy_2dHistMap;
  std::map<std::string,std::map<geo::PlaneID,TH2F*> > ClusterCompletenessHits_2dHistMap;
  std::map<std::string,std::map<geo::PlaneID,TH2F*> > ClusterPurityHits_2dHistMap;
  std::map<std::string,std::map<geo::PlaneID,TH2F*> > ClusterCompPurityEnergy_2dHistMap;
  std::map<std::string,std::map<geo::PlaneID,TH2F*> > ClusterCompPurityHits_2dHistMap;
  
  std::map<std::string,std::map<geo::PlaneID,TH2F*> > HitCompletenessEnergy_2dHistMap;
  
  std::map<std::string,TCanvas*> ShowerDirection_X_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerDirection_Y_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerDirection_Z_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerStart_X_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerStart_Y_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerStart_Z_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerLength_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerEnergyDiff_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerdEdx_2dCanvasMap;
  std::map<std::string,TCanvas*> EventSeggy_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerEnergyCompleteness_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerHitsPurity_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerHitsCompleteness_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerEnergyPurity_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerEnergy_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerHitNum_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerTotalEnergyDiff_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerMag_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerDirectionDiff_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerRecoEnergyVsTrueEnergyinRecoShower_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerTrueEnergy_2dCanvasMap;
  std::map<std::string,TCanvas*> TrueEnergy_2dCanvasMap;
  std::map<std::string,TCanvas*> TrueHitNum_2dCanvasMap;
  std::map<std::string,TCanvas*> ShowerBestPlane_2dCanvasMap;
  std::map<std::string,TCanvas*> GeoProjectionMatched_2dCanvasMap;

  std::map<std::string,TCanvas*> PosDir_2dCanvasMap;
  std::map<std::string,TH2F*>    PosDir_2dHistMap;

  std::map<std::string,std::map<geo::PlaneID,TCanvas*> > ClusterProjectionMatchedEnergy_2dCanvasMap;
  std::map<std::string,std::map<geo::PlaneID,TCanvas*> > ClusterCompletenessEnergy_2dCanvasMap;
  std::map<std::string,std::map<geo::PlaneID,TCanvas*> > ClusterPurityEnergy_2dCanvasMap;
  std::map<std::string,std::map<geo::PlaneID,TCanvas*> > ClusterCompletenessHits_2dCanvasMap;
  std::map<std::string,std::map<geo::PlaneID,TCanvas*> > ClusterPurityHits_2dCanvasMap;
  std::map<std::string,std::map<geo::PlaneID,TCanvas*> > ClusterCompPurityEnergy_2dCanvasMap;
  std::map<std::string,std::map<geo::PlaneID,TCanvas*> > ClusterCompPurityHits_2dCanvasMap;
  
  std::map<std::string,std::map<geo::PlaneID,TCanvas*> > HitCompletenessEnergy_2dCanvasMap;

  std::map<std::string,std::vector<float> > ShowerDirection_X_TreeVal;
  std::map<std::string,std::vector<float> > ShowerDirection_Y_TreeVal;
  std::map<std::string,std::vector<float> > ShowerDirection_Z_TreeVal;
  std::map<std::string,std::vector<float> > ShowerStart_X_TreeVal;
  std::map<std::string,std::vector<float> > ShowerStart_Y_TreeVal;
  std::map<std::string,std::vector<float> > ShowerStart_Z_TreeVal;
  std::map<std::string,std::vector<float> > ShowerLength_TreeVal;
  std::map<std::string,std::vector<float> > ShowerEnergyDiff_TreeVal;
  std::map<std::string,std::vector<float> > ShowerdEdx_TreeVal;
  std::map<std::string,std::vector<float> > EventSeggy_TreeVal;
  std::map<std::string,std::vector<float> > ShowerEnergyCompleteness_TreeVal;
  std::map<std::string,std::vector<float> > ShowerHitsPurity_TreeVal;
  std::map<std::string,std::vector<float> > ShowerHitsCompleteness_TreeVal;
  std::map<std::string,std::vector<float> > ShowerEnergyPurity_TreeVal;
  std::map<std::string,std::vector<float> > ShowerEnergy_TreeVal;
  std::map<std::string,std::vector<float> > ShowerHitNum_TreeVal;
  std::map<std::string,std::vector<float> > ShowerTotalEnergyDiff_TreeVal;
  std::map<std::string,std::vector<float> > ShowerMag_TreeVal;
  std::map<std::string,std::vector<float> > ShowerDirectionDiff_TreeVal;
  std::map<std::string,std::vector<float> > ShowerRecoEnergyVsTrueEnergyinRecoShower_TreeVal;
  std::map<std::string,std::vector<float> > ShowerTrueEnergy_TreeVal;
  std::map<std::string,std::vector<float> > TrueEnergy_TreeVal;
  std::map<std::string,std::vector<float> > TrueHitNum_TreeVal;
  std::map<std::string,std::vector<float> > ShowerBestPlane_TreeVal;
  std::map<std::string,std::vector<float> > GeoProjectionMatched_TreeVal;
  
  std::map<std::string,std::vector<float> > PosDir_TreeVal;

  std::map<std::string,std::vector<std::string> > ShowerStartEndProcess_TreeVal;


  //      fShowerModule RecoShower   NClusters         Plane  Value
  //  std::map<std::string,std::vector<std::vector<std::pair<int,float> > > > 

  std::map<std::string,std::vector<std::vector<float> > >  ClusterProjectionMatchedEnergy_TreeVal;
  std::map<std::string,std::vector<std::vector<float> > >  ClusterCompletenessEnergy_TreeVal;
  std::map<std::string,std::vector<std::vector<float> > >  ClusterPurityEnergy_TreeVal;
  std::map<std::string,std::vector<std::vector<float> > >  ClusterCompletenessHits_TreeVal;
  std::map<std::string,std::vector<std::vector<float> > >  ClusterPurityHits_TreeVal;
  std::map<std::string,std::vector<std::vector<float> > >  ClusterCompPurityEnergy_TreeVal;
  std::map<std::string,std::vector<std::vector<float> > >  ClusterCompPurityHits_TreeVal;

  std::map<std::string,std::vector<std::vector<float> > > HitCompletenessEnergy_TreeVal;


  TH1F* BiggestShowerMotherE_AfterCutHist;
  TH1F* BiggestShowerMotherE_BeforeCutHist;
  TH1F* AsscoiatedBiggestShowerMotherEHist;
  TCanvas* BiggestShowerMotherE_AfterCutCanvas; 
  TCanvas* BiggestShowerMotherE_BeforeCutCanvas;
  TCanvas* AsscoiatedBiggestShowerMotherECanvas;

  TH1* SmallestShowerMotherE_AfterCutHist;
  TH1* SmallestShowerMotherE_BeforeCutHist;
  TH1* AsscoiatedSmallestShowerMotherEHist;
  TCanvas* SmallestShowerMotherE_AfterCutCanvas; 
  TCanvas* SmallestShowerMotherE_BeforeCutCanvas;
  TCanvas* AsscoiatedSmallestShowerMotherECanvas;

  std::vector<float> SmallestShowerMotherE_BeforeCut_TreeVal;
  std::vector<float> BiggestShowerMotherE_BeforeCut_TreeVal;
  std::vector<float> SmallestShowerMotherE_AfterCut_TreeVal;
  std::vector<float> BiggestShowerMotherE_AfterCut_TreeVal;

  std::vector<float> EventNum_name_TreeVal;

  //TTree 
  TTree* Tree;

  //Service handles   
  art::ServiceHandle<cheat::BackTrackerService> backtracker;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<art::TFileService> tfs;


};

ana::ShowerValidation::ShowerValidation(const fhicl::ParameterSet& pset) : EDAnalyzer(pset){

  fGenieGenModuleLabel   = pset.get<std::string>("GenieGenModuleLabel");
  fLArGeantModuleLabel   = pset.get<std::string>("LArGeantModuleLabel"); 
  fHitsModuleLabel       = pset.get<std::string>("HitsModuleLabel");
  fTrackModuleLabel      = pset.get<std::string>("TrackModuleLabel");
  fShowerModuleLabels    = pset.get<std::vector<std::string> >("ShowerModuleLabels");
  fHitModuleLabels       = pset.get<std::vector<std::string> >("HitModuleLabels");
  fEnergies              = pset.get<std::vector<float> >("Energies");
  fUseBiggestShower      = pset.get<bool>("UseBiggestShower");
  fDrawCanvases          = pset.get<bool>("DrawCanvases"); 
  fFillOnlyClosestShower = pset.get<bool>("FillOnlyClosestShower");
  fVerbose               = pset.get<int>("Verbose");
  fMinHitSize            = pset.get<int>("MinHitSize");
  fEnergyWidth           = pset.get<float>("EnergyWidth");
  fSimEnergyCut          = pset.get<float>("SimEnergyCut");
  fDensityCut            = pset.get<float>("DensityCut");
  fMaxSimEnergy          = pset.get<float>("MaxSimEnergy");
}

//This is bad. try and think of a better way. 
void ana::ShowerValidation::InitaliseGraphs(std::string Name, std::string TitleName,std::map<std::string,TH1F*>& Name_HistMap, 
					    std::map<std::string,std::map<float,TH1F*> >& Energies_Name_HistMap, 
					    std::map<std::string,TGraphErrors*>& Energies_Mean_Name_GraphMap,
					    std::map<std::string,TGraphErrors*>& Energies_RMS_Name_GraphMap,
					    TMultiGraph*& Energies_Mean_Name_Multi,
					    TMultiGraph*& Energies_RMS_Name_Multi,
					    TCanvas*& Energies_Mean_Name_canvasMulti,
					    TCanvas*& Energies_RMS_Name_canvasMulti,
					    std::map<float,TCanvas*>& Energies_Name_canvasMap, 
					    TCanvas* & Name_canvas,
					    std::map<std::string,TH2F*>& Name_2dHistMap,
					    std::map<std::string,TCanvas*>& Name_2dCanvasMap,
					    int x_numbins, float x_start, float x_end,
					    int y_numbins, float y_start, float y_end,
					    std::map<std::string,std::vector<float> >& MetricVector){
  
  gStyle->SetOptStat(0);

  //Create The histogram canvas for Each Energy. 
  for(unsigned int i=0; i<fEnergies.size(); ++i){
    std::string Energies_string = Name + std::to_string(fEnergies[i]);
    std::string Energies_titlestring = ";" + TitleName + ";Enteries";
    const char* Energies_name   = Energies_string.c_str();
    const char* Energies_titlename   = Energies_titlestring.c_str();;  
    Energies_Name_canvasMap[fEnergies[i]] = tfs->makeAndRegister<TCanvas>(Energies_name,Energies_titlename);
  }

  for(unsigned int j=0; j<fShowerModuleLabels.size(); ++j){

    //Create the Tree Branch
    std::string Tree_string = Name + "_" + fShowerModuleLabels[j];
    const char* Tree_name   = Tree_string.c_str();
    Tree->Branch(Tree_name,"std::vector<float>", &MetricVector[fShowerModuleLabels[j]], 32000, 0);

    //Create the Hitogram which has all the energies in
    std::string Name_string      = Name + fShowerModuleLabels[j];
    std::string Name_titlestring = fShowerModuleLabels[j] + ";" + TitleName + ";Enteries";
    const char* Name_titlename   = Name_titlestring.c_str();
    const char* Name_name     =  Name_string.c_str();
    Name_HistMap[fShowerModuleLabels[j]]  = tfs->make<TH1F>(Name_name, Name_titlename,x_numbins,x_start,x_end);
    
    //Create the 2D Histogram of Energy and value 
    std::string Name_2d_string      = Name + "_2d_" + fShowerModuleLabels[j];
    std::string Name_2d_titlestring = ";" + TitleName + ";Energy (MeV);Enteries";
    const char* Name_2d_titlename   = Name_2d_titlestring.c_str();
    const char* Name_2d_name     =  Name_2d_string.c_str();
    Name_2dHistMap[fShowerModuleLabels[j]]  = tfs->make<TH2F>(Name_2d_name, Name_2d_titlename,x_numbins,x_start,x_end,y_numbins,y_start,y_end);

    //Create the 2D Histogram of Energy and value 
    std::string Name_2d_canvas_string      = Name + "_2d_canvas_" + fShowerModuleLabels[j];
    std::string Name_2d_canvas_titlestring = ";" + TitleName + ";Energy (MeV);Enteries";
    const char* Name_2d_canvas_titlename   = Name_2d_canvas_titlestring.c_str();
    const char* Name_2d_canvas_name     =  Name_2d_canvas_string.c_str();
    Name_2dCanvasMap[fShowerModuleLabels[j]]  = tfs->makeAndRegister<TCanvas>(Name_2d_canvas_name, Name_2d_canvas_titlename);

    for(unsigned int i=0; i<fEnergies.size(); ++i){
      
      //Create the Histogram maps name and string for the different energies. 
      std::string Energies_Name_string = Name + fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)";  
      std::string Energies_Name_titlestring  =  fShowerModuleLabels[j] +  ";" + TitleName + ";Enteries";
      const char* Energies_Name_titlename    = Energies_Name_titlestring.c_str();
      const char* Energies_Name_name         = Energies_Name_string.c_str();
      Energies_Name_HistMap[fShowerModuleLabels[j]][fEnergies[i]]   = tfs->make<TH1F>(Energies_Name_name, Energies_Name_titlename, x_numbins, x_start, x_end);
    }
  

    //Create the Independent TGraphs for the Metric vs Energy graphs    
    std::string meangrapherror_string      = Name + fShowerModuleLabels[j] + "_meangrapherror";
    std::string meangrapherror_titlestring = fShowerModuleLabels[j];//+ ";Energy (Mev);" + TitleName;
    const char* meangrapherror_name        = meangrapherror_string.c_str();
    const char* meangrapherror_titlename   = meangrapherror_titlestring.c_str();
    Energies_Mean_Name_GraphMap[fShowerModuleLabels[j]] = tfs->makeAndRegister<TGraphErrors>(meangrapherror_name,meangrapherror_titlename);
    
    //Create the Independent TGraphs for the Metric vs Energy graphs    
    std::string rmsgrapherror_string      = Name + fShowerModuleLabels[j] + "_rmsgrapherror";
    std::string rmsgrapherror_titlestring =  fShowerModuleLabels[j];//+ ";Energy (Mev);" + TitleName;
    const char* rmsgrapherror_name        = rmsgrapherror_string.c_str();
    const char* rmsgrapherror_titlename   = rmsgrapherror_titlestring.c_str();
    Energies_RMS_Name_GraphMap[fShowerModuleLabels[j]] = tfs->makeAndRegister<TGraphErrors>(rmsgrapherror_name,rmsgrapherror_titlename);
  
  }  
  
  //Canvas for all energy histogram
  std::string canvas_string      = Name + "_canvas";
  std::string canvas_titlestring = ";" + TitleName + ";Enteries";
  const char* canvas_name        = canvas_string.c_str();
  const char* canvas_titlename   = canvas_titlestring.c_str();
  Name_canvas                    = tfs->makeAndRegister<TCanvas>(canvas_name,canvas_titlename);
  
  //MultiGraph for Energy vs Mean metric
  std::string meanmulti_string      = Name + "_meanmulti";
  std::string meanmulti_titlestring = ";Energy (Mev);" + TitleName;
  const char* meanmulti_name        = meanmulti_string.c_str();
  const char* meanmulti_titlename   = meanmulti_titlestring.c_str();
  Energies_Mean_Name_Multi          = tfs->makeAndRegister<TMultiGraph>(meanmulti_name,meanmulti_titlename);
  
  //MultiGraph for Energy vs RMS metric
  std::string rmsmulti_string      = Name + "_rmsmulti";
  std::string rmsmulti_titlestring = ";Energy (Mev);" + TitleName;
  const char* rmsmulti_name        = rmsmulti_string.c_str();
  const char* rmsmulti_titlename   = rmsmulti_titlestring.c_str();
  Energies_RMS_Name_Multi          = tfs->makeAndRegister<TMultiGraph>(rmsmulti_name,rmsmulti_titlename);
  
  //Canvas for MultiGraph for Energy vs Mean metric
  std::string meancanvasmulti_string      = Name + "_meancanvasmulti";
  std::string meancanvasmulti_titlestring = ";Energy (Mev);" + TitleName;
  const char* meancanvasmulti_name        = meancanvasmulti_string.c_str();
  const char* meancanvasmulti_titlename   = meancanvasmulti_titlestring.c_str();
  Energies_Mean_Name_canvasMulti          = tfs->makeAndRegister<TCanvas>(meancanvasmulti_name,meancanvasmulti_titlename);
  
  //Canvasmulti for MultiGraph for Energy vs RMS metric
  std::string rmscanvasmulti_string      = Name + "_rmscanvasmulti";
  std::string rmscanvasmulti_titlestring = ";Energy (Mev);" + TitleName;
  const char* rmscanvasmulti_name        = rmscanvasmulti_string.c_str();
  const char* rmscanvasmulti_titlename   = rmscanvasmulti_titlestring.c_str();
  Energies_RMS_Name_canvasMulti          = tfs->makeAndRegister<TCanvas>(rmscanvasmulti_name,rmscanvasmulti_titlename);

}

void ana::ShowerValidation::InitaliseGraphs(std::string Name, std::string TitleName,
					    std::map<std::string,std::map<geo::PlaneID,TH1F*> >& Name_HistMap,
					    std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > >& Energies_Name_HistMap,
					    std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> >& Energies_Mean_Name_GraphMap,
					    std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> >& Energies_RMS_Name_GraphMap,
					    std::map<geo::PlaneID,TMultiGraph*>& Energies_Mean_Name_Multi,
					    std::map<geo::PlaneID,TMultiGraph*>& Energies_RMS_Name_Multi,
					    std::map<geo::PlaneID,TCanvas*>& Energies_Mean_Name_canvasMulti,
					    std::map<geo::PlaneID,TCanvas*>& Energies_RMS_Name_canvasMulti,
					    std::map<float,std::map<geo::PlaneID,TCanvas*> >& Energies_Name_canvasMap,
					    std::map<geo::PlaneID,TCanvas*>& Name_canvas,
					    std::map<std::string,std::map<geo::PlaneID,TH2F*> >& Name_2dHistMap,
					    std::map<std::string,std::map<geo::PlaneID,TCanvas*> >& Name_2dCanvasMap,
					    int x_numbins, float x_start, float x_end,
					    int y_numbins, float y_start, float y_end,
					    std::map<std::string,std::vector<std::vector<float > > >& MetricVector
					    ){
  
  gStyle->SetOptStat(0);

  for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){

    std::string Plane_string = " Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat); 

    //Create The histogram canvas for Each Energy. 
    for(unsigned int i=0; i<fEnergies.size(); ++i){
      std::string Energies_string = Name + std::to_string(fEnergies[i]) + Plane_string;
      std::string Energies_titlestring = ";" + TitleName + ";Enteries";
      const char* Energies_name   = Energies_string.c_str();
      const char* Energies_titlename   = Energies_titlestring.c_str();;  
      Energies_Name_canvasMap[fEnergies[i]][plane_id] = tfs->makeAndRegister<TCanvas>(Energies_name,Energies_titlename);
    }

    for(unsigned int j=0; j<fShowerModuleLabels.size(); ++j){
      
      //Create the Tree Branch
      std::string Tree_string = Name + "_" +fShowerModuleLabels[j];
      const char* Tree_name   = Tree_string.c_str();
      Tree->Branch(Tree_name,"std::vector<std::vector<float> >", &MetricVector[fShowerModuleLabels[j]],256000,0);
      
      //Create the Hitogram which has all the enrgies in
      std::string Name_string      = Name + fShowerModuleLabels[j] +  Plane_string;
      std::string Name_titlestring = fShowerModuleLabels[j] + ";" + TitleName + ";Enteries";
      const char* Name_titlename   = Name_titlestring.c_str();
      const char* Name_name     =  Name_string.c_str();
      Name_HistMap[fShowerModuleLabels[j]][plane_id]  = tfs->make<TH1F>(Name_name, Name_titlename,x_numbins,x_start,x_end);

      //Create the 2D Histogram of Energy and value 
      std::string Name_2d_string      = Name + "_2d" + fShowerModuleLabels[j] + Plane_string;
      std::string Name_2d_titlestring = ";" + TitleName + ";Energy (MeV);Enteries";
      const char* Name_2d_titlename   = Name_2d_titlestring.c_str();
      const char* Name_2d_name     =  Name_2d_string.c_str();
      Name_2dHistMap[fShowerModuleLabels[j]][plane_id]  = tfs->make<TH2F>(Name_2d_name, Name_2d_titlename,x_numbins,x_start,x_end,y_numbins,y_start,y_end);

      //Create the 2D Histogram of Energy and value 
      std::string Name_2d_canvas_string      = Name + "_2d_canvas_" + fShowerModuleLabels[j] + Plane_string;
      std::string Name_2d_canvas_titlestring = ";" + TitleName + ";Energy (MeV);Enteries";
      const char* Name_2d_canvas_titlename   = Name_2d_canvas_titlestring.c_str();
      const char* Name_2d_canvas_name     =  Name_2d_canvas_string.c_str();
      Name_2dCanvasMap[fShowerModuleLabels[j]][plane_id]  = tfs->makeAndRegister<TCanvas>(Name_2d_canvas_name, Name_2d_canvas_titlename);
      
      for(unsigned int i=0; i<fEnergies.size(); ++i){
	
	//Create the Histogram maps name and string for the different energies. 
	std::string Energies_Name_string = Name + fShowerModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)"  + " at Energy" + std::to_string(fEnergies[i]) +  Plane_string;
	std::string Energies_Name_titlestring  = fShowerModuleLabels[j] +  ";" + TitleName + ";Enteries";
	const char* Energies_Name_titlename    = Energies_Name_titlestring.c_str();
	const char* Energies_Name_name         = Energies_Name_string.c_str();
	Energies_Name_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]  = tfs->make<TH1F>(Energies_Name_name, Energies_Name_titlename, x_numbins,x_start,x_end);
      }
    
      
      //Create the Independent TGraphs for the Metric vs Energy graphs    
      std::string meangrapherror_string      = Name + fShowerModuleLabels[j] + Plane_string + "_meangrapherror";
      std::string meangrapherror_titlestring = fShowerModuleLabels[j];// + ";Energy (Mev);" + TitleName;
      const char* meangrapherror_name        = meangrapherror_string.c_str();
      const char* meangrapherror_titlename   = meangrapherror_titlestring.c_str();
      Energies_Mean_Name_GraphMap[fShowerModuleLabels[j]][plane_id] = tfs->makeAndRegister<TGraphErrors>(meangrapherror_name,meangrapherror_titlename);
      
      //Create the Independent TGraphs for the Metric vs Energy graphs    
      std::string rmsgrapherror_string      = Name +  fShowerModuleLabels[j] + Plane_string +  "_rmsgrapherror";
      std::string rmsgrapherror_titlestring = fShowerModuleLabels[j];// +  ";Energy (Mev);" + TitleName;
      const char* rmsgrapherror_name        = rmsgrapherror_string.c_str();
      const char* rmsgrapherror_titlename   = rmsgrapherror_titlestring.c_str();
      Energies_RMS_Name_GraphMap[fShowerModuleLabels[j]][plane_id] = tfs->makeAndRegister<TGraphErrors>(rmsgrapherror_name,rmsgrapherror_titlename);
    }
    
    //Canvas for all energy histogram
    std::string canvas_string      = Name + Plane_string + "_canvas";
    std::string canvas_titlestring = ";Energy (Mev);" + TitleName;
    const char* canvas_name        = canvas_string.c_str();
    const char* canvas_titlename   = canvas_titlestring.c_str();
    Name_canvas[plane_id]          = tfs->makeAndRegister<TCanvas>(canvas_name,canvas_titlename);
    
    //MultiGraph for Energy vs Mean metric
    std::string meanmulti_string       = Name + Plane_string + "_meanmulti";
    std::string meanmulti_titlestring  = ";Energy (Mev);" + TitleName;
    const char* meanmulti_name         = meanmulti_string.c_str();
    const char* meanmulti_titlename    = meanmulti_titlestring.c_str();
    Energies_Mean_Name_Multi[plane_id] = tfs->makeAndRegister<TMultiGraph>(meanmulti_name,meanmulti_titlename);
    
    //MultiGraph for Energy vs RMS metric
    std::string rmsmulti_string       = Name + Plane_string + "_rmsmulti";
    std::string rmsmulti_titlestring  = ";Energy (Mev);" + TitleName;
    const char* rmsmulti_name         = rmsmulti_string.c_str();
    const char* rmsmulti_titlename    = rmsmulti_titlestring.c_str();
    Energies_RMS_Name_Multi[plane_id] = tfs->makeAndRegister<TMultiGraph>(rmsmulti_name,rmsmulti_titlename);
    
    //Canvas for MultiGraph for Energy vs Mean metric
    std::string meancanvasmulti_string       = Name + Plane_string + "_meancanvasmulti";
    std::string meancanvasmulti_titlestring  = ";Energy (Mev);" + TitleName;
    const char* meancanvasmulti_name         = meancanvasmulti_string.c_str();
    const char* meancanvasmulti_titlename    = meancanvasmulti_titlestring.c_str();
    Energies_Mean_Name_canvasMulti[plane_id] = tfs->makeAndRegister<TCanvas>(meancanvasmulti_name,meancanvasmulti_titlename);
    
    //Canvasmulti for MultiGraph for Energy vs RMS metric
    std::string rmscanvasmulti_string       = Name + Plane_string + "_rmscanvasmulti";
    std::string rmscanvasmulti_titlestring  = ";Energy (Mev);" + TitleName;
    const char* rmscanvasmulti_name         = rmscanvasmulti_string.c_str();
    const char* rmscanvasmulti_titlename    = rmscanvasmulti_titlestring.c_str();
    Energies_RMS_Name_canvasMulti[plane_id] = tfs->makeAndRegister<TCanvas>(rmscanvasmulti_name,rmscanvasmulti_titlename);
  }
}

//For Hit Valiation
void ana::ShowerValidation::InitaliseHitGraphs(std::string Name, std::string TitleName,
					       std::map<std::string,std::map<geo::PlaneID,TH1F*> >& Name_HistMap,
					       std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > >& Energies_Name_HistMap,
					       std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> >& Energies_Mean_Name_GraphMap,
					       std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> >& Energies_RMS_Name_GraphMap,
					       std::map<geo::PlaneID,TMultiGraph*>& Energies_Mean_Name_Multi,
					       std::map<geo::PlaneID,TMultiGraph*>& Energies_RMS_Name_Multi,
					       std::map<geo::PlaneID,TCanvas*>& Energies_Mean_Name_canvasMulti,
					       std::map<geo::PlaneID,TCanvas*>& Energies_RMS_Name_canvasMulti,
					       std::map<float,std::map<geo::PlaneID,TCanvas*> >& Energies_Name_canvasMap,
					       std::map<geo::PlaneID,TCanvas*>& Name_canvas,
					       std::map<std::string,std::map<geo::PlaneID,TH2F*> >& Name_2dHistMap,
					       std::map<std::string,std::map<geo::PlaneID,TCanvas*> >& Name_2dCanvasMap,
					       int x_numbins, float x_start, float x_end,
					       int y_numbins, float y_start, float y_end,
					       std::map<std::string,std::vector<std::vector<float> > >& MetricVector
					       ){
  
  gStyle->SetOptStat(0);

  //Initialise Trees
  for(unsigned int j=0; j<fHitModuleLabels.size(); ++j){
    std::cout << geom->Nplanes() << std::endl;
    MetricVector[fHitModuleLabels[j]].resize(geom->Nplanes());
  }

  for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){

    std::string Plane_string = " Plane ID: P" + std::to_string(plane_id.Plane) + " T" + std::to_string(plane_id.TPC) + " C" + std::to_string(plane_id.Cryostat); 

    //Create The histogram canvas for Each Energy. 
    for(unsigned int i=0; i<fEnergies.size(); ++i){
      std::string Energies_string = Name + std::to_string(fEnergies[i]) + Plane_string;
      std::string Energies_titlestring = ";" + TitleName + ";Enteries";
      const char* Energies_name   = Energies_string.c_str();
      const char* Energies_titlename   = Energies_titlestring.c_str();;  
      Energies_Name_canvasMap[fEnergies[i]][plane_id] = tfs->makeAndRegister<TCanvas>(Energies_name,Energies_titlename);
    }

    for(unsigned int j=0; j<fHitModuleLabels.size(); ++j){

      //Create the Tree Branch 
      std::string Tree_string = Name + "_" + fHitModuleLabels[j];
      const char* Tree_name   = Tree_string.c_str();
      Tree->Branch(Tree_name,"std::vector<std::vector<float> >", &MetricVector[fHitModuleLabels[j]], 32000, 0);
      
      //Create the Hitogram which has all the enrgies in
      std::string Name_string      = Name + fHitModuleLabels[j] +  Plane_string;
      std::string Name_titlestring = fHitModuleLabels[j] + ";" + TitleName + ";Enteries";
      const char* Name_titlename   = Name_titlestring.c_str();
      const char* Name_name     =  Name_string.c_str();
      Name_HistMap[fHitModuleLabels[j]][plane_id]  = tfs->make<TH1F>(Name_name, Name_titlename,x_numbins,x_start,x_end);
      
      //Create the 2D Histogram of Energy and value 
      std::string Name_2d_string      = Name + "_2d_" + fHitModuleLabels[j] + Plane_string;
      std::string Name_2d_titlestring = ";" + TitleName + ";Energy (MeV);Enteries";
      const char* Name_2d_titlename   = Name_2d_titlestring.c_str();
      const char* Name_2d_name     =  Name_2d_string.c_str();
      Name_2dHistMap[fHitModuleLabels[j]][plane_id]  = tfs->make<TH2F>(Name_2d_name, Name_2d_titlename,x_numbins,x_start,x_end,y_numbins,y_start,y_end);
      
      //Create the 2D Histogram of Energy and value 
      std::string Name_2d_canvas_string      = Name + "_2d_canvas_" + fHitModuleLabels[j] + Plane_string;
      std::string Name_2d_canvas_titlestring = ";" + TitleName + ";Energy;Enteries";
      const char* Name_2d_canvas_titlename   = Name_2d_canvas_titlestring.c_str();
      const char* Name_2d_canvas_name     =  Name_2d_canvas_string.c_str();
      Name_2dCanvasMap[fHitModuleLabels[j]][plane_id]  = tfs->makeAndRegister<TCanvas>(Name_2d_canvas_name, Name_2d_canvas_titlename);

      for(unsigned int i=0; i<fEnergies.size(); ++i){
	
	//Create the Histogram maps name and string for the different energies. 
	std::string Energies_Name_string = Name + fHitModuleLabels[j] + " at Energy" + std::to_string(fEnergies[i]) + " (MeV)"  + " at Energy" + std::to_string(fEnergies[i]) +  Plane_string;
	std::string Energies_Name_titlestring  = fHitModuleLabels[j] + ";" + TitleName + ";Enteries";
	const char* Energies_Name_titlename    = Energies_Name_titlestring.c_str();
	const char* Energies_Name_name         = Energies_Name_string.c_str();
	Energies_Name_HistMap[fHitModuleLabels[j]][fEnergies[i]][plane_id]  = tfs->make<TH1F>(Energies_Name_name, Energies_Name_titlename, x_numbins,x_start,x_end);
      }
    
      
      //Create the Independent TGraphs for the Metric vs Energy graphs    
      std::string meangrapherror_string      = Name + fHitModuleLabels[j] + Plane_string + "_meangrapherror";
      std::string meangrapherror_titlestring = fShowerModuleLabels[j];// + ";Energy (Mev);" + TitleName;
      const char* meangrapherror_name        = meangrapherror_string.c_str();
      const char* meangrapherror_titlename   = meangrapherror_titlestring.c_str();
      Energies_Mean_Name_GraphMap[fHitModuleLabels[j]][plane_id] = tfs->makeAndRegister<TGraphErrors>(meangrapherror_name,meangrapherror_titlename);
      
      //Create the Independent TGraphs for the Metric vs Energy graphs    
      std::string rmsgrapherror_string      = Name +  fHitModuleLabels[j] + Plane_string +  "_rmsgrapherror";
      std::string rmsgrapherror_titlestring = fShowerModuleLabels[j];// + ";Energy (Mev);" + TitleName;
      const char* rmsgrapherror_name        = rmsgrapherror_string.c_str();
      const char* rmsgrapherror_titlename   = rmsgrapherror_titlestring.c_str();
      Energies_RMS_Name_GraphMap[fHitModuleLabels[j]][plane_id] = tfs->makeAndRegister<TGraphErrors>(rmsgrapherror_name,rmsgrapherror_titlename);
    }
    
    //Canvas for all energy histogram
    std::string canvas_string      = Name + Plane_string + "_canvas";
    std::string canvas_titlestring = ";Energy (Mev);" + TitleName;
    const char* canvas_name        = canvas_string.c_str();
    const char* canvas_titlename   = canvas_titlestring.c_str();
    Name_canvas[plane_id]          = tfs->makeAndRegister<TCanvas>(canvas_name,canvas_titlename);
    
    //MultiGraph for Energy vs Mean metric
    std::string meanmulti_string       = Name + Plane_string + "_meanmulti";
    std::string meanmulti_titlestring  = ";Energy (Mev);" + TitleName;
    const char* meanmulti_name         = meanmulti_string.c_str();
    const char* meanmulti_titlename    = meanmulti_titlestring.c_str();
    Energies_Mean_Name_Multi[plane_id] = tfs->makeAndRegister<TMultiGraph>(meanmulti_name,meanmulti_titlename);
    
    //MultiGraph for Energy vs RMS metric
    std::string rmsmulti_string       = Name + Plane_string + "_rmsmulti";
    std::string rmsmulti_titlestring  = ";Energy (Mev);" + TitleName;
    const char* rmsmulti_name         = rmsmulti_string.c_str();
    const char* rmsmulti_titlename    = rmsmulti_titlestring.c_str();
    Energies_RMS_Name_Multi[plane_id] = tfs->makeAndRegister<TMultiGraph>(rmsmulti_name,rmsmulti_titlename);
    
    //Canvas for MultiGraph for Energy vs Mean metric
    std::string meancanvasmulti_string       = Name + Plane_string + "_meancanvasmulti";
    std::string meancanvasmulti_titlestring  = ";Energy (Mev);" + TitleName;
    const char* meancanvasmulti_name         = meancanvasmulti_string.c_str();
    const char* meancanvasmulti_titlename    = meancanvasmulti_titlestring.c_str();
    Energies_Mean_Name_canvasMulti[plane_id] = tfs->makeAndRegister<TCanvas>(meancanvasmulti_name,meancanvasmulti_titlename);
    
    //Canvasmulti for MultiGraph for Energy vs RMS metric
    std::string rmscanvasmulti_string       = Name + Plane_string + "_rmscanvasmulti";
    std::string rmscanvasmulti_titlestring  = ";Energy (Mev);" + TitleName;
    const char* rmscanvasmulti_name         = rmscanvasmulti_string.c_str();
    const char* rmscanvasmulti_titlename    = rmscanvasmulti_titlestring.c_str();
    Energies_RMS_Name_canvasMulti[plane_id] = tfs->makeAndRegister<TCanvas>(rmscanvasmulti_name,rmscanvasmulti_titlename);
  }
}

  
void ana::ShowerValidation::beginJob() {
  
  Tree = tfs->make<TTree>("MetricTree", "Tree Holding all metric information");

  for(unsigned int j=0; j<fShowerModuleLabels.size(); ++j){

    //Initialise the Pos*Dir 
    std::string PosDir_string       = "Position and direction 2D Hist" + fShowerModuleLabels[j];
    std::string PosDir_titlestring  = ";Position;Direction;";
    const char* PosDir_name         = PosDir_string.c_str();
    const char* PosDir_titlename    = PosDir_titlestring.c_str();
    PosDir_2dHistMap[fShowerModuleLabels[j]] = tfs->make<TH2F>(PosDir_name, PosDir_titlename,200,0,50,100,-1,1);
    PosDir_2dCanvasMap[fShowerModuleLabels[j]] = tfs->makeAndRegister<TCanvas>(PosDir_name,PosDir_titlename);
    
  }

  BiggestShowerMotherE_AfterCutHist  = tfs->make<TH1F>("BiggestShowerMotherE_AfterCutHist", ";Mother Energy (MeV);Enteries",100,0,fMaxSimEnergy);
  BiggestShowerMotherE_BeforeCutHist = tfs->make<TH1F>("BiggestShowerMotherE_BeforeCutHist", ";Mother Energy (MeV);Enteries",100,0,fMaxSimEnergy);
  AsscoiatedBiggestShowerMotherEHist = tfs->make<TH1F>("AsscoiatedBiggestShowerMotherEHist", ";Mother Energy (MeV);Enteries",100,0,fMaxSimEnergy);

  BiggestShowerMotherE_AfterCutCanvas  = tfs->makeAndRegister<TCanvas>("BiggestShowerMotherE_AfterCutCanvas", ";Mother Energy (MeV);Enteries");
  BiggestShowerMotherE_BeforeCutCanvas = tfs->makeAndRegister<TCanvas>("BiggestShowerMotherE_BeforeCutCanvas", ";Mother Energy (MeV);Enteries");
  AsscoiatedBiggestShowerMotherECanvas = tfs->makeAndRegister<TCanvas>("AsscoiatedBiggestShowerMotherECanvas", ";Mother Energy (MeV);Enteries");

  SmallestShowerMotherE_AfterCutHist  = tfs->make<TH1F>("SmallestShowerMotherE_AfterCutHist", ";Mother Energy (MeV);Enteries",100,0,fMaxSimEnergy);
  SmallestShowerMotherE_BeforeCutHist = tfs->make<TH1F>("SmallestShowerMotherE_BeforeCutHist", ";Mother Energy (MeV);Enteries",100,0,fMaxSimEnergy);
  AsscoiatedSmallestShowerMotherEHist = tfs->make<TH1F>("AsscoiatedSmallestShowerMotherEHist", ";Mother Energy (MeV);Enteries",100,0,fMaxSimEnergy);

  SmallestShowerMotherE_AfterCutCanvas  = tfs->makeAndRegister<TCanvas>("SmallestShowerMotherE_AfterCutCanvas", ";Mother Energy (MeV);Enteries");
  SmallestShowerMotherE_BeforeCutCanvas = tfs->makeAndRegister<TCanvas>("SmallestShowerMotherE_AfterCutCanvas", ";Mother Energy (MeV);Enteries");
  AsscoiatedSmallestShowerMotherECanvas = tfs->makeAndRegister<TCanvas>("SmallestShowerMotherE_AfterCutCanvas", ";Mother Energy (MeV);Enteries");

  Tree->Branch("SmallestShowerMotherE_BeforeCut","std::vector<float>", &SmallestShowerMotherE_BeforeCut_TreeVal, 32000, 0);
  Tree->Branch("BiggestShowerMotherE_BeforeCut","std::vector<float>", &BiggestShowerMotherE_BeforeCut_TreeVal, 32000, 0);
  Tree->Branch("SmallestShowerMotherE_AfterCut","std::vector<float>", &SmallestShowerMotherE_AfterCut_TreeVal, 32000, 0);
  Tree->Branch("BiggestShowerMotherE_AfterCut","std::vector<float>", &BiggestShowerMotherE_AfterCut_TreeVal, 32000, 0);
  
  Tree->Branch("EventNumber","std::vectorstd::float>",&EventNum_name_TreeVal,32000,0);


  for(unsigned int j=0; j<fShowerModuleLabels.size(); ++j){
    std::string Process_string = "ShowerStartEndProcess_" + fShowerModuleLabels[j];
    const char* Process_name   = Process_string.c_str();
    Tree->Branch(Process_name,"std::vector<std::string>",&ShowerStartEndProcess_TreeVal[fShowerModuleLabels[j]],32000,0);
  }

  ana::ShowerValidation::InitaliseGraphs("ShowerDirection_X","cos(x)",ShowerDirection_X_HistMap,Energies_ShowerDirection_X_HistMap,Energies_Mean_ShowerDirection_X_GraphMap,Energies_RMS_ShowerDirection_X_GraphMap,Energies_Mean_ShowerDirection_X_Multi,Energies_RMS_ShowerDirection_X_Multi,Energies_Mean_ShowerDirection_X_canvasMulti,Energies_RMS_ShowerDirection_X_canvasMulti,Energies_ShowerDirection_X_canvasMap, ShowerDirection_X_canvas,ShowerDirection_X_2dHistMap,ShowerDirection_X_2dCanvasMap,100,-1,1,200,0,fMaxSimEnergy,ShowerDirection_X_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ShowerDirection_Y","cos(y)",ShowerDirection_Y_HistMap,Energies_ShowerDirection_Y_HistMap,Energies_Mean_ShowerDirection_Y_GraphMap,Energies_RMS_ShowerDirection_Y_GraphMap,Energies_Mean_ShowerDirection_Y_Multi,Energies_RMS_ShowerDirection_Y_Multi,Energies_Mean_ShowerDirection_Y_canvasMulti,Energies_RMS_ShowerDirection_Y_canvasMulti,Energies_ShowerDirection_Y_canvasMap, ShowerDirection_Y_canvas,ShowerDirection_Y_2dHistMap,ShowerDirection_Y_2dCanvasMap,100,-1,1,200,0,fMaxSimEnergy,ShowerDirection_Y_TreeVal);
  
  ana::ShowerValidation::InitaliseGraphs("ShowerDirection_Z","cos(z)",ShowerDirection_Z_HistMap,Energies_ShowerDirection_Z_HistMap,Energies_Mean_ShowerDirection_Z_GraphMap,Energies_RMS_ShowerDirection_Z_GraphMap,Energies_Mean_ShowerDirection_Z_Multi,Energies_RMS_ShowerDirection_Z_Multi,Energies_Mean_ShowerDirection_Z_canvasMulti,Energies_RMS_ShowerDirection_Z_canvasMulti,Energies_ShowerDirection_Z_canvasMap, ShowerDirection_Z_canvas,ShowerDirection_Z_2dHistMap,ShowerDirection_Z_2dCanvasMap,100,-1,1,200,0,fMaxSimEnergy,ShowerDirection_Z_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ShowerStart_X","|True Shower Start X - Reco Shower Start X|",ShowerStart_X_HistMap,Energies_ShowerStart_X_HistMap,Energies_Mean_ShowerStart_X_GraphMap,Energies_RMS_ShowerStart_X_GraphMap,Energies_Mean_ShowerStart_X_Multi,Energies_RMS_ShowerStart_X_Multi,Energies_Mean_ShowerStart_X_canvasMulti,Energies_RMS_ShowerStart_X_canvasMulti,Energies_ShowerStart_X_canvasMap, ShowerStart_X_canvas,ShowerStart_X_2dHistMap,ShowerStart_X_2dCanvasMap,200,0,50,100,0,fMaxSimEnergy,ShowerStart_X_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ShowerStart_Y","|True Shower Start Y - Reco Shower Start Y|",ShowerStart_Y_HistMap,Energies_ShowerStart_Y_HistMap,Energies_Mean_ShowerStart_Y_GraphMap,Energies_RMS_ShowerStart_Y_GraphMap,Energies_Mean_ShowerStart_Y_Multi,Energies_RMS_ShowerStart_Y_Multi,Energies_Mean_ShowerStart_Y_canvasMulti,Energies_RMS_ShowerStart_Y_canvasMulti,Energies_ShowerStart_Y_canvasMap, ShowerStart_Y_canvas,ShowerStart_Y_2dHistMap,ShowerStart_Y_2dCanvasMap,200,0,50,100,0,fMaxSimEnergy,ShowerStart_Y_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ShowerStart_Z","|True Shower Start Z - Reco Shower Start Z|",ShowerStart_Z_HistMap,Energies_ShowerStart_Z_HistMap,Energies_Mean_ShowerStart_Z_GraphMap,Energies_RMS_ShowerStart_Z_GraphMap,Energies_Mean_ShowerStart_Z_Multi,Energies_RMS_ShowerStart_Z_Multi,Energies_Mean_ShowerStart_Z_canvasMulti,Energies_RMS_ShowerStart_Z_canvasMulti,Energies_ShowerStart_Z_canvasMap, ShowerStart_Z_canvas,ShowerStart_Z_2dHistMap,ShowerStart_Z_2dCanvasMap,200,0,50,100,0,fMaxSimEnergy,ShowerStart_Z_TreeVal);
   
  ana::ShowerValidation::InitaliseGraphs("ShowerLength","Shower Length (cm)",ShowerLength_HistMap,Energies_ShowerLength_HistMap,Energies_Mean_ShowerLength_GraphMap,Energies_RMS_ShowerLength_GraphMap,Energies_Mean_ShowerLength_Multi,Energies_RMS_ShowerLength_Multi,Energies_Mean_ShowerLength_canvasMulti,Energies_RMS_ShowerLength_canvasMulti,Energies_ShowerLength_canvasMap,ShowerLength_canvas,ShowerLength_2dHistMap,ShowerLength_2dCanvasMap,100,0,500,100,0,fMaxSimEnergy,ShowerLength_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ShowerEnergyDiff","Reco Energy/True Energy",ShowerEnergyDiff_HistMap,Energies_ShowerEnergyDiff_HistMap,Energies_Mean_ShowerEnergyDiff_GraphMap,Energies_RMS_ShowerEnergyDiff_GraphMap,Energies_Mean_ShowerEnergyDiff_Multi,Energies_RMS_ShowerEnergyDiff_Multi,Energies_Mean_ShowerEnergyDiff_canvasMulti,Energies_RMS_ShowerEnergyDiff_canvasMulti,Energies_ShowerEnergyDiff_canvasMap, ShowerEnergyDiff_canvas,ShowerEnergyDiff_2dHistMap,ShowerEnergyDiff_2dCanvasMap,100,0,1.2,100,0,fMaxSimEnergy,ShowerEnergyDiff_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ShowerdEdx","dEdx (meV)",ShowerdEdx_HistMap,Energies_ShowerdEdx_HistMap,Energies_Mean_ShowerdEdx_GraphMap,Energies_RMS_ShowerdEdx_GraphMap,Energies_Mean_ShowerdEdx_Multi,Energies_RMS_ShowerdEdx_Multi,Energies_Mean_ShowerdEdx_canvasMulti,Energies_RMS_ShowerdEdx_canvasMulti,Energies_ShowerdEdx_canvasMap, ShowerdEdx_canvas,ShowerdEdx_2dHistMap,ShowerdEdx_2dCanvasMap,100,0,10,100,0,fMaxSimEnergy,ShowerdEdx_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("EventSeggy","Shower Segmenation",EventSeggy_HistMap,Energies_EventSeggy_HistMap,Energies_Mean_EventSeggy_GraphMap,Energies_RMS_EventSeggy_GraphMap,Energies_Mean_EventSeggy_Multi,Energies_RMS_EventSeggy_Multi,Energies_Mean_EventSeggy_canvasMulti,Energies_RMS_EventSeggy_canvasMulti,Energies_EventSeggy_canvasMap, EventSeggy_canvas,EventSeggy_2dHistMap,EventSeggy_2dCanvasMap,20,-0.5,20.5,100,0,fMaxSimEnergy,EventSeggy_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ShowerEnergyCompleteness"," Energy Completeness",ShowerEnergyCompleteness_HistMap,Energies_ShowerEnergyCompleteness_HistMap,Energies_Mean_ShowerEnergyCompleteness_GraphMap,Energies_RMS_ShowerEnergyCompleteness_GraphMap,Energies_Mean_ShowerEnergyCompleteness_Multi,Energies_RMS_ShowerEnergyCompleteness_Multi,Energies_Mean_ShowerEnergyCompleteness_canvasMulti,Energies_RMS_ShowerEnergyCompleteness_canvasMulti,Energies_ShowerEnergyCompleteness_canvasMap, ShowerEnergyCompleteness_canvas,ShowerEnergyCompleteness_2dHistMap,ShowerEnergyCompleteness_2dCanvasMap,100,-0.2,1.2,100,0,fMaxSimEnergy,ShowerEnergyCompleteness_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ShowerEnergyPurity","Energy Purity",ShowerEnergyPurity_HistMap,Energies_ShowerEnergyPurity_HistMap,Energies_Mean_ShowerEnergyPurity_GraphMap,Energies_RMS_ShowerEnergyPurity_GraphMap,Energies_Mean_ShowerEnergyPurity_Multi,Energies_RMS_ShowerEnergyPurity_Multi,Energies_Mean_ShowerEnergyPurity_canvasMulti,Energies_RMS_ShowerEnergyPurity_canvasMulti,Energies_ShowerEnergyPurity_canvasMap, ShowerEnergyPurity_canvas,ShowerEnergyPurity_2dHistMap,ShowerEnergyPurity_2dCanvasMap,100,0,1.2,100,0,fMaxSimEnergy,ShowerEnergyPurity_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ShowerHitsCompleteness","Hit Completeness",ShowerHitsCompleteness_HistMap,Energies_ShowerHitsCompleteness_HistMap,Energies_Mean_ShowerHitsCompleteness_GraphMap,Energies_RMS_ShowerHitsCompleteness_GraphMap,Energies_Mean_ShowerHitsCompleteness_Multi,Energies_RMS_ShowerHitsCompleteness_Multi,Energies_Mean_ShowerHitsCompleteness_canvasMulti,Energies_RMS_ShowerHitsCompleteness_canvasMulti,Energies_ShowerHitsCompleteness_canvasMap, ShowerHitsCompleteness_canvas,ShowerHitsCompleteness_2dHistMap,ShowerHitsCompleteness_2dCanvasMap,100,-0.2,1.2,100,0,fMaxSimEnergy,ShowerHitsCompleteness_TreeVal);
  
  ana::ShowerValidation::InitaliseGraphs("ShowerHitsPurity","Hit Purity",ShowerHitsPurity_HistMap,Energies_ShowerHitsPurity_HistMap,Energies_Mean_ShowerHitsPurity_GraphMap,Energies_RMS_ShowerHitsPurity_GraphMap,Energies_Mean_ShowerHitsPurity_Multi,Energies_RMS_ShowerHitsPurity_Multi,Energies_Mean_ShowerHitsPurity_canvasMulti,Energies_RMS_ShowerHitsPurity_canvasMulti,Energies_ShowerHitsPurity_canvasMap, ShowerHitsPurity_canvas,ShowerHitsPurity_2dHistMap,ShowerHitsPurity_2dCanvasMap,100,0,1.2,100,0,fMaxSimEnergy,ShowerHitsPurity_TreeVal);
  
  ana::ShowerValidation::InitaliseGraphs("ShowerEnergy","Shower Energy (MeV)",ShowerEnergy_HistMap,Energies_ShowerEnergy_HistMap,Energies_Mean_ShowerEnergy_GraphMap,Energies_RMS_ShowerEnergy_GraphMap,Energies_Mean_ShowerEnergy_Multi,Energies_RMS_ShowerEnergy_Multi,Energies_Mean_ShowerEnergy_canvasMulti,Energies_RMS_ShowerEnergy_canvasMulti,Energies_ShowerEnergy_canvasMap,ShowerEnergy_canvas,ShowerEnergy_2dHistMap,ShowerEnergy_2dCanvasMap,100,0,fMaxSimEnergy,100,0,fMaxSimEnergy,ShowerEnergy_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ShowerHitNum","Hit Number in reco showers",ShowerHitNum_HistMap,Energies_ShowerHitNum_HistMap,Energies_Mean_ShowerHitNum_GraphMap,Energies_RMS_ShowerHitNum_GraphMap,Energies_Mean_ShowerHitNum_Multi,Energies_RMS_ShowerHitNum_Multi,Energies_Mean_ShowerHitNum_canvasMulti,Energies_RMS_ShowerHitNum_canvasMulti,Energies_ShowerHitNum_canvasMap, ShowerHitNum_canvas,ShowerHitNum_2dHistMap,ShowerHitNum_2dCanvasMap,200,0,3000,100,0,fMaxSimEnergy,ShowerHitNum_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ShowerTotalEnergyDiff","(Reco Energy - True Energy)/True Energy",ShowerTotalEnergyDiff_HistMap,Energies_ShowerTotalEnergyDiff_HistMap,Energies_Mean_ShowerTotalEnergyDiff_GraphMap,Energies_RMS_ShowerTotalEnergyDiff_GraphMap,Energies_Mean_ShowerTotalEnergyDiff_Multi,Energies_RMS_ShowerTotalEnergyDiff_Multi,Energies_Mean_ShowerTotalEnergyDiff_canvasMulti,Energies_RMS_ShowerTotalEnergyDiff_canvasMulti,Energies_ShowerTotalEnergyDiff_canvasMap, ShowerTotalEnergyDiff_canvas,ShowerTotalEnergyDiff_2dHistMap,ShowerTotalEnergyDiff_2dCanvasMap,100,-1.2,0.2,100,0,fMaxSimEnergy,ShowerTotalEnergyDiff_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ShowerMag","|True Start - Reco Start| (cm)",ShowerMag_HistMap,Energies_ShowerMag_HistMap,Energies_Mean_ShowerMag_GraphMap,Energies_RMS_ShowerMag_GraphMap,Energies_Mean_ShowerMag_Multi,Energies_RMS_ShowerMag_Multi,Energies_Mean_ShowerMag_canvasMulti,Energies_RMS_ShowerMag_canvasMulti,Energies_ShowerMag_canvasMap, ShowerMag_canvas,ShowerMag_2dHistMap,ShowerMag_2dCanvasMap,100,0,50,100,0,fMaxSimEnergy,ShowerMag_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ShowerDirectionDiff","cos(theta)",ShowerDirectionDiff_HistMap,Energies_ShowerDirectionDiff_HistMap,Energies_Mean_ShowerDirectionDiff_GraphMap,Energies_RMS_ShowerDirectionDiff_GraphMap,Energies_Mean_ShowerDirectionDiff_Multi,Energies_RMS_ShowerDirectionDiff_Multi,Energies_Mean_ShowerDirectionDiff_canvasMulti,Energies_RMS_ShowerDirectionDiff_canvasMulti,Energies_ShowerDirectionDiff_canvasMap, ShowerDirectionDiff_canvas,ShowerDirectionDiff_2dHistMap,ShowerDirectionDiff_2dCanvasMap,100,-1,1,100,0,fMaxSimEnergy,ShowerDirectionDiff_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ShowerRecoEnergyVsTrueEnergyinRecoShower","TrueEnergyinRecoShower/ShowerRecoEnergy",ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap,Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap,Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap,Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap,Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_Multi,Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_Multi,Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMulti,Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMulti,Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMap, ShowerRecoEnergyVsTrueEnergyinRecoShower_canvas,ShowerRecoEnergyVsTrueEnergyinRecoShower_2dHistMap,ShowerRecoEnergyVsTrueEnergyinRecoShower_2dCanvasMap,100,0,3,100,0,fMaxSimEnergy,ShowerRecoEnergyVsTrueEnergyinRecoShower_TreeVal);
  
  ana::ShowerValidation::InitaliseGraphs("ShowerTrueEnergy","True Energy (MeV)",ShowerTrueEnergy_HistMap,Energies_ShowerTrueEnergy_HistMap,Energies_Mean_ShowerTrueEnergy_GraphMap,Energies_RMS_ShowerTrueEnergy_GraphMap,Energies_Mean_ShowerTrueEnergy_Multi,Energies_RMS_ShowerTrueEnergy_Multi,Energies_Mean_ShowerTrueEnergy_canvasMulti,Energies_RMS_ShowerTrueEnergy_canvasMulti,Energies_ShowerTrueEnergy_canvasMap, ShowerTrueEnergy_canvas,ShowerTrueEnergy_2dHistMap,ShowerTrueEnergy_2dCanvasMap,100,0,fMaxSimEnergy,100,0,fMaxSimEnergy,ShowerTrueEnergy_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("TrueEnergy","True Energy (MeV)",TrueEnergy_HistMap,Energies_TrueEnergy_HistMap,Energies_Mean_TrueEnergy_GraphMap,Energies_RMS_TrueEnergy_GraphMap,Energies_Mean_TrueEnergy_Multi,Energies_RMS_TrueEnergy_Multi,Energies_Mean_TrueEnergy_canvasMulti,Energies_RMS_TrueEnergy_canvasMulti,Energies_TrueEnergy_canvasMap, TrueEnergy_canvas,TrueEnergy_2dHistMap,TrueEnergy_2dCanvasMap,100,0,fMaxSimEnergy,100,0,fMaxSimEnergy,TrueEnergy_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("TrueHitNum","True number of hits in the true showers",TrueHitNum_HistMap,Energies_TrueHitNum_HistMap,Energies_Mean_TrueHitNum_GraphMap,Energies_RMS_TrueHitNum_GraphMap,Energies_Mean_TrueHitNum_Multi,Energies_RMS_TrueHitNum_Multi,Energies_Mean_TrueHitNum_canvasMulti,Energies_RMS_TrueHitNum_canvasMulti,Energies_TrueHitNum_canvasMap, TrueHitNum_canvas,TrueHitNum_2dHistMap,TrueHitNum_2dCanvasMap,200,0,5000,100,0,fMaxSimEnergy,TrueHitNum_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("GeoProjectionMatched","Projection matching score",GeoProjectionMatched_HistMap,Energies_GeoProjectionMatched_HistMap,Energies_Mean_GeoProjectionMatched_GraphMap,Energies_RMS_GeoProjectionMatched_GraphMap,Energies_Mean_GeoProjectionMatched_Multi,Energies_RMS_GeoProjectionMatched_Multi,Energies_Mean_GeoProjectionMatched_canvasMulti,Energies_RMS_GeoProjectionMatched_canvasMulti,Energies_GeoProjectionMatched_canvasMap, GeoProjectionMatched_canvas,GeoProjectionMatched_2dHistMap,GeoProjectionMatched_2dCanvasMap,200,0,1,100,0,fMaxSimEnergy,GeoProjectionMatched_TreeVal);


  ana::ShowerValidation::InitaliseGraphs("ClusterProjectionMatchedEnergy","Projection Matched",ClusterProjectionMatchedEnergy_HistMap,Energies_ClusterProjectionMatchedEnergy_HistMap,Energies_Mean_ClusterProjectionMatchedEnergy_GraphMap,Energies_RMS_ClusterProjectionMatchedEnergy_GraphMap,Energies_Mean_ClusterProjectionMatchedEnergy_Multi,Energies_RMS_ClusterProjectionMatchedEnergy_Multi,Energies_Mean_ClusterProjectionMatchedEnergy_canvasMulti,Energies_RMS_ClusterProjectionMatchedEnergy_canvasMulti,Energies_ClusterProjectionMatchedEnergy_canvasMap, ClusterProjectionMatchedEnergy_canvas,ClusterProjectionMatchedEnergy_2dHistMap,ClusterProjectionMatchedEnergy_2dCanvasMap,4,-1.5,2.5,100,0,fMaxSimEnergy,ClusterProjectionMatchedEnergy_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ClusterCompletenessEnergy","Energy Completness ",ClusterCompletenessEnergy_HistMap,Energies_ClusterCompletenessEnergy_HistMap,Energies_Mean_ClusterCompletenessEnergy_GraphMap,Energies_RMS_ClusterCompletenessEnergy_GraphMap,Energies_Mean_ClusterCompletenessEnergy_Multi,Energies_RMS_ClusterCompletenessEnergy_Multi,Energies_Mean_ClusterCompletenessEnergy_canvasMulti,Energies_RMS_ClusterCompletenessEnergy_canvasMulti,Energies_ClusterCompletenessEnergy_canvasMap, ClusterCompletenessEnergy_canvas,ClusterCompletenessEnergy_2dHistMap,ClusterCompletenessEnergy_2dCanvasMap,200,-0.2,1.2,200,0,fMaxSimEnergy,ClusterCompletenessEnergy_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ClusterCompletenessHits","Hit Completeness",ClusterCompletenessHits_HistMap,Energies_ClusterCompletenessHits_HistMap,Energies_Mean_ClusterCompletenessHits_GraphMap,Energies_RMS_ClusterCompletenessHits_GraphMap,Energies_Mean_ClusterCompletenessHits_Multi,Energies_RMS_ClusterCompletenessHits_Multi,Energies_Mean_ClusterCompletenessHits_canvasMulti,Energies_RMS_ClusterCompletenessHits_canvasMulti,Energies_ClusterCompletenessHits_canvasMap, ClusterCompletenessHits_canvas,ClusterCompletenessHits_2dHistMap,ClusterCompletenessHits_2dCanvasMap,200,-0.2,1.2,200,0,fMaxSimEnergy,ClusterCompletenessHits_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ClusterPurityEnergy","Energy Purity",ClusterPurityEnergy_HistMap,Energies_ClusterPurityEnergy_HistMap,Energies_Mean_ClusterPurityEnergy_GraphMap,Energies_RMS_ClusterPurityEnergy_GraphMap,Energies_Mean_ClusterPurityEnergy_Multi,Energies_RMS_ClusterPurityEnergy_Multi,Energies_Mean_ClusterPurityEnergy_canvasMulti,Energies_RMS_ClusterPurityEnergy_canvasMulti,Energies_ClusterPurityEnergy_canvasMap, ClusterPurityEnergy_canvas,ClusterPurityEnergy_2dHistMap,ClusterPurityEnergy_2dCanvasMap,200,0,1.2,200,0,fMaxSimEnergy,ClusterPurityEnergy_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ClusterPurityHits","Hit Purity",ClusterPurityHits_HistMap,Energies_ClusterPurityHits_HistMap,Energies_Mean_ClusterPurityHits_GraphMap,Energies_RMS_ClusterPurityHits_GraphMap,Energies_Mean_ClusterPurityHits_Multi,Energies_RMS_ClusterPurityHits_Multi,Energies_Mean_ClusterPurityHits_canvasMulti,Energies_RMS_ClusterPurityHits_canvasMulti,Energies_ClusterPurityHits_canvasMap, ClusterPurityHits_canvas,ClusterPurityHits_2dHistMap,ClusterPurityHits_2dCanvasMap,200,0,1.2,200,0,fMaxSimEnergy,ClusterPurityHits_TreeVal);
  
  ana::ShowerValidation::InitaliseGraphs("ClusterCompPurityEnergy","Energy Completeness*Purity",ClusterCompPurityEnergy_HistMap,Energies_ClusterCompPurityEnergy_HistMap,Energies_Mean_ClusterCompPurityEnergy_GraphMap,Energies_RMS_ClusterCompPurityEnergy_GraphMap,Energies_Mean_ClusterCompPurityEnergy_Multi,Energies_RMS_ClusterCompPurityEnergy_Multi,Energies_Mean_ClusterCompPurityEnergy_canvasMulti,Energies_RMS_ClusterCompPurityEnergy_canvasMulti,Energies_ClusterCompPurityEnergy_canvasMap, ClusterCompPurityEnergy_canvas,ClusterCompPurityEnergy_2dHistMap,ClusterCompPurityEnergy_2dCanvasMap,100,-0.2,1.2,200,0,fMaxSimEnergy,ClusterCompPurityEnergy_TreeVal);

  ana::ShowerValidation::InitaliseGraphs("ClusterCompPurityHits","Hit Completeness*Purity",ClusterCompPurityHits_HistMap,Energies_ClusterCompPurityHits_HistMap,Energies_Mean_ClusterCompPurityHits_GraphMap,Energies_RMS_ClusterCompPurityHits_GraphMap,Energies_Mean_ClusterCompPurityHits_Multi,Energies_RMS_ClusterCompPurityHits_Multi,Energies_Mean_ClusterCompPurityHits_canvasMulti,Energies_RMS_ClusterCompPurityHits_canvasMulti,Energies_ClusterCompPurityHits_canvasMap, ClusterCompPurityHits_canvas,ClusterCompPurityHits_2dHistMap,ClusterCompPurityHits_2dCanvasMap,200,-0.2,1.2,200,0,fMaxSimEnergy,ClusterCompPurityHits_TreeVal);

  ana::ShowerValidation::InitaliseHitGraphs("HitCompletenessEnergy","Completeness",HitCompletenessEnergy_HistMap,Energies_HitCompletenessEnergy_HistMap,Energies_Mean_HitCompletenessEnergy_GraphMap,Energies_RMS_HitCompletenessEnergy_GraphMap,Energies_Mean_HitCompletenessEnergy_Multi,Energies_RMS_HitCompletenessEnergy_Multi,Energies_Mean_HitCompletenessEnergy_canvasMulti,Energies_RMS_HitCompletenessEnergy_canvasMulti,Energies_HitCompletenessEnergy_canvasMap, HitCompletenessEnergy_canvas,HitCompletenessEnergy_2dHistMap,HitCompletenessEnergy_2dCanvasMap,200,-0.2,1.2,200,0,fMaxSimEnergy,HitCompletenessEnergy_TreeVal);

}


void ana::ShowerValidation::analyze(const art::Event& evt) {

  //Fill the Event Number info 
  EventNum_name_TreeVal.push_back(evt.run());
  EventNum_name_TreeVal.push_back(evt.subRun());
  EventNum_name_TreeVal.push_back(evt.event());

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

    if(fVerbose > 1){std::cout << "True Particle with track ID: " << particle->TrackId() << " Has code of: " << particle->PdgCode() << " and Energy of: " << particle->E() << " With Mother: " << particle->Mother() << " Proccess: " << particle->Process() << " End Process: "  << particle->EndProcess() << std::endl;}


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
    if(TMath::Abs(pdgcode) == 11 || pdgcode == 22){
  
      if(motherparticle->E() > fSimEnergyCut){
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
	  if(Hit_Density > fDensityCut){
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

  //If there are no true showers we can't validate 
  if(ShowersMothers.size() == 0){return;}

  //If we are looking at a pion then fill the smallest and largest showers
  for (sim::ParticleList::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt){
    
    const simb::MCParticle *particle = particleIt->second;

    float SmallestShowerE_Before = 999;
    float BiggestShowerE_Before  = -999;
    float SmallestShowerE_After = 999;
    float BiggestShowerE_After  = -999;

    if(particle->PdgCode() == 111){      
      
      int NumDaughters = particle->NumberDaughters();
      for(int daughter=0; daughter<NumDaughters; ++daughter){
	
	//Get the daughter track ID
	int DaughterID = particle->Daughter(daughter);

	//Get the largest and smallest energies before the cuts
	if(trueParticles[DaughterID]->E() < SmallestShowerE_Before){SmallestShowerE_Before= trueParticles[DaughterID]->E();}
	if(trueParticles[DaughterID]->E() > BiggestShowerE_Before) {BiggestShowerE_Before = trueParticles[DaughterID]->E();}

	//Get the largest and smallest energies after the cuts 
	if(ShowersMothers.find(DaughterID) != ShowersMothers.end()){
	  if(trueParticles[DaughterID]->E() < SmallestShowerE_After){SmallestShowerE_After = trueParticles[DaughterID]->E();}
	  if(trueParticles[DaughterID]->E() > BiggestShowerE_After) {BiggestShowerE_After  = trueParticles[DaughterID]->E();}
	}
      }
    }

    if(particle->PdgCode() == 11 && particle->TrackId() == 0){
      SmallestShowerE_Before  = particle->E();
      BiggestShowerE_Before = particle->E();
      SmallestShowerE_After = particle->E();
      BiggestShowerE_After  = particle->E();
    }

      //Fill the plots
      BiggestShowerMotherE_BeforeCutHist ->Fill(BiggestShowerE_Before*1000);
      SmallestShowerMotherE_BeforeCutHist->Fill(SmallestShowerE_Before*1000);
      BiggestShowerMotherE_BeforeCut_TreeVal.push_back(BiggestShowerE_Before*1000);
      SmallestShowerMotherE_BeforeCut_TreeVal.push_back(SmallestShowerE_Before*1000);

      BiggestShowerMotherE_AfterCutHist ->Fill(BiggestShowerE_After*1000);
      SmallestShowerMotherE_AfterCutHist->Fill(SmallestShowerE_After*1000);
      BiggestShowerMotherE_AfterCut_TreeVal.push_back(BiggestShowerE_After*1000);
      SmallestShowerMotherE_AfterCut_TreeVal.push_back(SmallestShowerE_After*1000);
      //}  
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

  //######################
  //### Hit Validation ###
  //######################
 
  //Loop over the hit handles 
  for(unsigned int hitlab_it=0; hitlab_it<fHitModuleLabels.size(); ++hitlab_it){
    
    //Getting the Hit Information                                                       
    art::Handle<std::vector<recob::Shower> > hitValHandle;
    std::vector<art::Ptr<recob::Shower> > hits_fromValhandle;
    std::string fHitModuleLabel = fHitModuleLabels[hitlab_it];
    if(evt.getByLabel(fHitModuleLabel,hitValHandle))
      {art::fill_ptr_vector(hits_fromValhandle,hitValHandle);}
    
    //Calculate the total energy deposited
    float TotalEnergyDeposited = 0;
    for(auto track_iter: MCTrack_Energy_map){
      TotalEnergyDeposited += track_iter.second;
    }
    
    //Calculate how much energy was deposited in the hits 
    for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){
      
      float TotalEnergyDepinHits = RecoUtils::TotalEnergyDepinHits(hits,plane_id.Plane);
      if(TotalEnergyDepinHits != 0){ 
	HitCompletenessEnergy_HistMap[fHitModuleLabel][plane_id]->Fill(( TotalEnergyDepinHits)/TotalEnergyDeposited);
	HitCompletenessEnergy_2dHistMap[fHitModuleLabel][plane_id]->Fill((TotalEnergyDepinHits)/TotalEnergyDeposited,simenergy*1000);
	HitCompletenessEnergy_TreeVal[fHitModuleLabel][plane_id.Plane].push_back((TotalEnergyDepinHits)/TotalEnergyDeposited);
	if(fVerbose > 1){std::cout << "Hit Completeness:"  << (TotalEnergyDepinHits)/TotalEnergyDeposited << std::endl;}
	
	if(fEnergies.size() != 0){
	  for(unsigned int i=0; i<fEnergies.size(); ++i){
	    if(TMath::Abs(simenergy*1000 - fEnergies[i])< fEnergyWidth){
	      Energies_HitCompletenessEnergy_HistMap[fHitModuleLabel][fEnergies[i]][plane_id]->Fill((TotalEnergyDepinHits)/TotalEnergyDeposited);
	    }
	  }
	}
      }
      else{
	HitCompletenessEnergy_TreeVal[fHitModuleLabel][plane_id.Plane].push_back(-99999);
      }
    }
  }


  if(fVerbose > 1){std::cout << "Hit Validation Complete" << std::endl;} 
		
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

    std::vector< art::Ptr<recob::Hit> > showerhits; //hits in the shower 

    //Get the ID of the shower hit module
    art::ProductID showerhit_productid = fmh.at(0).front().id();
   
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

    std::map<int,float> MinStartDiff;
    std::map<int,float> MinStartX;
    std::map<int,float> MinStartY;
    std::map<int,float> MinStartZ;
    std::map<int,float> MinShowerDirection_Xdiff;
    std::map<int,float> MinShowerDirection_Ydiff;
    std::map<int,float> MinShowerDirection_Zdiff;
    std::map<int,float> MinShowerDirection_diff;
    std::map<int,int>   MinEvaluateShowerDirection;
    std::map<int,int>   MinEvaluateShowerStart;

    //Initialise the minimum values for each shower mother
    for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end();++showermother){
         MinStartDiff[showermother->first] = 99999;
	 MinStartX[showermother->first] = -99999;
	 MinStartY[showermother->first] = -99999;
	 MinStartZ[showermother->first] = -99999;
	 MinShowerDirection_Xdiff[showermother->first] = -99999;
	 MinShowerDirection_Ydiff[showermother->first] = -99999;
	 MinShowerDirection_Zdiff[showermother->first] = -99999;
	 MinShowerDirection_diff[showermother->first] = -99999;
	 MinEvaluateShowerDirection[showermother->first] = 0;
	 MinEvaluateShowerStart[showermother->first] = 0;
    }

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


      if((int) showerhits.size() < fMinHitSize) {continue;}

      int ShowerBest_Plane = shower->best_plane();
      if(std::isnan(ShowerBest_Plane)){ShowerBest_Plane = 0;}

      //Function from RecoUtils, finds the most probable track ID associated with the set of hits from there true energy depositons. The pair returns the energy as well. 
      std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain(ShowersMothers,showerhits,ShowerBest_Plane); 

      int ShowerTrackID = ShowerTrackInfo.first;
      double TrueEnergyDepWithinShower_FromTrueShower = ShowerTrackInfo.second;

      //Check to see if the shower was correctly matched.
      if(TrueEnergyDepWithinShower_FromTrueShower == -99999){
	continue;
	std::cout << "Reco Shower not matched to true shower. Think of reducing the energy threshold" << std::endl; 
      }

      //Get the number of hits associated to the true particle. 
      int TrueHitDep_FromTrueShower = 0;
      for(std::vector<int>::iterator track_id=ShowersMothers[ShowerTrackID].begin(); track_id!=ShowersMothers[ShowerTrackID].end(); ++track_id){
	for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){
	  TrueHitDep_FromTrueShower +=  MCTrack_hit_map[showerhit_productid][*track_id][plane_id];
	}
      }

      int NumberofHitsinRecoShower = showerhits.size();

      //Get the number of hits in the reco shower from the true shower.
      int TrueHitsDep_WithinRecoShower = 0;
      std::map<int,std::map<geo::PlaneID,int> >  MCTrack_showerhit_map = RecoUtils::NumberofPlaneHitsPerTrack(showerhits);
      for(auto const& planehit_map : MCTrack_showerhit_map){
	if(std::find(ShowersMothers[ShowerTrackID].begin(), ShowersMothers[ShowerTrackID].end(),planehit_map.first) == ShowersMothers[ShowerTrackID].end()){continue;}
	for(auto const& hit_num : planehit_map.second){
	  TrueHitsDep_WithinRecoShower += hit_num.second;
	}
      }

      
      double TrueEnergyDep_FromShower = 0;
      //Calculate the true Energy deposited By Shower 
      for(std::vector<int>::iterator daughterID=ShowersMothers[ShowerTrackID].begin(); daughterID!=ShowersMothers[ShowerTrackID].end(); ++daughterID){
	TrueEnergyDep_FromShower += MCTrack_Energy_map[*daughterID];
      }

      double hitcompleteness = 0;
      if(TrueEnergyDep_FromShower != 0){
	hitcompleteness = ((double) TrueHitsDep_WithinRecoShower)/(double) TrueHitDep_FromTrueShower;
	ShowerHitsCompleteness_HistMap[fShowerModuleLabel]->Fill(hitcompleteness);
	ShowerHitsCompleteness_2dHistMap[fShowerModuleLabel]->Fill(hitcompleteness,simenergy*1000);
	ShowerHitsCompleteness_TreeVal[fShowerModuleLabel].push_back(hitcompleteness);
      }
      
      double hitpurity = 0;
      if(TrueHitsDep_WithinRecoShower != 0){
	hitpurity   =   (double) TrueHitsDep_WithinRecoShower/(double) NumberofHitsinRecoShower;
	ShowerHitsPurity_HistMap[fShowerModuleLabel]->Fill(hitpurity);
	ShowerHitsPurity_2dHistMap[fShowerModuleLabel]->Fill(hitpurity,simenergy*1000);
	ShowerHitsPurity_TreeVal[fShowerModuleLabel].push_back(hitpurity);
      }
      
      //Energy deposited within the set of Hits associated to the shower.
      double TrueEnergyDep_WithinRecoShower = 0; 

      //Loop over the hits and find the IDEs
      if(fVerbose > 1){std::cout << "shower hits size: " << showerhits.size() << std::endl;}
      for(std::vector< art::Ptr<recob::Hit> >::iterator hitIt=showerhits.begin(); hitIt!=showerhits.end(); ++hitIt){
	
	//Get the plane ID 
	geo::WireID wireid = (*hitIt)->WireID();
	int PlaneID = wireid.Plane;

	if(PlaneID != ShowerBest_Plane){continue;}

	//Split the Hit into its IDE for each track it associates with.                               
	std::vector<sim::TrackIDE> trackIDEs = backtracker->HitToTrackIDEs((*hitIt));
	for (unsigned int idIt = 0; idIt < trackIDEs.size(); ++idIt){ 
	  
	  //Find the true total energy deposited in a set of hits. 
	  TrueEnergyDep_WithinRecoShower += trackIDEs.at(idIt).energy;
	}
      }//Hit Loop
      
      double energycompleteness = 0;
      if(TrueEnergyDep_FromShower != 0){
	energycompleteness =  (TrueEnergyDepWithinShower_FromTrueShower)/TrueEnergyDep_FromShower;
	ShowerEnergyCompleteness_HistMap[fShowerModuleLabel]->Fill(energycompleteness);
	ShowerEnergyCompleteness_2dHistMap[fShowerModuleLabel]->Fill(energycompleteness,simenergy*1000);
	ShowerEnergyCompleteness_TreeVal[fShowerModuleLabel].push_back(energycompleteness);
      }
      
      double energypurity = 0;
      if(TrueEnergyDep_WithinRecoShower != 0){
	energypurity       =  TrueEnergyDepWithinShower_FromTrueShower/TrueEnergyDep_WithinRecoShower;
	ShowerEnergyPurity_HistMap[fShowerModuleLabel]->Fill(energypurity);
	ShowerEnergyPurity_2dHistMap[fShowerModuleLabel]->Fill(energypurity,simenergy*1000);
	ShowerEnergyPurity_TreeVal[fShowerModuleLabel].push_back(energypurity);
      }
      
      //Find the MCParticle this shower associates to
      const simb::MCParticle* MCShowerParticle = trueParticles.at(ShowerTrackID);
      
      //Find the Energy of the particle: 
      float Energy = MCShowerParticle->E();
      ShowerTrueEnergy_HistMap[fShowerModuleLabel]->Fill(Energy*1000);
      ShowerTrueEnergy_TreeVal[fShowerModuleLabel].push_back(Energy*1000);
      ShowerTrueEnergy_2dHistMap[fShowerModuleLabel]->Fill(Energy,simenergy*1000);

      ShowerStartEndProcess_TreeVal[fShowerModuleLabel].push_back(MCShowerParticle->EndProcess());

      if(std::find(SmallestShowerMotherE_AfterCut_TreeVal.begin(),SmallestShowerMotherE_AfterCut_TreeVal.end(),Energy*1000) != SmallestShowerMotherE_AfterCut_TreeVal.end()){AsscoiatedSmallestShowerMotherEHist->Fill(Energy*1000);}
    else if(std::find(BiggestShowerMotherE_AfterCut_TreeVal.begin(),BiggestShowerMotherE_AfterCut_TreeVal.end(),Energy*1000) != BiggestShowerMotherE_AfterCut_TreeVal.end()){AsscoiatedBiggestShowerMotherEHist->Fill(Energy*1000);}
      

      //Get the number of Traj points to loop over                                                    
      unsigned int TrajPoints = MCShowerParticle->NumberTrajectoryPoints();
      
      //Find the start and end points of the initial particle. 
      TLorentzVector PositionTrajStart =  MCShowerParticle->Position(0);
      TLorentzVector PositionTrajEnd   =  MCShowerParticle->Position(TrajPoints-1);
      
      //The Start of position for the electron shower is PositionTrajStart but the start position for photon showers is at the end of the shower (the start of the e+- track) 
      if(MCShowerParticle->PdgCode() == 22){
	PositionTrajStart = PositionTrajEnd;
      }
      
      //The three vector for track length is the shower direction 
      //TVector3  TrueShowerDirection = (PositionTrajEnd - PositionTrajStart).Vect();
      TVector3  TrueShowerDirection(MCShowerParticle->Px(), MCShowerParticle->Py(),MCShowerParticle->Pz());

      //Initial track lentgh of the shower.
      double TrueTrackLength = TrueShowerDirection.Mag();
      
      //Get the information for the shower  
      TVector3  ShowerDirection                  = shower->Direction();
      TVector3  ShowerStart                      = shower->ShowerStart();//cm
      double ShowerTrackLength                   = shower->Length();//cm        
      std::vector< double >  ShowerEnergyPlanes  = shower->Energy();//MeV
      std::vector< double >  ShowerdEdX_vec      = shower->dEdx();//MeV/cm  

      //Bools to fill metric histrograms wheen needed.
      bool EvaluateShowerDirection       = false;
      bool EvaluateShowerStart           = false;
      bool EvaluateShowerLength          = false;
      bool EvaluateShowerEnergy          = false;
      bool EvaluateShowerdEdx            = false;
      bool EvalulateGeoProjectionMatched = false;

      //Evaulate 3D Shower Reconstruction Dependent Metrics
      if(!std::isnan(ShowerDirection.X()) || ShowerDirection.Mag() == 0) {EvaluateShowerDirection = true; ++MinEvaluateShowerDirection[ShowerTrackID];}
      if(!std::isnan(ShowerStart.X()))                                   {EvaluateShowerStart     = true; ++MinEvaluateShowerStart[ShowerTrackID];}
      if(!std::isnan(ShowerTrackLength) || shower->has_length())         {EvaluateShowerLength    = true;}
      if(ShowerEnergyPlanes.size() != 0){
	if(!std::isnan(ShowerEnergyPlanes.at(0))){
	  EvaluateShowerEnergy = true;
	}
      }
      if(ShowerdEdX_vec.size() != 0){
	if(!std::isnan(ShowerdEdX_vec.at(0))){
	  EvaluateShowerdEdx = true;
	}
      }
   
      //Get the angles between the direction
      float ShowerDirection_Xdiff = -99999;
      float ShowerDirection_Ydiff = -99999;
      float ShowerDirection_Zdiff = -99999;
      float ShowerDirection_diff  = -99999;

      std::vector<double> ShowerEnergyPlanes_remove(3);
      
	ShowerEnergyPlanes_remove[0] = ShowerEnergyPlanes[0];
      	ShowerEnergyPlanes_remove[1] = ShowerEnergyPlanes[1];
      	ShowerEnergyPlanes_remove[2] = (ShowerEnergyPlanes[2] - 0.00155171)*0.00155171/4.39964 + 4.39964;
      
      if(TMath::Sqrt((TrueShowerDirection.Y()*TrueShowerDirection.Y() + TrueShowerDirection.Z()*TrueShowerDirection.Z()))*TMath::Sqrt((ShowerDirection.Y()*ShowerDirection.Y() + ShowerDirection.Z()*ShowerDirection.Z())) !=0){
	ShowerDirection_Xdiff = (TrueShowerDirection.Y()*ShowerDirection.Y() + TrueShowerDirection.Z()*ShowerDirection.Z())/(TMath::Sqrt((TrueShowerDirection.Y()*TrueShowerDirection.Y() + TrueShowerDirection.Z()*TrueShowerDirection.Z()))*TMath::Sqrt((ShowerDirection.Y()*ShowerDirection.Y() + ShowerDirection.Z()*ShowerDirection.Z())));
      }

      if(TMath::Sqrt((TrueShowerDirection.X()*TrueShowerDirection.X() + TrueShowerDirection.Z()*TrueShowerDirection.Z()))*TMath::Sqrt((ShowerDirection.X()*ShowerDirection.X() + ShowerDirection.Z()*ShowerDirection.Z())) !=0){
	ShowerDirection_Ydiff = (TrueShowerDirection.X()*ShowerDirection.X() + TrueShowerDirection.Z()*ShowerDirection.Z())/(TMath::Sqrt((TrueShowerDirection.X()*TrueShowerDirection.X() + TrueShowerDirection.Z()*TrueShowerDirection.Z()))*TMath::Sqrt((ShowerDirection.X()*ShowerDirection.X() + ShowerDirection.Z()*ShowerDirection.Z())));
      }
      
      if(TMath::Sqrt((TrueShowerDirection.Y()*TrueShowerDirection.Y() + TrueShowerDirection.X()*TrueShowerDirection.X()))*TMath::Sqrt((ShowerDirection.Y()*ShowerDirection.Y() + ShowerDirection.X()*ShowerDirection.X())) !=0){
	ShowerDirection_Zdiff = (TrueShowerDirection.Y()*ShowerDirection.Y() + TrueShowerDirection.X()*ShowerDirection.X())/(TMath::Sqrt((TrueShowerDirection.Y()*TrueShowerDirection.Y() + TrueShowerDirection.X()*TrueShowerDirection.X()))*TMath::Sqrt((ShowerDirection.Y()*ShowerDirection.Y() + ShowerDirection.X()*ShowerDirection.X())));
      }
	
      if(TrueShowerDirection.Mag() != 0 || ShowerDirection.Mag() !=0){ 
	ShowerDirection_diff  = TrueShowerDirection.Dot(ShowerDirection)/(TrueShowerDirection.Mag()*ShowerDirection.Mag());
      }

      //Get the Error in the position. intialised as 0,0,0 this is a problem here. 
      double Start_diff = TMath::Sqrt(TMath::Power(PositionTrajStart.X()-ShowerStart.X(),2) + TMath::Power(PositionTrajStart.Y()-ShowerStart.Y(),2) + TMath::Power(PositionTrajStart.Z()-ShowerStart.Z(),2));
   
      //Fill the histograms.
      if(fFillOnlyClosestShower){
	if(Start_diff < MinStartDiff[ShowerTrackID]){
	  if(EvaluateShowerStart){
	    MinStartDiff[ShowerTrackID] = Start_diff;
	    MinStartX[ShowerTrackID]    = TMath::Abs(PositionTrajStart.X()-ShowerStart.X());
	    MinStartY[ShowerTrackID]    = TMath::Abs(PositionTrajStart.Y()-ShowerStart.Y());
	    MinStartZ[ShowerTrackID]    = TMath::Abs(PositionTrajStart.Z()-ShowerStart.Z());
	  }
	    if(EvaluateShowerDirection){
	    MinShowerDirection_Xdiff[ShowerTrackID] = ShowerDirection_Xdiff;
	    MinShowerDirection_Ydiff[ShowerTrackID] = ShowerDirection_Ydiff;
	    MinShowerDirection_Zdiff[ShowerTrackID] = ShowerDirection_Zdiff;
	    MinShowerDirection_diff[ShowerTrackID]  = ShowerDirection_diff;
	  }
	}
      }
      else{
	if(EvaluateShowerDirection){
	  ShowerDirection_X_HistMap[fShowerModuleLabel]->Fill(ShowerDirection_Xdiff);
	  ShowerDirection_Y_HistMap[fShowerModuleLabel]->Fill(ShowerDirection_Ydiff);
	  ShowerDirection_Z_HistMap[fShowerModuleLabel]->Fill(ShowerDirection_Zdiff);
	  ShowerDirectionDiff_HistMap[fShowerModuleLabel]->Fill(ShowerDirection_diff);
	  ShowerDirection_X_2dHistMap[fShowerModuleLabel]->Fill(ShowerDirection_Xdiff,simenergy*1000);
	  ShowerDirection_Y_2dHistMap[fShowerModuleLabel]->Fill(ShowerDirection_Ydiff,simenergy*1000);
	  ShowerDirection_Z_2dHistMap[fShowerModuleLabel]->Fill(ShowerDirection_Zdiff,simenergy*1000);
	  ShowerDirectionDiff_2dHistMap[fShowerModuleLabel]->Fill(ShowerDirection_diff,simenergy*1000);
	}
	
	if(EvaluateShowerStart){
	  ShowerStart_X_HistMap[fShowerModuleLabel]->Fill(TMath::Abs(PositionTrajStart.X()-ShowerStart.X()));
	  ShowerStart_Y_HistMap[fShowerModuleLabel]->Fill(TMath::Abs(PositionTrajStart.Y()-ShowerStart.Y()));
	  ShowerStart_Z_HistMap[fShowerModuleLabel]->Fill(TMath::Abs(PositionTrajStart.Z()-ShowerStart.Z()));
	  ShowerMag_HistMap[fShowerModuleLabel]->Fill(Start_diff);
	  ShowerStart_X_2dHistMap[fShowerModuleLabel]->Fill(TMath::Abs(PositionTrajStart.X()-ShowerStart.X()),simenergy*1000);
	  ShowerStart_Y_2dHistMap[fShowerModuleLabel]->Fill(TMath::Abs(PositionTrajStart.Y()-ShowerStart.Y()),simenergy*1000);
	  ShowerStart_Z_2dHistMap[fShowerModuleLabel]->Fill(TMath::Abs(PositionTrajStart.Z()-ShowerStart.Z()),simenergy*1000);
	  ShowerMag_2dHistMap[fShowerModuleLabel]->Fill(Start_diff,simenergy*1000);
  	}
      }      
    
      if(EvaluateShowerDirection){
	ShowerDirection_X_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_Xdiff);
	ShowerDirection_Y_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_Ydiff);
	ShowerDirection_Z_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_Zdiff);
	ShowerDirectionDiff_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_diff);
      }
      else{
	ShowerDirection_X_TreeVal[fShowerModuleLabel].push_back(-99999);
	ShowerDirection_Y_TreeVal[fShowerModuleLabel].push_back(-99999);
	ShowerDirection_Z_TreeVal[fShowerModuleLabel].push_back(-99999);
	ShowerDirectionDiff_TreeVal[fShowerModuleLabel].push_back(-99999);
      }
      
      if(EvaluateShowerStart){
	ShowerStart_X_TreeVal[fShowerModuleLabel].push_back(TMath::Abs(PositionTrajStart.X()-ShowerStart.X()));
	ShowerStart_Y_TreeVal[fShowerModuleLabel].push_back(TMath::Abs(PositionTrajStart.Y()-ShowerStart.Y()));
	ShowerStart_Z_TreeVal[fShowerModuleLabel].push_back(TMath::Abs(PositionTrajStart.Z()-ShowerStart.Z()));
	ShowerMag_TreeVal[fShowerModuleLabel].push_back(Start_diff);
      }
      else{
	ShowerStart_X_TreeVal[fShowerModuleLabel].push_back(-99999);
	ShowerStart_Y_TreeVal[fShowerModuleLabel].push_back(-99999);
	ShowerStart_Z_TreeVal[fShowerModuleLabel].push_back(-99999);
	ShowerMag_TreeVal[fShowerModuleLabel].push_back(-99999);
      }
   
      if(EvaluateShowerStart && EvaluateShowerDirection){
	PosDir_2dHistMap[fShowerModuleLabel]->Fill(Start_diff,TrueShowerDirection.Dot(ShowerDirection)/(TrueShowerDirection.Mag()*ShowerDirection.Mag()));
      }

      if(EvaluateShowerLength){
	ShowerLength_HistMap[fShowerModuleLabel]->Fill(TMath::Abs(TrueTrackLength-ShowerTrackLength));
	ShowerLength_TreeVal[fShowerModuleLabel].push_back(TMath::Abs(TrueTrackLength-ShowerTrackLength));
      }
      else{      
	ShowerLength_TreeVal[fShowerModuleLabel].push_back(-99999);
      }

      if(TrueEnergyDep_FromShower != 0 && EvaluateShowerEnergy){
	ShowerEnergyDiff_HistMap[fShowerModuleLabel]->Fill(ShowerEnergyPlanes_remove[ShowerBest_Plane]/TrueEnergyDep_FromShower);
 	ShowerTotalEnergyDiff_HistMap[fShowerModuleLabel]->Fill((ShowerEnergyPlanes_remove[ShowerBest_Plane] - TrueEnergyDep_FromShower)/TrueEnergyDep_FromShower);
 	ShowerEnergyDiff_TreeVal[fShowerModuleLabel].push_back(ShowerEnergyPlanes_remove[ShowerBest_Plane]/TrueEnergyDep_FromShower);
 	ShowerTotalEnergyDiff_TreeVal[fShowerModuleLabel].push_back((ShowerEnergyPlanes_remove[ShowerBest_Plane] - TrueEnergyDep_FromShower)/TrueEnergyDep_FromShower);
      }
      else{
	ShowerEnergyDiff_TreeVal[fShowerModuleLabel].push_back(-99999);
 	ShowerTotalEnergyDiff_TreeVal[fShowerModuleLabel].push_back(-99999);
      }
      
      if(TrueEnergyDepWithinShower_FromTrueShower != 0 && EvaluateShowerEnergy){
 	ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap[fShowerModuleLabel]->Fill(ShowerEnergyPlanes_remove[ShowerBest_Plane]/TrueEnergyDepWithinShower_FromTrueShower);
	ShowerRecoEnergyVsTrueEnergyinRecoShower_TreeVal[fShowerModuleLabel].push_back(ShowerEnergyPlanes_remove[ShowerBest_Plane]/TrueEnergyDepWithinShower_FromTrueShower);
      }
      else{
	ShowerRecoEnergyVsTrueEnergyinRecoShower_TreeVal[fShowerModuleLabel].push_back(-99999);
      }

      if(EvaluateShowerEnergy){
	ShowerEnergy_HistMap[fShowerModuleLabel]->Fill(ShowerEnergyPlanes_remove[ShowerBest_Plane]);
	ShowerEnergy_TreeVal[fShowerModuleLabel].push_back(ShowerEnergyPlanes_remove[ShowerBest_Plane]);
      }
      else{
	ShowerEnergy_TreeVal[fShowerModuleLabel].push_back(-99999);
      }

      ShowerHitNum_HistMap[fShowerModuleLabel]->Fill(showerhits.size());
      ShowerHitNum_TreeVal[fShowerModuleLabel].push_back(showerhits.size());

      if(EvaluateShowerdEdx){
 	ShowerdEdx_HistMap[fShowerModuleLabel]->Fill((ShowerdEdX_vec[ShowerBest_Plane]));
	ShowerdEdx_TreeVal[fShowerModuleLabel].push_back((ShowerdEdX_vec[ShowerBest_Plane]));
      }
      else{
	ShowerdEdx_TreeVal[fShowerModuleLabel].push_back(-99999);
      }
    
      //Fill the 2D histograms
      if(EvaluateShowerLength){
	ShowerLength_2dHistMap[fShowerModuleLabel]->Fill(TMath::Abs(TrueTrackLength-ShowerTrackLength),simenergy*1000);
      }

      if(TrueEnergyDep_FromShower != 0 && EvaluateShowerEnergy){
	ShowerEnergyDiff_2dHistMap[fShowerModuleLabel]->Fill(ShowerEnergyPlanes_remove[ShowerBest_Plane]/TrueEnergyDep_FromShower,simenergy*1000);
	ShowerTotalEnergyDiff_2dHistMap[fShowerModuleLabel]->Fill((ShowerEnergyPlanes_remove[ShowerBest_Plane] - TrueEnergyDep_FromShower)/TrueEnergyDep_FromShower,simenergy*1000);
      }

      if(TrueEnergyDepWithinShower_FromTrueShower != 0 && EvaluateShowerEnergy){
	ShowerRecoEnergyVsTrueEnergyinRecoShower_2dHistMap[fShowerModuleLabel]->Fill(ShowerEnergyPlanes_remove[ShowerBest_Plane]/TrueEnergyDepWithinShower_FromTrueShower,simenergy*1000);
      }

      if(EvaluateShowerEnergy){
	ShowerEnergy_2dHistMap[fShowerModuleLabel]->Fill(ShowerEnergyPlanes_remove[ShowerBest_Plane],simenergy*1000);
      }
     
      ShowerHitNum_2dHistMap[fShowerModuleLabel]->Fill(showerhits.size(),simenergy*1000);
      
      if(EvaluateShowerdEdx){
	ShowerdEdx_2dHistMap[fShowerModuleLabel]->Fill((ShowerdEdX_vec[ShowerBest_Plane]),simenergy*1000);
      }
    
      //Fill the Energy dependent Histograms. 
      if(fEnergies.size() != 0){
	for(unsigned int i=0; i<fEnergies.size(); ++i){
	  if(TMath::Abs(simenergy*1000 - fEnergies[i])< fEnergyWidth){
	    
	    if(!fFillOnlyClosestShower){
	      
	      if(EvaluateShowerDirection){
		Energies_ShowerDirection_X_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(ShowerDirection_Xdiff);
		Energies_ShowerDirection_Y_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(ShowerDirection_Ydiff);
		Energies_ShowerDirection_Z_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(ShowerDirection_Zdiff);
		Energies_ShowerDirectionDiff_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(ShowerDirection_diff);
	      }
	      if(EvaluateShowerStart){
		Energies_ShowerStart_X_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(TMath::Abs(PositionTrajStart.X()-ShowerStart.X()));
		Energies_ShowerStart_Y_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(TMath::Abs(PositionTrajStart.Y()-ShowerStart.Y()));
		Energies_ShowerStart_Z_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(TMath::Abs(PositionTrajStart.Z()-ShowerStart.Z()));
		Energies_ShowerMag_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(Start_diff );
		Energies_ShowerDirectionDiff_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(ShowerDirection_diff);
	      }
	    }
	    
	    if(EvaluateShowerEnergy){
	      Energies_ShowerEnergy_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(ShowerEnergyPlanes_remove[ShowerBest_Plane]);
	    }

	    if(EvaluateShowerLength){
	      Energies_ShowerLength_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(TMath::Abs(TrueTrackLength-ShowerTrackLength));
	    }

	    if(TrueEnergyDep_FromShower != 0 && EvaluateShowerEnergy){
	      Energies_ShowerEnergyDiff_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(ShowerEnergyPlanes_remove[ShowerBest_Plane]/TrueEnergyDep_FromShower);
	      Energies_ShowerTotalEnergyDiff_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill((ShowerEnergyPlanes_remove[ShowerBest_Plane] - TrueEnergyDep_FromShower)/TrueEnergyDep_FromShower,simenergy*1000);
	    }
	    
	    if(TrueEnergyDepWithinShower_FromTrueShower != 0 && EvaluateShowerEnergy){
	      Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(ShowerEnergyPlanes_remove[ShowerBest_Plane]/3/TrueEnergyDepWithinShower_FromTrueShower);
	    }
	 
	    Energies_ShowerHitNum_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(showerhits.size());
	    Energies_ShowerTrueEnergy_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(Energy);

	    if(TrueEnergyDep_FromShower != 0)      {Energies_ShowerEnergyCompleteness_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(energycompleteness);}
	    if(TrueEnergyDep_WithinRecoShower != 0){Energies_ShowerEnergyPurity_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(energypurity);}
	    if(TrueEnergyDep_FromShower != 0)      {Energies_ShowerHitsCompleteness_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(hitcompleteness);}
	    if(TrueHitsDep_WithinRecoShower != 0)  {Energies_ShowerHitsPurity_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(hitpurity);}
	  
	    if(EvaluateShowerdEdx){
	      Energies_ShowerdEdx_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill((ShowerdEdX_vec[ShowerBest_Plane]));
	    }
	  }
	}
      }
        
      if(fVerbose > 0){
	std::cout << "#################################################" << std::endl;
	std::cout << "Global Event Information" << std::endl;
	std::cout << "#################################################" << std::endl;
	std::cout << "Shower Label: " <<  fShowerModuleLabel << std::endl;
	std::cout << "Number of Reco showers: " << showers.size() << std::endl;
	std::cout << "Energy Simulated: " << simenergy << std::endl;
	std::cout << "Number of True Showers That pass the E cut: " <<  num_of_showers_viaEcut << std::endl;
	std::cout << "Number of True Showers That pass the Density cut: " << num_of_showers_viaDensitycut << std::endl;
	std::cout << "Number of True Showers in the Event: " << ShowersMothers.size() << std::endl;
	std::cout << "#################################################" << std::endl;
	std::cout << "Reconstructed/Associated Truth Information" << std::endl;
	std::cout << "#################################################" << std::endl;
	std::cout << "Hit Size: " << showerhits.size() << std::endl;
	std::cout << "ShowerBest_Plane: " << ShowerBest_Plane << std::endl;
	if(EvaluateShowerDirection){
	  std::cout << "True Start: " << PositionTrajStart.X() << " Shower Start: " << ShowerStart.X() << std::endl;
	  std::cout << "X Poisition: " <<  ShowerStart.X() << "Y Position " << ShowerStart.Y() << " Z Poistion: " << ShowerStart.Z() << std::endl;
	}
	if(EvaluateShowerLength){
	  std::cout << "TrueTrackLength: " << TrueTrackLength << " ShowerTrackLength: " << ShowerTrackLength << std::endl;
	}
	if(EvaluateShowerEnergy){
	  std::cout << "Best Plane Reco Shower Energy " << ShowerEnergyPlanes_remove[ShowerBest_Plane] << std::endl;
	}
	std::cout << "True Energy Deposited from true shower within Shower: " << TrueEnergyDepWithinShower_FromTrueShower << std::endl;
	std::cout << "True Energy Deposited From true associated Shower: " << TrueEnergyDep_FromShower << std::endl;
	std::cout << "True Energy Deposited by hits in shower: " <<  TrueEnergyDep_WithinRecoShower << std::endl;
	std::cout << "Energy Purity: " << energypurity << " Energy completeness: " << energycompleteness << std::endl;
	std::cout << "Hit Purity: " << hitpurity << "Hit Completeness: " << hitcompleteness << std::endl;
	std::cout << "#################################################" <<std::endl;
      }

      if(fVerbose > 1){std::cout << "Shower Validation Complete" << std::endl;} 
      

      //##########################
      //### Cluster Validation ###
      //##########################
      
      //Firstly look inf the hits in the shower match in all given planes
      std::map<std::vector<double>, int> HitCoord_map;
      float numclusters = -9999;
      float geoprojectionmatched_score = -99999;
      for(std::vector< art::Ptr<recob::Hit> >::iterator hitIt=showerhits.begin(); hitIt!=showerhits.end(); ++hitIt){
	try{
	  const std::vector<sim::IDE> hitIDEs = backtracker->HitToAvgSimIDEs(*hitIt);
	  for(unsigned int hitIDE=0; hitIDE<hitIDEs.size(); ++hitIDE){
	    std::vector<double> hitcoord = {hitIDEs[hitIDE].x, hitIDEs[hitIDE].y, hitIDEs[hitIDE].z}; 
	    //std::cout << "x: " << hitIDEs[hitIDE].x << ", y: " << hitIDEs[hitIDE].y << ", z: " << hitIDEs[hitIDE].z << " plane: " << (*hitIt)->WireID().Plane << std::endl;;
	    ++HitCoord_map[hitcoord];
	  }
	}
	catch(...){
	  if(fVerbose>1){std::cout << "Noise Hit" << std::endl;}
	}
      }
    
      //Get the clusters associated to the shower.
      if(fmch.isValid()){
	art::Handle<std::vector<recob::Cluster > > clusterHandle;
	evt.get(fmch.at(shower.key()).front().id(),clusterHandle);
	if(clusterHandle.isValid()){
	  std::vector<art::Ptr<recob::Cluster> > showerclusters = fmch.at(shower.key());
	  ana::ShowerValidation::ClusterValidation(showerclusters,evt,clusterHandle,ShowersMothers,MCTrack_Energy_map,MCTrack_hit_map,ShowerTrackID,simenergy,fShowerModuleLabel);


	  //Loop over the hit coordinate map where there there is a hit on every plane give a 1.
	  numclusters = showerclusters.size();
	  int geoprojectionmatched = 0;  
	  for(std::map<std::vector<double>, int>::iterator coord=HitCoord_map.begin(); coord!=HitCoord_map.end(); ++coord){
	    if(coord->second == numclusters){++geoprojectionmatched;}
	  }
	  if(numclusters > 0){ 
	    EvalulateGeoProjectionMatched = true;
	    geoprojectionmatched_score = (float) geoprojectionmatched/(float) HitCoord_map.size();
	  } 
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
	      
	      //Loop over the hit coordinate map where there there is a hit on every plane give a 1.
	      numclusters = showerclusters.size();
	      int geoprojectionmatched = 0;  
	      for(std::map<std::vector<double>, int>::iterator coord=HitCoord_map.begin(); coord!=HitCoord_map.end(); ++coord){
		if(coord->second == numclusters){++geoprojectionmatched;}
	      }

	      if(numclusters > 0){ 
		EvalulateGeoProjectionMatched = true;
		geoprojectionmatched_score = (float) geoprojectionmatched/(float) HitCoord_map.size();
	      }
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
      
      if(EvalulateGeoProjectionMatched){
	std::cout << "numclusters: " << numclusters << std::endl;
	std::cout << "geoprojectionmatched_score: " << geoprojectionmatched_score << std::endl;
	GeoProjectionMatched_HistMap[fShowerModuleLabel]->Fill(geoprojectionmatched_score);
	GeoProjectionMatched_2dHistMap[fShowerModuleLabel]->Fill(geoprojectionmatched_score,simenergy*1000);
	GeoProjectionMatched_TreeVal[fShowerModuleLabel].push_back(geoprojectionmatched_score);

	//Fill the Energy dependent Histograms. 
	if(fEnergies.size() != 0){
	  for(unsigned int i=0; i<fEnergies.size(); ++i){
	    if(TMath::Abs(simenergy*1000 - fEnergies[i])< fEnergyWidth){
	      Energies_GeoProjectionMatched_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(geoprojectionmatched_score);
	    }
	  }
	}
      }
      else{
	GeoProjectionMatched_TreeVal[fShowerModuleLabel].push_back(-99999);
      }
      
      if(fVerbose > 1){std::cout << "Cluster Validation Complete" << std::endl;}     

    }//Shower Loop 
   
    //#########################
    //###   Event Metric    ### 
    //#########################
    
    //Calculate the True Hit number 
    for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end(); ++showermother){
      int TrueHitDep_FromTrueShowers = 0;
      for(std::vector<int>::iterator track_id=(showermother->second).begin(); track_id!=(showermother->second).end(); ++track_id){
	for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){
	  TrueHitDep_FromTrueShowers +=  MCTrack_hit_map[showerhit_productid][*track_id][plane_id];
	}
      }
      TrueHitNum_HistMap[fShowerModuleLabel]  ->Fill(TrueHitDep_FromTrueShowers);
      TrueHitNum_2dHistMap[fShowerModuleLabel]->Fill(TrueHitDep_FromTrueShowers,simenergy*1000);
      
      if(fEnergies.size() != 0){
	for(unsigned int i=0; i<fEnergies.size(); ++i){
	  if(TMath::Abs(simenergy*1000 - fEnergies[i])< fEnergyWidth){
	    Energies_TrueHitNum_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(TrueHitDep_FromTrueShowers);
	  }
	}
      }
    }

    //Whats the segementyness of the event.
    EventSeggy_HistMap[fShowerModuleLabel]  ->Fill(showers.size()/num_of_showers_viaDensitycut);
    TrueEnergy_HistMap[fShowerModuleLabel]  ->Fill(simenergy*1000);
    EventSeggy_2dHistMap[fShowerModuleLabel]->Fill(showers.size()/num_of_showers_viaDensitycut,simenergy*1000);
    
    EventSeggy_TreeVal[fShowerModuleLabel].push_back(showers.size()/num_of_showers_viaDensitycut);
    TrueEnergy_TreeVal[fShowerModuleLabel].push_back(simenergy*1000);

    if(fEnergies.size() != 0){
      for(unsigned int i=0; i<fEnergies.size(); ++i){
	if(TMath::Abs(simenergy*1000 - fEnergies[i])< fEnergyWidth){
	  Energies_EventSeggy_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(showers.size()/num_of_showers_viaDensitycut);
	  Energies_TrueEnergy_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(simenergy*1000);
	}
      }
    }
  
    //As starting position shouldn't be penalised for segmentation then here is the option to fill the closest shower to the true start
    if(fFillOnlyClosestShower){
      for(std::map<int,float>::iterator trueshwr_iter=MinStartDiff.begin() ; trueshwr_iter!=MinStartDiff.end(); ++trueshwr_iter){
	if(MinEvaluateShowerDirection[trueshwr_iter->first] > 0){
	  ShowerDirection_X_HistMap[fShowerModuleLabel]   ->Fill(MinShowerDirection_Xdiff[trueshwr_iter->first]);
	  ShowerDirection_Y_HistMap[fShowerModuleLabel]   ->Fill(MinShowerDirection_Ydiff[trueshwr_iter->first]);
	  ShowerDirection_Z_HistMap[fShowerModuleLabel]   ->Fill(MinShowerDirection_Zdiff[trueshwr_iter->first]);
	  ShowerDirectionDiff_HistMap[fShowerModuleLabel] ->Fill(MinShowerDirection_diff[trueshwr_iter->first]);
	  ShowerDirection_X_2dHistMap[fShowerModuleLabel] ->Fill(MinShowerDirection_Xdiff[trueshwr_iter->first],simenergy*1000);
	  ShowerDirection_Y_2dHistMap[fShowerModuleLabel] ->Fill(MinShowerDirection_Ydiff[trueshwr_iter->first],simenergy*1000);
	  ShowerDirection_Z_2dHistMap[fShowerModuleLabel] ->Fill(MinShowerDirection_Zdiff[trueshwr_iter->first],simenergy*1000);
	  ShowerDirectionDiff_2dHistMap[fShowerModuleLabel]->Fill(MinShowerDirection_diff[trueshwr_iter->first],simenergy*1000);
	}
      }

      for(std::map<int,float>::iterator trueshwr_iter=MinStartDiff.begin() ; trueshwr_iter!=MinStartDiff.end(); ++trueshwr_iter){
	if(MinEvaluateShowerStart[trueshwr_iter->first] > 0){
	  ShowerStart_X_HistMap[fShowerModuleLabel]      ->Fill(MinStartX[trueshwr_iter->first]);
	  ShowerStart_Y_HistMap[fShowerModuleLabel]      ->Fill(MinStartY[trueshwr_iter->first]);
	  ShowerStart_Z_HistMap[fShowerModuleLabel]      ->Fill(MinStartZ[trueshwr_iter->first]);
	  ShowerMag_HistMap[fShowerModuleLabel]          ->Fill(MinStartDiff[trueshwr_iter->first]);
	  ShowerStart_X_2dHistMap[fShowerModuleLabel]    ->Fill(MinStartX[trueshwr_iter->first],simenergy*1000);
	  ShowerStart_Y_2dHistMap[fShowerModuleLabel]    ->Fill(MinStartY[trueshwr_iter->first],simenergy*1000);
	  ShowerStart_Z_2dHistMap[fShowerModuleLabel]    ->Fill(MinStartZ[trueshwr_iter->first],simenergy*1000);
	  ShowerMag_2dHistMap[fShowerModuleLabel]        ->Fill(MinStartDiff[trueshwr_iter->first],simenergy*1000);
	}	
      }
	
      if(fEnergies.size() != 0){
	for(unsigned int i=0; i<fEnergies.size(); ++i){
	  if(TMath::Abs(simenergy*1000 - fEnergies[i])< fEnergyWidth){
	    for(std::map<int,float>::iterator trueshwr_iter=MinStartDiff.begin() ; trueshwr_iter!=MinStartDiff.end(); ++trueshwr_iter){
	      if(MinEvaluateShowerDirection[trueshwr_iter->first] > 0){
		Energies_ShowerDirection_X_HistMap[fShowerModuleLabel][fEnergies[i]]  ->Fill(MinShowerDirection_Xdiff[trueshwr_iter->first]);
		Energies_ShowerDirection_Y_HistMap[fShowerModuleLabel][fEnergies[i]]  ->Fill(MinShowerDirection_Ydiff[trueshwr_iter->first]);
		Energies_ShowerDirection_Z_HistMap[fShowerModuleLabel][fEnergies[i]]  ->Fill(MinShowerDirection_Zdiff[trueshwr_iter->first]);
		Energies_ShowerDirectionDiff_HistMap[fShowerModuleLabel][fEnergies[i]]->Fill(MinShowerDirection_diff[trueshwr_iter->first]);
	      }
	    }

	    for(std::map<int,float>::iterator trueshwr_iter=MinStartDiff.begin() ; trueshwr_iter!=MinStartDiff.end(); ++trueshwr_iter){
	      if(MinEvaluateShowerStart[trueshwr_iter->first] > 0){
		Energies_ShowerStart_X_HistMap[fShowerModuleLabel][fEnergies[i]]      ->Fill(MinStartX[trueshwr_iter->first]);
		Energies_ShowerStart_Y_HistMap[fShowerModuleLabel][fEnergies[i]]      ->Fill(MinStartY[trueshwr_iter->first]);
		Energies_ShowerStart_Z_HistMap[fShowerModuleLabel][fEnergies[i]]      ->Fill(MinStartZ[trueshwr_iter->first]);
		Energies_ShowerMag_HistMap[fShowerModuleLabel][fEnergies[i]]          ->Fill(MinStartDiff[trueshwr_iter->first]);
	      }
	    }
	  }
	}
      }
    }
  }//Shower Module labels
  
  //Fill the tree
  Tree->Fill();

  for(unsigned int shwrlab_it=0; shwrlab_it<fShowerModuleLabels.size(); ++shwrlab_it){
    ShowerDirection_X_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerDirection_Y_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerDirection_Z_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerStart_X_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerStart_Y_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerStart_Z_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerLength_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerEnergyDiff_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerdEdx_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    EventSeggy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerEnergyCompleteness_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerHitsPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerHitsCompleteness_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerEnergyPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerEnergy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerHitNum_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerTotalEnergyDiff_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerMag_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerDirectionDiff_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerRecoEnergyVsTrueEnergyinRecoShower_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ShowerTrueEnergy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    TrueEnergy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    TrueHitNum_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    PosDir_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();

    ClusterProjectionMatchedEnergy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ClusterCompletenessEnergy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ClusterPurityEnergy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ClusterCompletenessHits_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ClusterPurityHits_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ClusterCompPurityEnergy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    ClusterCompPurityHits_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    
  }

  for(unsigned int hitlab_it=0; hitlab_it<fHitModuleLabels.size(); ++hitlab_it){
    for(unsigned int plane_it=0; plane_it<geom->Nplanes(); ++plane_it){
        HitCompletenessEnergy_TreeVal[fHitModuleLabels[hitlab_it]][plane_it].clear();
    }
  }

  SmallestShowerMotherE_BeforeCut_TreeVal.clear();
  BiggestShowerMotherE_BeforeCut_TreeVal.clear();
  SmallestShowerMotherE_AfterCut_TreeVal.clear();
  BiggestShowerMotherE_AfterCut_TreeVal.clear();


  return; 
}


void ana::ShowerValidation::ClusterValidation(std::vector< art::Ptr<recob::Cluster> >& clusters, 
					      const art::Event& evt, 
					      art::Handle<std::vector<recob::Cluster> >& clusterHandle, 
					      std::map<int,std::vector<int> >& ShowerMotherTrackIDs, 
					      std::map<int,float>& MCTrack_Energy_map, 
					      std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > >& MCTrack_hit_map, 
					      int& TrueShowerID, 
					      float& simenergy,
					      std::string& fShowerModuleLabel){

  //Initialise Trees                                                                                                                                                               
  for(unsigned int j=0; j<fShowerModuleLabels.size(); ++j){
    ClusterCompletenessEnergy_TreeVal     [fShowerModuleLabel].push_back(std::vector<float>(clusters.size(),-99999));
    //    ClusterPurityEnergy_TreeVal           [fShowerModuleLabel].push_back(std::vector<std::vector<float> >((int)geom->Nplanes()));
    //ClusterCompletenessHits_TreeVal       [fShowerModuleLabel].push_back(std::vector<std::vector<float> >((int)geom->Nplanes()));
    //ClusterPurityHits_TreeVal             [fShowerModuleLabel].push_back(std::vector<std::vector<float> >((int)geom->Nplanes()));
    //ClusterCompPurityEnergy_TreeVal       [fShowerModuleLabel].push_back(std::vector<std::vector<float> >((int)geom->Nplanes()));
    //ClusterCompPurityHits_TreeVal         [fShowerModuleLabel].push_back(std::vector<std::vector<float> >((int)geom->Nplanes()));
  }

  //Get the associated hits 
  art::FindManyP<recob::Hit> fmhc(clusterHandle, evt, clusterHandle.provenance()->moduleLabel());

  //Holder for cluster hits
  std::vector< art::Ptr<recob::Hit> > clusterhits; 

  //Get the Hits Handle used for this cluster type 
  art::Handle<std::vector<recob::Hit > > hitHandle;
  evt.get(fmhc.at(clusters.at(0).key()).front().id(),hitHandle);
  
  if(!hitHandle.isValid()){   
    mf::LogError("ShowerValidation") << "Hits handle is stale. No clustering validation done" << std::endl;
    return; 
  }

  //Get the hits vector from the shower
  for(auto const& cluster : clusters){

    clusterhits = fmhc.at(cluster.key());
    
    //Function from RecoUtils, finds the most probable track ID associated with the set of hits from there true energy depositons. The pair returns the energy in the hits as well.
    std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain(ShowerMotherTrackIDs,clusterhits, cluster->Plane().Plane);

    //Make sure the cluster has been matched. 
    if(ShowerTrackInfo.second == -99999){
      std::cout << "Reco cluster not matched to a True shower" << std::endl;
      continue;
    }

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
    float completeness_hits   = 0;
    float purity_hits         = 0;
    float completeness_energy = 0;
    float purity_energy       = 0;
    
    float TotalEnergyDepinHits = RecoUtils::TotalEnergyDepinHits(clusterhits,cluster->Plane().Plane);
    
    ClusterProjectionMatchedEnergy_HistMap[fShowerModuleLabel][cluster->Plane()]  ->Fill(projection_match);
    ClusterProjectionMatchedEnergy_2dHistMap[fShowerModuleLabel][cluster->Plane()]->Fill(projection_match,simenergy*1000);
    //    ClusterProjectionMatchedEnergy_TreeVal[fShowerModuleLabel][ClusterProjectionMatchedEnergy_TreeVal.size()-1][cluster->Plane().Plane].push_back(projection_match);

    if(totalhits != 0){
      completeness_hits = (signalhits)/totalhits;
      ClusterCompletenessHits_HistMap[fShowerModuleLabel][cluster->Plane()]  ->Fill(completeness_hits);
      ClusterCompletenessHits_2dHistMap[fShowerModuleLabel][cluster->Plane()]->Fill(completeness_hits,simenergy*1000);
      //      ClusterCompletenessHits_TreeVal[fShowerModuleLabel][ClusterCompletenessHits_TreeVal.size()-1][cluster->Plane().Plane].push_back(completeness_hits);
    }

    if(clusterhits.size() != 0){
      purity_hits = signalhits/clusterhits.size();
      ClusterPurityHits_HistMap[fShowerModuleLabel][cluster->Plane()]  ->Fill(purity_hits);
      ClusterPurityHits_2dHistMap[fShowerModuleLabel][cluster->Plane()]->Fill(purity_hits,simenergy*1000);
      //ClusterPurityHits_TreeVal[fShowerModuleLabel][ClusterPurityHits_TreeVal.size()-1][cluster->Plane().Plane].push_back(purity_hits);
    }

    if(TotalTrueEnergy != 0){
      completeness_energy = (ShowerTrackInfo.second)/TotalTrueEnergy; 
      ClusterCompletenessEnergy_HistMap[fShowerModuleLabel][cluster->Plane()]  ->Fill(completeness_energy);
      ClusterCompletenessEnergy_2dHistMap[fShowerModuleLabel][cluster->Plane()]->Fill(completeness_energy,simenergy*1000);
      //      ClusterCompletenessEnergy_TreeVal[fShowerModuleLabel][ClusterCompletenessEnergy_TreeVal.size()-1][cluster->Plane().Plane].push_back(completeness_energy);
    }

    if(TotalEnergyDepinHits != 0){
      purity_energy = ShowerTrackInfo.second/TotalEnergyDepinHits;
      ClusterPurityEnergy_HistMap[fShowerModuleLabel][cluster->Plane()]  ->Fill(purity_energy);
      ClusterPurityEnergy_2dHistMap[fShowerModuleLabel][cluster->Plane()]->Fill(purity_energy,simenergy*1000);
      //      ClusterPurityEnergy_TreeVal[fShowerModuleLabel][ClusterPurityEnergy_TreeVal.size()-1][cluster->Plane().Plane].push_back(purity_energy);
    }

    if(totalhits != 0 && clusterhits.size() != 0){
      ClusterCompPurityHits_HistMap[fShowerModuleLabel][cluster->Plane()]  ->Fill(completeness_hits*purity_hits);
      ClusterCompPurityHits_2dHistMap[fShowerModuleLabel][cluster->Plane()]->Fill(completeness_hits*purity_hits,simenergy*1000);
      //ClusterCompPurityHits_TreeVal[fShowerModuleLabel][ClusterCompPurityHits_TreeVal.size()-1][cluster->Plane().Plane].push_back(completeness_hits*purity_hits);
    }

    if(TotalTrueEnergy != 0 && TotalEnergyDepinHits != 0){
      ClusterCompPurityEnergy_HistMap[fShowerModuleLabel][cluster->Plane()]  ->Fill(completeness_energy*purity_energy);
      ClusterCompPurityEnergy_2dHistMap[fShowerModuleLabel][cluster->Plane()]->Fill(completeness_energy*purity_energy,simenergy*1000);
      // ClusterCompPurityEnergy_TreeVal[fShowerModuleLabel][ClusterCompPurityEnergy_TreeVal.size()-1][cluster->Plane().Plane].push_back(completeness_energy*purity_energy);
    }

    if(fVerbose>1){
      std::cout << "#################################################"   << std::endl;
      std::cout << "               Cluster Metrics                   "   << std::endl;
      std::cout << "#################################################"   << std::endl;
      std::cout << "Projection matched:          " << projection_match   << std::endl;
      std::cout << "Cluster hit completeness:    " << completeness_hits  << std::endl;
      std::cout << "Cluster hit purity:          " << purity_hits        << std::endl;
      std::cout << "Cluster energy completeness: " << completeness_energy << std::endl;
      std::cout << "CLuster energy purity:       " << purity_energy      << std::endl;
      std::cout << "#################################################"   << std::endl;
    }
    
    //Add to the Energy dependent histograms
    if(fEnergies.size() != 0){
      for(unsigned int i=0; i<fEnergies.size(); ++i){
        if(TMath::Abs(simenergy*1000 - fEnergies[i])< fEnergyWidth){
	  Energies_ClusterProjectionMatchedEnergy_HistMap[fShowerModuleLabel][fEnergies[i]][cluster->Plane()]->Fill(projection_match);
	  if(TotalTrueEnergy != 0){Energies_ClusterCompletenessEnergy_HistMap[fShowerModuleLabel][fEnergies[i]][cluster->Plane()]->Fill(completeness_energy);}
	  if(TotalEnergyDepinHits != 0){Energies_ClusterPurityEnergy_HistMap[fShowerModuleLabel][fEnergies[i]][cluster->Plane()] ->Fill(purity_energy);}
	  if(totalhits != 0){Energies_ClusterCompletenessHits_HistMap[fShowerModuleLabel][fEnergies[i]][cluster->Plane()]        ->Fill(completeness_hits);}
	  if(clusterhits.size() != 0){Energies_ClusterPurityHits_HistMap[fShowerModuleLabel][fEnergies[i]][cluster->Plane()]     ->Fill(purity_hits);}
	  if(TotalTrueEnergy != 0 && TotalEnergyDepinHits != 0){Energies_ClusterCompPurityEnergy_HistMap[fShowerModuleLabel][fEnergies[i]][cluster->Plane()]->Fill(completeness_energy*purity_energy);}
	  if(totalhits != 0 && clusterhits.size() != 0){Energies_ClusterCompPurityHits_HistMap[fShowerModuleLabel][fEnergies[i]][cluster->Plane()]         ->Fill(completeness_hits*purity_hits);}
	}
      }
    }

  }//Cluster Loop
  return;
}

void ana::ShowerValidation::DrawGraphs(std::map<std::string,TH1F*>& Name_HistMap, 
				       std::map<std::string,std::map<float,TH1F*> >& Energies_Name_HistMap, 
				       std::map<std::string,TGraphErrors*>& Energies_Mean_Name_GraphMap,
				       std::map<std::string,TGraphErrors*>& Energies_RMS_Name_GraphMap,
				       TMultiGraph*& Energies_Mean_Name_Multi,TMultiGraph*& Energies_RMS_Name_Multi,
				       TCanvas*& Energies_Mean_Name_canvasMulti,TCanvas*& Energies_RMS_Name_canvasMulti,
				       std::map<float,TCanvas*>& Energies_Name_canvasMap, TCanvas*& Name_canvas,
				       std::map<std::string,TH2F*>& Name_2dHistMap,
				       std::map<std::string,TCanvas*>& Name_2dCanvasMap
				       ){

  std::vector<int> colours = {600,632,800,880,840};
  if((int)fShowerModuleLabels.size() > 5){
    for(unsigned int j=0; j<fShowerModuleLabels.size() -5 ; ++j){
      colours.push_back(j);
    }
  }
  
  
  int Name_HistMap_Max = 0;
  for(unsigned int j=0; j<fShowerModuleLabels.size(); ++j){
    if(Name_HistMap[fShowerModuleLabels[j]]->GetMaximum() > Name_HistMap_Max){
      Name_HistMap_Max = Name_HistMap[fShowerModuleLabels[j]]->GetMaximum();
    }
  }
  

  for(unsigned int j=0; j<fShowerModuleLabels.size(); ++j){
    
    Name_canvas->cd();
    Name_HistMap[fShowerModuleLabels[j]]->SetLineColor(colours[j]);
    Name_HistMap[fShowerModuleLabels[j]]->SetLineWidth(2);
    Name_HistMap[fShowerModuleLabels[j]]->Draw("SAME");
    Name_HistMap[fShowerModuleLabels[j]]->SetMinimum(0);
    Name_HistMap[fShowerModuleLabels[j]]->SetMaximum(Name_HistMap_Max+1);
    Name_canvas->Update();

    Name_2dCanvasMap[fShowerModuleLabels[j]]->cd();
    Name_2dHistMap[fShowerModuleLabels[j]]->Draw("COLZ");

    if(fDrawCanvases){
      std::string canvasname(Name_2dCanvasMap[fShowerModuleLabels[j]]->GetName());
      std::string pngstring = canvasname + ".png";
      const char* pngname = pngstring.c_str();
      Name_2dCanvasMap[fShowerModuleLabels[j]]->Print(pngname);
    }


    Energies_RMS_Name_GraphMap[fShowerModuleLabels[j]] ->SetMarkerStyle(8);
    Energies_Mean_Name_GraphMap[fShowerModuleLabels[j]]->SetMarkerStyle(8);
    Energies_Mean_Name_GraphMap[fShowerModuleLabels[j]]->SetMarkerColor(colours[j]);    
    Energies_RMS_Name_GraphMap[fShowerModuleLabels[j]] ->SetMarkerColor(colours[j]);


    for(unsigned int i=0; i<fEnergies.size(); ++i){

      float Enteries = Energies_Name_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetEntries();
      
      if(Enteries == 0){continue;}
            
      Energies_Name_canvasMap[fEnergies[i]]->cd();
      Energies_Name_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetLineWidth(2);
      Energies_Name_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetLineColor(colours[j]);
      Energies_Name_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->Draw("SAME");
      Energies_Name_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetMinimum(0);
      Energies_Name_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->SetMaximum(Name_HistMap_Max+1);

      float Energies_TrueEnergy_Mean = Energies_TrueEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();
      float Energies_TrueEnergy_RMS = Energies_TrueEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
      float Error_on_True_Energy = Energies_TrueEnergy_RMS;
    
      float Energies_Name_Mean = Energies_Name_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();
      float Energies_Name_RMS = Energies_Name_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
      Energies_Mean_Name_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_Mean_Name_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_Name_Mean);
      Energies_RMS_Name_GraphMap[fShowerModuleLabels[j]]->SetPoint(Energies_RMS_Name_GraphMap[fShowerModuleLabels[j]]->GetN(),Energies_TrueEnergy_Mean,Energies_Name_RMS);
      Energies_Mean_Name_GraphMap[fShowerModuleLabels[j]]->SetPointError(Energies_Mean_Name_GraphMap[fShowerModuleLabels[j]]->GetN()-1,Error_on_True_Energy,Energies_Name_RMS);

    }//Energy Loop
    
    Energies_Mean_Name_Multi->Add(Energies_Mean_Name_GraphMap[fShowerModuleLabels[j]]);
    Energies_RMS_Name_Multi->Add(Energies_RMS_Name_GraphMap[fShowerModuleLabels[j]]);

  }//Shower Labels Loop


  Name_canvas->cd();
  Name_canvas->BuildLegend();
  if(fDrawCanvases){
    std::string canvasname(Name_canvas->GetName());
    std::string pngstring = canvasname + ".png";
    std::cout << "pngstring: " << pngstring << std::endl;
    const char* pngname = pngstring.c_str();
    Name_canvas->Print(pngname);
  }
  
  for(unsigned int i=0; i<fEnergies.size(); ++i){
    Energies_Name_canvasMap[fEnergies[i]]->cd();
    Energies_Name_canvasMap[fEnergies[i]]->BuildLegend();
  }

  Energies_Mean_Name_canvasMulti->cd();
  Energies_Mean_Name_Multi->Draw("AP");
  Energies_Mean_Name_canvasMulti->Update();
  Energies_Mean_Name_canvasMulti->BuildLegend(0.74, 0.7, 0.9, 0.9, "", "PL");
  
  if(fDrawCanvases){
    std::string canvasname(Energies_Mean_Name_canvasMulti->GetName());
    std::string pngstring = canvasname + ".png";
    const char* pngname = pngstring.c_str();
    Energies_Mean_Name_canvasMulti->Print(pngname);
  }

  Energies_RMS_Name_canvasMulti->cd();
  Energies_RMS_Name_Multi->Draw("AP");
  Energies_RMS_Name_canvasMulti->Update();
  Energies_RMS_Name_canvasMulti->BuildLegend(0.74, 0.7, 0.9, 0.9, "", "PL");

  if(fDrawCanvases){
    std::string canvasname(Energies_RMS_Name_canvasMulti->GetName());
    std::string pngstring = canvasname + ".png";
    const char* pngname = pngstring.c_str();
    Energies_RMS_Name_canvasMulti->Print(pngname);
  }
}


void ana::ShowerValidation::DrawGraphs(std::map<std::string,std::map<geo::PlaneID,TH1F*> >& Name_HistMap,
				       std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > >& Energies_Name_HistMap,
				       std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> >& Energies_Mean_Name_GraphMap,
				       std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> >& Energies_RMS_Name_GraphMap,
				       std::map<geo::PlaneID,TMultiGraph*>& Energies_Mean_Name_Multi,
				       std::map<geo::PlaneID,TMultiGraph*>& Energies_RMS_Name_Multi,
				       std::map<geo::PlaneID,TCanvas*>& Energies_Mean_Name_canvasMulti,
				       std::map<geo::PlaneID,TCanvas*>& Energies_RMS_Name_canvasMulti,
				       std::map<float,std::map<geo::PlaneID,TCanvas*> >& Energies_Name_canvasMap,
				       std::map<geo::PlaneID,TCanvas*>& Name_canvas,
				       std::map<std::string,std::map<geo::PlaneID,TH2F*> >& Name_2dHistMap,
				       std::map<std::string,std::map<geo::PlaneID,TCanvas*> >& Name_2dCanvasMap
				       ){


  std::vector<int> colours = {600,632,800,880,840};
  if((int)fShowerModuleLabels.size() > 5){
    for(unsigned int j=0; j<fShowerModuleLabels.size() -5 ; ++j){
      colours.push_back(j);
    }
  }

  for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){
    
    int Name_HistMap_Max = 0;
    for(unsigned int j=0; j<fShowerModuleLabels.size(); ++j){
      if(Name_HistMap[fShowerModuleLabels[j]][plane_id]->GetMaximum() > Name_HistMap_Max){
	Name_HistMap_Max = Name_HistMap[fShowerModuleLabels[j]][plane_id]->GetMaximum();
      }
    }

    
    for(unsigned int j=0; j<fShowerModuleLabels.size(); ++j){

      Name_canvas[plane_id]->cd();
      gStyle->SetPalette(kBird);
      Name_HistMap[fShowerModuleLabels[j]][plane_id]->SetLineWidth(2);
      Name_HistMap[fShowerModuleLabels[j]][plane_id]->SetLineColor(colours[j]);
      Name_HistMap[fShowerModuleLabels[j]][plane_id]->Draw("SAME");
      Name_HistMap[fShowerModuleLabels[j]][plane_id]->SetMinimum(0);
      Name_HistMap[fShowerModuleLabels[j]][plane_id]->SetMaximum(Name_HistMap_Max +1);
          
      Name_2dCanvasMap[fShowerModuleLabels[j]][plane_id]->cd();
      gStyle->SetPalette(kBird);
      Name_2dHistMap[fShowerModuleLabels[j]][plane_id]->Draw("COLZ");
      
      if(fDrawCanvases){
	std::string canvasname(Name_2dCanvasMap[fShowerModuleLabels[j]][plane_id]->GetName());
	std::string pngstring = canvasname + ".png";
        const char* pngname = pngstring.c_str();
        Name_2dCanvasMap[fShowerModuleLabels[j]][plane_id]->Print(pngname);
      }

      Energies_Mean_Name_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerStyle(8);
      Energies_RMS_Name_GraphMap[fShowerModuleLabels[j]][plane_id] ->SetMarkerStyle(8);
      Energies_Mean_Name_GraphMap[fShowerModuleLabels[j]][plane_id]->SetMarkerColor(colours[j]);    
      Energies_RMS_Name_GraphMap[fShowerModuleLabels[j]][plane_id] ->SetMarkerColor(colours[j]);
      
      for(unsigned int i=0; i<fEnergies.size(); ++i){
	
	float Enteries = Energies_Name_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->GetEntries();
	
	if(Enteries == 0){continue;}
	
	
	Energies_Name_canvasMap[fEnergies[i]][plane_id]->cd();
	gStyle->SetPalette(kBird);
	Energies_Name_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetLineWidth(2);
	Energies_Name_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetLineColor(colours[j]);
	Energies_Name_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->Draw("SAME");
	Energies_Name_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetMinimum(0);
        Energies_Name_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->SetMaximum(Name_HistMap_Max+1);
	
	float Energies_TrueEnergy_Mean = Energies_TrueEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetMean();   
	float Energies_TrueEnergy_RMS = Energies_TrueEnergy_HistMap[fShowerModuleLabels[j]][fEnergies[i]]->GetRMS();
	float Error_on_True_Energy = Energies_TrueEnergy_RMS;
	
	float Energies_Name_Mean = Energies_Name_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->GetMean();
	float Energies_Name_RMS = Energies_Name_HistMap[fShowerModuleLabels[j]][fEnergies[i]][plane_id]->GetRMS();
	Energies_Mean_Name_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPoint(Energies_Mean_Name_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN(),Energies_TrueEnergy_Mean,Energies_Name_Mean);
	Energies_RMS_Name_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPoint(Energies_RMS_Name_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN(),Energies_TrueEnergy_Mean,Energies_Name_RMS);
	Energies_Mean_Name_GraphMap[fShowerModuleLabels[j]][plane_id]->SetPointError(Energies_Mean_Name_GraphMap[fShowerModuleLabels[j]][plane_id]->GetN()-1,Error_on_True_Energy,Energies_Name_RMS);
	
      }//Energy Loop
      
      Energies_Mean_Name_Multi[plane_id]->Add(Energies_Mean_Name_GraphMap[fShowerModuleLabels[j]][plane_id]);
      Energies_RMS_Name_Multi[plane_id]->Add(Energies_RMS_Name_GraphMap[fShowerModuleLabels[j]][plane_id]);
    }//Shower Labels Loop

    Name_canvas[plane_id]->cd();
    Name_canvas[plane_id]->BuildLegend();
    if(fDrawCanvases){
      std::string canvasname(Name_canvas[plane_id]->GetName());
      std::string pngstring = canvasname + ".png";
      const char* pngname = pngstring.c_str();
      Name_canvas[plane_id]->Print(pngname);
    }

    for(unsigned int i=0; i<fEnergies.size(); ++i){
      Energies_Name_canvasMap[fEnergies[i]][plane_id]->cd();
      Energies_Name_canvasMap[fEnergies[i]][plane_id]->BuildLegend();
      if(fDrawCanvases){
	std::string canvasname(Energies_Name_canvasMap[fEnergies[i]][plane_id]->GetName());
	std::string pngstring = canvasname + ".png";
	const char* pngname = pngstring.c_str();
	Energies_Name_canvasMap[fEnergies[i]][plane_id]->Print(pngname);
      }
    }
    
    Energies_Mean_Name_canvasMulti[plane_id]->cd();
    Energies_Mean_Name_Multi[plane_id]->Draw("AP");
    Energies_Mean_Name_canvasMulti[plane_id]->Update();
    Energies_Mean_Name_canvasMulti[plane_id]->BuildLegend(0.74, 0.7, 0.9, 0.9, "", "PL");

    if(fDrawCanvases){
      std::string canvasname(Energies_Mean_Name_canvasMulti[plane_id]->GetName());
      std::string pngstring = canvasname + ".png";
      const char* pngname = pngstring.c_str();
      Energies_Mean_Name_canvasMulti[plane_id]->Print(pngname);
    }

    Energies_RMS_Name_canvasMulti[plane_id]->cd();
    Energies_RMS_Name_Multi[plane_id]->Draw("AP");
    Energies_RMS_Name_canvasMulti[plane_id]->Update();
    Energies_RMS_Name_canvasMulti[plane_id]->BuildLegend(0.74, 0.7, 0.9, 0.9, "", "PL");

    if(fDrawCanvases){
      std::string canvasname(Energies_RMS_Name_canvasMulti[plane_id]->GetName());
      std::string pngstring = canvasname + ".png";
      const char* pngname = pngstring.c_str();
      Energies_RMS_Name_canvasMulti[plane_id]->Print(pngname);
    }
  }
}

void ana::ShowerValidation::DrawHitGraphs(std::map<std::string,std::map<geo::PlaneID,TH1F*> >& Name_HistMap,
					  std::map<std::string,std::map<float,std::map<geo::PlaneID,TH1F*> > >& Energies_Name_HistMap,
					  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> >& Energies_Mean_Name_GraphMap,
					  std::map<std::string,std::map<geo::PlaneID,TGraphErrors*> >& Energies_RMS_Name_GraphMap,
					  std::map<geo::PlaneID,TMultiGraph*>& Energies_Mean_Name_Multi,
					  std::map<geo::PlaneID,TMultiGraph*>& Energies_RMS_Name_Multi,
					  std::map<geo::PlaneID,TCanvas*>& Energies_Mean_Name_canvasMulti,
					  std::map<geo::PlaneID,TCanvas*>& Energies_RMS_Name_canvasMulti,
					  std::map<float,std::map<geo::PlaneID,TCanvas*> >& Energies_Name_canvasMap,
					  std::map<geo::PlaneID,TCanvas*>& Name_canvas,
					  std::map<std::string,std::map<geo::PlaneID,TH2F*> >& Name_2dHistMap,
					  std::map<std::string,std::map<geo::PlaneID,TCanvas*> >& Name_2dCanvasMap
					  ){

  //Create the legend  
  TLegend *leg = new TLegend();
  for(unsigned int j=0; j<fHitModuleLabels.size(); ++j){
    const char* name_legend = fHitModuleLabels[j].c_str();
    leg->AddEntry(Name_HistMap[fHitModuleLabels[j]][geo::PlaneID(0,0,2)], name_legend);
  }


  std::vector<int> colours = {600,632,800,880,840};
  if((int)fHitModuleLabels.size() > 5){
    for(unsigned int j=0; j<fShowerModuleLabels.size() -5 ; ++j){
      colours.push_back(j);
    }
  }

  for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){
    
    int Name_HistMap_Max = 0;
    for(unsigned int j=0; j<fHitModuleLabels.size(); ++j){
      if(Name_HistMap[fHitModuleLabels[j]][plane_id]->GetMaximum() > Name_HistMap_Max){
	Name_HistMap_Max = Name_HistMap[fHitModuleLabels[j]][plane_id]->GetMaximum();
      }
    }

    for(unsigned int j=0; j<fHitModuleLabels.size(); ++j){
      
      Name_canvas[plane_id]->cd();
      Name_HistMap[fHitModuleLabels[j]][plane_id]->SetLineWidth(2);
      Name_HistMap[fHitModuleLabels[j]][plane_id]->SetLineColor(colours[j]);
      Name_HistMap[fHitModuleLabels[j]][plane_id]->Draw("SAME");
      Name_HistMap[fHitModuleLabels[j]][plane_id]->SetMinimum(0);
      Name_HistMap[fHitModuleLabels[j]][plane_id]->SetMaximum(Name_HistMap_Max+1); 

      Name_2dCanvasMap[fHitModuleLabels[j]][plane_id]->cd();
      gStyle->SetPalette(kBird);
      Name_2dHistMap[fHitModuleLabels[j]][plane_id]->Draw("COLZ");

      if(fDrawCanvases){
	std::string canvasname(Name_2dCanvasMap[fHitModuleLabels[j]][plane_id]->GetName());
	std::string pngstring = canvasname + ".png";
        const char* pngname = pngstring.c_str();
        Name_2dCanvasMap[fHitModuleLabels[j]][plane_id]->Print(pngname);
      }

      Energies_Mean_Name_GraphMap[fHitModuleLabels[j]][plane_id]->SetMarkerStyle(8);
      Energies_RMS_Name_GraphMap[fHitModuleLabels[j]][plane_id] ->SetMarkerStyle(8);
      Energies_Mean_Name_GraphMap[fHitModuleLabels[j]][plane_id]->SetMarkerColor(colours[j]);    
      Energies_RMS_Name_GraphMap[fHitModuleLabels[j]][plane_id] ->SetMarkerColor(colours[j]);
      
      for(unsigned int i=0; i<fEnergies.size(); ++i){
	
	float Enteries = Energies_Name_HistMap[fHitModuleLabels[j]][fEnergies[i]][plane_id]->GetEntries();
	
	if(Enteries == 0){continue;}
	
	Energies_Name_canvasMap[fEnergies[i]][plane_id]->cd();
	Energies_Name_HistMap[fHitModuleLabels[j]][fEnergies[i]][plane_id]->SetLineWidth(2);
	Energies_Name_HistMap[fHitModuleLabels[j]][fEnergies[i]][plane_id]->SetLineColor(colours[j]);
	Energies_Name_HistMap[fHitModuleLabels[j]][fEnergies[i]][plane_id]->Draw("SAME");
	Energies_Name_HistMap[fHitModuleLabels[j]][fEnergies[i]][plane_id]->SetMinimum(0);
	Energies_Name_HistMap[fHitModuleLabels[j]][fEnergies[i]][plane_id]->SetMaximum(Name_HistMap_Max +1);

	float Energies_TrueEnergy_Mean = Energies_TrueEnergy_HistMap[fShowerModuleLabels[0]][fEnergies[i]]->GetMean();   
	float Energies_TrueEnergy_RMS = Energies_TrueEnergy_HistMap[fShowerModuleLabels[0]][fEnergies[i]]->GetRMS();
	float Error_on_True_Energy = Energies_TrueEnergy_RMS;

	float Energies_Name_Mean = Energies_Name_HistMap[fHitModuleLabels[j]][fEnergies[i]][plane_id]->GetMean();
	float Energies_Name_RMS = Energies_Name_HistMap[fHitModuleLabels[j]][fEnergies[i]][plane_id]->GetRMS();
	Energies_Mean_Name_GraphMap[fHitModuleLabels[j]][plane_id]->SetPoint(Energies_Mean_Name_GraphMap[fHitModuleLabels[j]][plane_id]->GetN(),Energies_TrueEnergy_Mean,Energies_Name_Mean);
	Energies_RMS_Name_GraphMap[fHitModuleLabels[j]][plane_id]->SetPoint(Energies_RMS_Name_GraphMap[fHitModuleLabels[j]][plane_id]->GetN(),Energies_TrueEnergy_Mean,Energies_Name_RMS);
	Energies_Mean_Name_GraphMap[fHitModuleLabels[j]][plane_id]->SetPointError(Energies_Mean_Name_GraphMap[fHitModuleLabels[j]][plane_id]->GetN()-1,Error_on_True_Energy,Energies_Name_RMS);
	
      }//Energy Loop

      Energies_Mean_Name_Multi[plane_id]->Add(Energies_Mean_Name_GraphMap[fHitModuleLabels[j]][plane_id]);
      Energies_RMS_Name_Multi[plane_id]->Add(Energies_RMS_Name_GraphMap[fHitModuleLabels[j]][plane_id]);
    }//Shower Labels Loop

    Name_canvas[plane_id]->cd();
    Name_canvas[plane_id]->BuildLegend();

    if(fDrawCanvases){
	std::string canvasname(Name_canvas[plane_id]->GetName());
	std::string pngstring = canvasname + ".png";
        const char* pngname = pngstring.c_str();
        Name_canvas[plane_id]->Print(pngname);
      }

    for(unsigned int i=0; i<fEnergies.size(); ++i){
      Energies_Name_canvasMap[fEnergies[i]][plane_id]->cd();
      Energies_Name_canvasMap[fEnergies[i]][plane_id]->BuildLegend();
      if(fDrawCanvases){
	std::string canvasname(Energies_Name_canvasMap[fEnergies[i]][plane_id]->GetName());
	std::string pngstring = canvasname + ".png";
        const char* pngname = pngstring.c_str();
        Energies_Name_canvasMap[fEnergies[i]][plane_id]->Print(pngname);
      }
    }
    
    Energies_Mean_Name_canvasMulti[plane_id]->cd();
    Energies_Mean_Name_Multi[plane_id]->Draw("AP");
    Energies_Mean_Name_canvasMulti[plane_id]->Update();
    Energies_Mean_Name_canvasMulti[plane_id]->BuildLegend(0.74, 0.7, 0.9, 0.9, "", "PL");

    if(fDrawCanvases){
      std::string canvasname(Energies_Mean_Name_canvasMulti[plane_id]->GetName());
      std::string pngstring = canvasname + ".png";
      const char* pngname = pngstring.c_str();
      Energies_Mean_Name_canvasMulti[plane_id]->Print(pngname);
    }

    Energies_RMS_Name_canvasMulti[plane_id]->cd();
    Energies_RMS_Name_Multi[plane_id]->Draw("AP");
    Energies_RMS_Name_canvasMulti[plane_id]->Update();
    Energies_RMS_Name_canvasMulti[plane_id]->BuildLegend();
    Energies_RMS_Name_canvasMulti[plane_id]->BuildLegend(0.74, 0.7, 0.9, 0.9, "", "PL");

    if(fDrawCanvases){
      std::string canvasname(Energies_RMS_Name_canvasMulti[plane_id]->GetName());
      std::string pngstring = canvasname + ".png";
      const char* pngname = pngstring.c_str();
      Energies_RMS_Name_canvasMulti[plane_id]->Print(pngname);
    }
  }
}



void ana::ShowerValidation::endJob() {
  
  gStyle->SetOptTitle(0);

  for(unsigned int j=0; j<fShowerModuleLabels.size(); ++j){
    PosDir_2dCanvasMap[fShowerModuleLabels[j]]->cd();
    gStyle->SetPalette(kBird);
    PosDir_2dHistMap[fShowerModuleLabels[j]]->Draw("COLZ");
    
    if(fDrawCanvases){
      std::string canvasname(PosDir_2dCanvasMap[fShowerModuleLabels[j]]->GetName());
      std::string pngstring = canvasname + ".png";
      const char* pngname = pngstring.c_str();
      PosDir_2dCanvasMap[fShowerModuleLabels[j]]->Print(pngname);
    }
  }
    

  BiggestShowerMotherE_AfterCutCanvas->cd();
  BiggestShowerMotherE_AfterCutHist->Draw();
  if(fDrawCanvases){BiggestShowerMotherE_AfterCutCanvas->Print("BiggestShowerMotherE_AfterCutCanvas.png");}

  BiggestShowerMotherE_BeforeCutCanvas->cd();
  BiggestShowerMotherE_BeforeCutHist->Draw();
  if(fDrawCanvases){BiggestShowerMotherE_BeforeCutCanvas->Print("BiggestShowerMotherE_BeforeCutHist.png");}

  AsscoiatedBiggestShowerMotherECanvas->cd();
  AsscoiatedBiggestShowerMotherEHist->Draw();
  if(fDrawCanvases){AsscoiatedBiggestShowerMotherECanvas->Print("AsscoiatedBiggestShowerMotherEHist.png");}

  SmallestShowerMotherE_AfterCutCanvas->cd();
  SmallestShowerMotherE_AfterCutHist->Draw();
  if(fDrawCanvases){SmallestShowerMotherE_AfterCutCanvas->Print("SmallestShowerMotherE_AfterCutHist.png");}
  
  SmallestShowerMotherE_BeforeCutCanvas->cd();
  SmallestShowerMotherE_BeforeCutHist->Draw();
  if(fDrawCanvases){SmallestShowerMotherE_BeforeCutCanvas->Print("SmallestShowerMotherE_BeforeCutHist.png");}

  AsscoiatedSmallestShowerMotherECanvas->cd();
  AsscoiatedSmallestShowerMotherEHist->Draw();
  if(fDrawCanvases){AsscoiatedSmallestShowerMotherECanvas->Print("AsscoiatedSmallestShowerMotherEHist.png");}



  gStyle->SetPalette(kBird);

  ana::ShowerValidation::DrawGraphs(ShowerDirection_X_HistMap,Energies_ShowerDirection_X_HistMap,Energies_Mean_ShowerDirection_X_GraphMap,Energies_RMS_ShowerDirection_X_GraphMap,Energies_Mean_ShowerDirection_X_Multi,Energies_RMS_ShowerDirection_X_Multi,Energies_Mean_ShowerDirection_X_canvasMulti,Energies_RMS_ShowerDirection_X_canvasMulti,Energies_ShowerDirection_X_canvasMap, ShowerDirection_X_canvas,ShowerDirection_X_2dHistMap,ShowerDirection_X_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ShowerDirection_Y_HistMap,Energies_ShowerDirection_Y_HistMap,Energies_Mean_ShowerDirection_Y_GraphMap,Energies_RMS_ShowerDirection_Y_GraphMap,Energies_Mean_ShowerDirection_Y_Multi,Energies_RMS_ShowerDirection_Y_Multi,Energies_Mean_ShowerDirection_Y_canvasMulti,Energies_RMS_ShowerDirection_Y_canvasMulti,Energies_ShowerDirection_Y_canvasMap, ShowerDirection_Y_canvas,ShowerDirection_Y_2dHistMap,ShowerDirection_Y_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ShowerDirection_Z_HistMap,Energies_ShowerDirection_Z_HistMap,Energies_Mean_ShowerDirection_Z_GraphMap,Energies_RMS_ShowerDirection_Z_GraphMap,Energies_Mean_ShowerDirection_Z_Multi,Energies_RMS_ShowerDirection_Z_Multi,Energies_Mean_ShowerDirection_Z_canvasMulti,Energies_RMS_ShowerDirection_Z_canvasMulti,Energies_ShowerDirection_Z_canvasMap, ShowerDirection_Z_canvas,ShowerDirection_Z_2dHistMap,ShowerDirection_Z_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ShowerStart_X_HistMap,Energies_ShowerStart_X_HistMap,Energies_Mean_ShowerStart_X_GraphMap,Energies_RMS_ShowerStart_X_GraphMap,Energies_Mean_ShowerStart_X_Multi,Energies_RMS_ShowerStart_X_Multi,Energies_Mean_ShowerStart_X_canvasMulti,Energies_RMS_ShowerStart_X_canvasMulti,Energies_ShowerStart_X_canvasMap, ShowerStart_X_canvas,ShowerStart_X_2dHistMap,ShowerStart_X_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ShowerStart_Y_HistMap,Energies_ShowerStart_Y_HistMap,Energies_Mean_ShowerStart_Y_GraphMap,Energies_RMS_ShowerStart_Y_GraphMap,Energies_Mean_ShowerStart_Y_Multi,Energies_RMS_ShowerStart_Y_Multi,Energies_Mean_ShowerStart_Y_canvasMulti,Energies_RMS_ShowerStart_Y_canvasMulti,Energies_ShowerStart_Y_canvasMap, ShowerStart_Y_canvas,ShowerStart_Y_2dHistMap,ShowerStart_Y_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ShowerStart_Z_HistMap,Energies_ShowerStart_Z_HistMap,Energies_Mean_ShowerStart_Z_GraphMap,Energies_RMS_ShowerStart_Z_GraphMap,Energies_Mean_ShowerStart_Z_Multi,Energies_RMS_ShowerStart_Z_Multi,Energies_Mean_ShowerStart_Z_canvasMulti,Energies_RMS_ShowerStart_Z_canvasMulti,Energies_ShowerStart_Z_canvasMap, ShowerStart_Z_canvas,ShowerStart_Z_2dHistMap,ShowerStart_Z_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ShowerLength_HistMap,Energies_ShowerLength_HistMap,Energies_Mean_ShowerLength_GraphMap,Energies_RMS_ShowerLength_GraphMap,Energies_Mean_ShowerLength_Multi,Energies_RMS_ShowerLength_Multi,Energies_Mean_ShowerLength_canvasMulti,Energies_RMS_ShowerLength_canvasMulti,Energies_ShowerLength_canvasMap, ShowerLength_canvas,ShowerLength_2dHistMap,ShowerLength_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ShowerEnergyDiff_HistMap,Energies_ShowerEnergyDiff_HistMap,Energies_Mean_ShowerEnergyDiff_GraphMap,Energies_RMS_ShowerEnergyDiff_GraphMap,Energies_Mean_ShowerEnergyDiff_Multi,Energies_RMS_ShowerEnergyDiff_Multi,Energies_Mean_ShowerEnergyDiff_canvasMulti,Energies_RMS_ShowerEnergyDiff_canvasMulti,Energies_ShowerEnergyDiff_canvasMap, ShowerEnergyDiff_canvas,ShowerEnergyDiff_2dHistMap,ShowerEnergyDiff_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ShowerdEdx_HistMap,Energies_ShowerdEdx_HistMap,Energies_Mean_ShowerdEdx_GraphMap,Energies_RMS_ShowerdEdx_GraphMap,Energies_Mean_ShowerdEdx_Multi,Energies_RMS_ShowerdEdx_Multi,Energies_Mean_ShowerdEdx_canvasMulti,Energies_RMS_ShowerdEdx_canvasMulti,Energies_ShowerdEdx_canvasMap, ShowerdEdx_canvas,ShowerdEdx_2dHistMap,ShowerdEdx_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(EventSeggy_HistMap,Energies_EventSeggy_HistMap,Energies_Mean_EventSeggy_GraphMap,Energies_RMS_EventSeggy_GraphMap,Energies_Mean_EventSeggy_Multi,Energies_RMS_EventSeggy_Multi,Energies_Mean_EventSeggy_canvasMulti,Energies_RMS_EventSeggy_canvasMulti,Energies_EventSeggy_canvasMap, EventSeggy_canvas,EventSeggy_2dHistMap,EventSeggy_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ShowerEnergyCompleteness_HistMap,Energies_ShowerEnergyCompleteness_HistMap,Energies_Mean_ShowerEnergyCompleteness_GraphMap,Energies_RMS_ShowerEnergyCompleteness_GraphMap,Energies_Mean_ShowerEnergyCompleteness_Multi,Energies_RMS_ShowerEnergyCompleteness_Multi,Energies_Mean_ShowerEnergyCompleteness_canvasMulti,Energies_RMS_ShowerEnergyCompleteness_canvasMulti,Energies_ShowerEnergyCompleteness_canvasMap, ShowerEnergyCompleteness_canvas,ShowerEnergyCompleteness_2dHistMap,ShowerEnergyCompleteness_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ShowerEnergyPurity_HistMap,Energies_ShowerEnergyPurity_HistMap,Energies_Mean_ShowerEnergyPurity_GraphMap,Energies_RMS_ShowerEnergyPurity_GraphMap,Energies_Mean_ShowerEnergyPurity_Multi,Energies_RMS_ShowerEnergyPurity_Multi,Energies_Mean_ShowerEnergyPurity_canvasMulti,Energies_RMS_ShowerEnergyPurity_canvasMulti,Energies_ShowerEnergyPurity_canvasMap, ShowerEnergyPurity_canvas,ShowerEnergyPurity_2dHistMap,ShowerEnergyPurity_2dCanvasMap);

    ana::ShowerValidation::DrawGraphs(ShowerHitsCompleteness_HistMap,Energies_ShowerHitsCompleteness_HistMap,Energies_Mean_ShowerHitsCompleteness_GraphMap,Energies_RMS_ShowerHitsCompleteness_GraphMap,Energies_Mean_ShowerHitsCompleteness_Multi,Energies_RMS_ShowerHitsCompleteness_Multi,Energies_Mean_ShowerHitsCompleteness_canvasMulti,Energies_RMS_ShowerHitsCompleteness_canvasMulti,Energies_ShowerHitsCompleteness_canvasMap, ShowerHitsCompleteness_canvas,ShowerHitsCompleteness_2dHistMap,ShowerHitsCompleteness_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ShowerHitsPurity_HistMap,Energies_ShowerHitsPurity_HistMap,Energies_Mean_ShowerHitsPurity_GraphMap,Energies_RMS_ShowerHitsPurity_GraphMap,Energies_Mean_ShowerHitsPurity_Multi,Energies_RMS_ShowerHitsPurity_Multi,Energies_Mean_ShowerHitsPurity_canvasMulti,Energies_RMS_ShowerHitsPurity_canvasMulti,Energies_ShowerHitsPurity_canvasMap, ShowerHitsPurity_canvas,ShowerHitsPurity_2dHistMap,ShowerHitsPurity_2dCanvasMap);


  ana::ShowerValidation::DrawGraphs(ShowerEnergy_HistMap,Energies_ShowerEnergy_HistMap,Energies_Mean_ShowerEnergy_GraphMap,Energies_RMS_ShowerEnergy_GraphMap,Energies_Mean_ShowerEnergy_Multi,Energies_RMS_ShowerEnergy_Multi,Energies_Mean_ShowerEnergy_canvasMulti,Energies_RMS_ShowerEnergy_canvasMulti,Energies_ShowerEnergy_canvasMap, ShowerEnergy_canvas,ShowerEnergy_2dHistMap,ShowerEnergy_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ShowerHitNum_HistMap,Energies_ShowerHitNum_HistMap,Energies_Mean_ShowerHitNum_GraphMap,Energies_RMS_ShowerHitNum_GraphMap,Energies_Mean_ShowerHitNum_Multi,Energies_RMS_ShowerHitNum_Multi,Energies_Mean_ShowerHitNum_canvasMulti,Energies_RMS_ShowerHitNum_canvasMulti,Energies_ShowerHitNum_canvasMap, ShowerHitNum_canvas,ShowerHitNum_2dHistMap,ShowerHitNum_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(GeoProjectionMatched_HistMap,Energies_GeoProjectionMatched_HistMap,Energies_Mean_GeoProjectionMatched_GraphMap,Energies_RMS_GeoProjectionMatched_GraphMap,Energies_Mean_GeoProjectionMatched_Multi,Energies_RMS_GeoProjectionMatched_Multi,Energies_Mean_GeoProjectionMatched_canvasMulti,Energies_RMS_GeoProjectionMatched_canvasMulti,Energies_GeoProjectionMatched_canvasMap, GeoProjectionMatched_canvas,GeoProjectionMatched_2dHistMap,GeoProjectionMatched_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ShowerTotalEnergyDiff_HistMap,Energies_ShowerTotalEnergyDiff_HistMap,Energies_Mean_ShowerTotalEnergyDiff_GraphMap,Energies_RMS_ShowerTotalEnergyDiff_GraphMap,Energies_Mean_ShowerTotalEnergyDiff_Multi,Energies_RMS_ShowerTotalEnergyDiff_Multi,Energies_Mean_ShowerTotalEnergyDiff_canvasMulti,Energies_RMS_ShowerTotalEnergyDiff_canvasMulti,Energies_ShowerTotalEnergyDiff_canvasMap,ShowerTotalEnergyDiff_canvas,ShowerTotalEnergyDiff_2dHistMap,ShowerTotalEnergyDiff_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ShowerMag_HistMap,Energies_ShowerMag_HistMap,Energies_Mean_ShowerMag_GraphMap,Energies_RMS_ShowerMag_GraphMap,Energies_Mean_ShowerMag_Multi,Energies_RMS_ShowerMag_Multi,Energies_Mean_ShowerMag_canvasMulti,Energies_RMS_ShowerMag_canvasMulti,Energies_ShowerMag_canvasMap, ShowerMag_canvas,ShowerMag_2dHistMap,ShowerMag_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ShowerDirectionDiff_HistMap,Energies_ShowerDirectionDiff_HistMap,Energies_Mean_ShowerDirectionDiff_GraphMap,Energies_RMS_ShowerDirectionDiff_GraphMap,Energies_Mean_ShowerDirectionDiff_Multi,Energies_RMS_ShowerDirectionDiff_Multi,Energies_Mean_ShowerDirectionDiff_canvasMulti,Energies_RMS_ShowerDirectionDiff_canvasMulti,Energies_ShowerDirectionDiff_canvasMap, ShowerDirectionDiff_canvas,ShowerDirectionDiff_2dHistMap,ShowerDirectionDiff_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap,Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_HistMap,Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap,Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_GraphMap,Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_Multi,Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_Multi,Energies_Mean_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMulti,Energies_RMS_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMulti,Energies_ShowerRecoEnergyVsTrueEnergyinRecoShower_canvasMap, ShowerRecoEnergyVsTrueEnergyinRecoShower_canvas,ShowerRecoEnergyVsTrueEnergyinRecoShower_2dHistMap,ShowerRecoEnergyVsTrueEnergyinRecoShower_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ShowerTrueEnergy_HistMap,Energies_ShowerTrueEnergy_HistMap,Energies_Mean_ShowerTrueEnergy_GraphMap,Energies_RMS_ShowerTrueEnergy_GraphMap,Energies_Mean_ShowerTrueEnergy_Multi,Energies_RMS_ShowerTrueEnergy_Multi,Energies_Mean_ShowerTrueEnergy_canvasMulti,Energies_RMS_ShowerTrueEnergy_canvasMulti,Energies_ShowerTrueEnergy_canvasMap, ShowerTrueEnergy_canvas,ShowerTrueEnergy_2dHistMap,ShowerTrueEnergy_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(TrueEnergy_HistMap,Energies_TrueEnergy_HistMap,Energies_Mean_TrueEnergy_GraphMap,Energies_RMS_TrueEnergy_GraphMap,Energies_Mean_TrueEnergy_Multi,Energies_RMS_TrueEnergy_Multi,Energies_Mean_TrueEnergy_canvasMulti,Energies_RMS_TrueEnergy_canvasMulti,Energies_TrueEnergy_canvasMap, TrueEnergy_canvas,TrueEnergy_2dHistMap,TrueEnergy_2dCanvasMap);


  ana::ShowerValidation::DrawGraphs(TrueHitNum_HistMap,Energies_TrueHitNum_HistMap,Energies_Mean_TrueHitNum_GraphMap,Energies_RMS_TrueHitNum_GraphMap,Energies_Mean_TrueHitNum_Multi,Energies_RMS_TrueHitNum_Multi,Energies_Mean_TrueHitNum_canvasMulti,Energies_RMS_TrueHitNum_canvasMulti,Energies_TrueHitNum_canvasMap, TrueHitNum_canvas,TrueHitNum_2dHistMap,TrueHitNum_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ClusterProjectionMatchedEnergy_HistMap,Energies_ClusterProjectionMatchedEnergy_HistMap,Energies_Mean_ClusterProjectionMatchedEnergy_GraphMap,Energies_RMS_ClusterProjectionMatchedEnergy_GraphMap,Energies_Mean_ClusterProjectionMatchedEnergy_Multi,Energies_RMS_ClusterProjectionMatchedEnergy_Multi,Energies_Mean_ClusterProjectionMatchedEnergy_canvasMulti,Energies_RMS_ClusterProjectionMatchedEnergy_canvasMulti,Energies_ClusterProjectionMatchedEnergy_canvasMap, ClusterProjectionMatchedEnergy_canvas,ClusterProjectionMatchedEnergy_2dHistMap,ClusterProjectionMatchedEnergy_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ClusterCompletenessEnergy_HistMap,Energies_ClusterCompletenessEnergy_HistMap,Energies_Mean_ClusterCompletenessEnergy_GraphMap,Energies_RMS_ClusterCompletenessEnergy_GraphMap,Energies_Mean_ClusterCompletenessEnergy_Multi,Energies_RMS_ClusterCompletenessEnergy_Multi,Energies_Mean_ClusterCompletenessEnergy_canvasMulti,Energies_RMS_ClusterCompletenessEnergy_canvasMulti,Energies_ClusterCompletenessEnergy_canvasMap, ClusterCompletenessEnergy_canvas,ClusterCompletenessEnergy_2dHistMap,ClusterCompletenessEnergy_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ClusterCompletenessHits_HistMap,Energies_ClusterCompletenessHits_HistMap,Energies_Mean_ClusterCompletenessHits_GraphMap,Energies_RMS_ClusterCompletenessHits_GraphMap,Energies_Mean_ClusterCompletenessHits_Multi,Energies_RMS_ClusterCompletenessHits_Multi,Energies_Mean_ClusterCompletenessHits_canvasMulti,Energies_RMS_ClusterCompletenessHits_canvasMulti,Energies_ClusterCompletenessHits_canvasMap, ClusterCompletenessHits_canvas,ClusterCompletenessHits_2dHistMap,ClusterCompletenessHits_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ClusterPurityEnergy_HistMap,Energies_ClusterPurityEnergy_HistMap,Energies_Mean_ClusterPurityEnergy_GraphMap,Energies_RMS_ClusterPurityEnergy_GraphMap,Energies_Mean_ClusterPurityEnergy_Multi,Energies_RMS_ClusterPurityEnergy_Multi,Energies_Mean_ClusterPurityEnergy_canvasMulti,Energies_RMS_ClusterPurityEnergy_canvasMulti,Energies_ClusterPurityEnergy_canvasMap, ClusterPurityEnergy_canvas,ClusterPurityEnergy_2dHistMap,ClusterPurityEnergy_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ClusterPurityHits_HistMap,Energies_ClusterPurityHits_HistMap,Energies_Mean_ClusterPurityHits_GraphMap,Energies_RMS_ClusterPurityHits_GraphMap,Energies_Mean_ClusterPurityHits_Multi,Energies_RMS_ClusterPurityHits_Multi,Energies_Mean_ClusterPurityHits_canvasMulti,Energies_RMS_ClusterPurityHits_canvasMulti,Energies_ClusterPurityHits_canvasMap, ClusterPurityHits_canvas,ClusterPurityHits_2dHistMap,ClusterPurityHits_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ClusterCompPurityEnergy_HistMap,Energies_ClusterCompPurityEnergy_HistMap,Energies_Mean_ClusterCompPurityEnergy_GraphMap,Energies_RMS_ClusterCompPurityEnergy_GraphMap,Energies_Mean_ClusterCompPurityEnergy_Multi,Energies_RMS_ClusterCompPurityEnergy_Multi,Energies_Mean_ClusterCompPurityEnergy_canvasMulti,Energies_RMS_ClusterCompPurityEnergy_canvasMulti,Energies_ClusterCompPurityEnergy_canvasMap, ClusterCompPurityEnergy_canvas,ClusterCompPurityEnergy_2dHistMap,ClusterCompPurityEnergy_2dCanvasMap);

  ana::ShowerValidation::DrawGraphs(ClusterCompPurityHits_HistMap,Energies_ClusterCompPurityHits_HistMap,Energies_Mean_ClusterCompPurityHits_GraphMap,Energies_RMS_ClusterCompPurityHits_GraphMap,Energies_Mean_ClusterCompPurityHits_Multi,Energies_RMS_ClusterCompPurityHits_Multi,Energies_Mean_ClusterCompPurityHits_canvasMulti,Energies_RMS_ClusterCompPurityHits_canvasMulti,Energies_ClusterCompPurityHits_canvasMap, ClusterCompPurityHits_canvas,ClusterCompPurityHits_2dHistMap,ClusterCompPurityHits_2dCanvasMap);

  ana::ShowerValidation::DrawHitGraphs(HitCompletenessEnergy_HistMap,Energies_HitCompletenessEnergy_HistMap,Energies_Mean_HitCompletenessEnergy_GraphMap,Energies_RMS_HitCompletenessEnergy_GraphMap,Energies_Mean_HitCompletenessEnergy_Multi,Energies_RMS_HitCompletenessEnergy_Multi,Energies_Mean_HitCompletenessEnergy_canvasMulti,Energies_RMS_HitCompletenessEnergy_canvasMulti,Energies_HitCompletenessEnergy_canvasMap, HitCompletenessEnergy_canvas,HitCompletenessEnergy_2dHistMap,HitCompletenessEnergy_2dCanvasMap);

}


DEFINE_ART_MODULE(ana::ShowerValidation)
