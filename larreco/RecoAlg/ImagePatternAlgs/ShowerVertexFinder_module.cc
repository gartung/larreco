/////////////////////////////////////////////////////////////////////////////////
// Class:       ShowerVertexFinder
// Module Type: producer
// File:        ShowerVertexFinder_module.cc
// Authors:     L. Whitehead (leigh.howard.whitehead@cern.ch)
//
// Module applies CNN to 2D image made of deconvoluted wire waveforms in order
// to tag hits that are likely to be the starting point of showers. The CNN 
// returns a single value between 0 and 1, where values close to 1 represent 
// hits most likely to be the shower vertex position.
//
// Potential vertices in each view are determined and matched between the 2 or
// 3 available views. These matched vertices are then added to the art::Event
// as a vector of recob::Vertex objects.
//
/////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"

#include "larreco/RecoAlg/ImagePatternAlgs/PointIdAlg/PointIdAlg.h"
#include "larreco/RecoAlg/ImagePatternAlgs/ShowerVertex/VertexCandidate.h"
#include "larreco/RecoAlg/ImagePatternAlgs/ShowerVertex/HitCNNOutput.h"
#include "lardata/ArtDataHelper/MVAReader.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "TH2D.h"

#include <memory>

namespace nnet {

  class ShowerVertexFinder : public art::EDProducer {
    public:

      // these types to be replaced with use of feature proposed in redmine #12602
      typedef std::unordered_map< unsigned int, std::vector< size_t > > view_keymap;
      typedef std::unordered_map< unsigned int, view_keymap > tpc_view_keymap;
      typedef std::unordered_map< unsigned int, tpc_view_keymap > cryo_tpc_view_keymap;

      struct PiZero {

        PiZero(){};
        PiZero(unsigned int piZeroIndex){fPiZeroIndex = piZeroIndex;};
        PiZero(const PiZero& rhs){
          fPiZeroIndex = rhs.fPiZeroIndex;
          fPhotonIndex = rhs.fPhotonIndex;
          fMatchedPhoton = rhs.fMatchedPhoton;
          fShouldReco = rhs.fShouldReco;
        };
        
        unsigned int fPiZeroIndex;
        std::vector<unsigned int> fPhotonIndex;
        std::vector<bool> fMatchedPhoton;
        std::vector<bool> fShouldReco;
      };

      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Table<nnet::PointIdAlg::Config> PointIdAlg {
          Name("PointIdAlg")
        };

        fhicl::Atom<art::InputTag> WireLabel {
          Name("WireLabel"),
            Comment("tag of deconvoluted ADC on wires (recob::Wire)")
        };

        fhicl::Atom<art::InputTag> HitModuleLabel {
          Name("HitModuleLabel"),
            Comment("tag of hits to be EM/track tagged")
        };

        fhicl::Atom<art::InputTag> InputMVALabel {
          Name("InputMVALabel"),
            Comment("Name of the MVA used to produce the tagged clusters")
        };

        fhicl::Atom<float> InputMVACut {
          Name("InputMVACut"),
            Comment("Cut value below which to accept a cluster as shower like")
        };
        
        fhicl::Atom<float> VertexMVACut {
          Name("VertexMVACut"),
            Comment("Cut value above which to accept hit as a vertex")
        };

        fhicl::Sequence<int> Views {
          Name("Views"),
            Comment("tag clusters in selected views only, or in all views if empty list")
        };

        fhicl::Atom<float> DriftLimit {
          Name("DriftLimit"),
          Comment("Tolerence in drift coordinate to associate vertices in different views")
        };       

        fhicl::Atom<float> Match3DLimit {
          Name("Match3DLimit"),
          Comment("Distance tolerence between two predicted 3D vertices using different pairs of view positions to allow a match in 3D")
        };      

        fhicl::Atom<bool> IgnoreCosmics {
          Name("IgnoreCosmics"),
          Comment("Ignore tracks tagged as cosmics")
        };

        fhicl::Atom<bool> IgnoreLongTracks {
          Name("IgnoreLongTracks"),
          Comment("Ignore long tracks even if not tagged as cosmic")
        };

        fhicl::Atom<float> LongTrackCut {
          Name("LongTrackCut"),
          Comment("Maximum track length before it is ignored")
        };

        fhicl::Atom<art::InputTag> CosmicModuleLabel {
          Name("CosmicModuleLabel"),
          Comment("Name of the module that produced anab::Cosmic tags")
        };

        fhicl::Atom<bool> CalcPiZeroEfficiency {
          Name("CalcPiZeroEfficiency"),
          Comment("Calculate the efficiency of reconstructing pi-zero decay photon conversion vertices")
        };

        fhicl::Atom<float> TruthMatchDist {
          Name("TruthMatchDist"),
          Comment("Distance between true and reconstructed photon conversion vertex allowed")
        };

        fhicl::Atom<bool> SaveHitMaps {
          Name("SaveHitMaps"),
          Comment("Save the (wire,drift) hit maps for each event")
        };

      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit ShowerVertexFinder(Parameters const & p);

      ShowerVertexFinder(ShowerVertexFinder const &) = delete;
      ShowerVertexFinder(ShowerVertexFinder &&) = delete;
      ShowerVertexFinder & operator = (ShowerVertexFinder const &) = delete;
      ShowerVertexFinder & operator = (ShowerVertexFinder &&) = delete;

      void produce(art::Event & e) override;
      void beginJob() override;
      void endJob() override;

      void GetVertexCandidates(const art::Event& evt, const std::map<unsigned int,std::vector<std::pair<geo::WireID,float> > > &input, std::vector<VertexCandidate> &output); 
      void GetFinalCandidates(std::vector<VertexCandidate> &input, std::map<float,VertexCandidate> &output);

    private:
      template<size_t N> void ReadMVA(std::vector<float> &weights, art::Event& evt);

      bool Check3DPosition(const art::Event &evt, TVector3 pos, float drift, unsigned int missingView);

      void GetEfficiency(art::ValidHandle< std::vector<simb::MCParticle> > p, std::vector<VertexCandidate> const &v);

      bool IsFiducial(float x, float y, float z, float dx, float dy, float dz);
      void GetDetectorLimits();

      bool isViewSelected(int view) const;

      PointIdAlg fPointIdAlg;

      art::InputTag fWireProducerLabel;
      art::InputTag fHitModuleLabel;
      std::vector< int > fViews;

      art::InputTag fInputMVALabel; // Track / Shower MVA name
      float fInputMVACut;  // Cut to select shower-like clusters
      float fVertexMVACut;

      float fDriftLimit;
      float fMatch3DLimit;

      bool fIgnoreCosmics;
      bool fIgnoreLongTracks;
      float fLongTrackCut;
      art::InputTag fCosmicModuleLabel;
      bool fCalcPiZeroEff;
      float fTruthMatchDist;
      bool fSaveHitMaps;

      std::string fMVAClusterLabel; // Get this from the MVA itself.
      art::Handle<std::vector<recob::Cluster> > fMVAClusters;

      // Values for debugging performance
      unsigned int fNTrueConvs; // Count all photon conversions with at least 40 MeV KE
      unsigned int fNRecoConvs; // Count those true conversions matched to reco

      std::vector<float> fDimensionsMin;
      std::vector<float> fDimensionsMax;

      TH1D* fVtxDist;
      TH1D* fDirXAll;
      TH1D* fDirXVtx;
      TH1D* fDirYAll;
      TH1D* fDirYVtx;
      TH1D* fPhotDistAll;
      TH1D* fPhotDistVtx;
      TH1D* fPhotAngAll;
      TH1D* fPhotAngVtx;
      TH1D* fEnAll;
      TH1D* fEnVtx;
 
      ShowerVertexFinder::cryo_tpc_view_keymap fHitMap;
  };
  // ------------------------------------------------------

  ShowerVertexFinder::ShowerVertexFinder(ShowerVertexFinder::Parameters const& config) :
    fPointIdAlg(config().PointIdAlg()),
    fWireProducerLabel(config().WireLabel()),
    fHitModuleLabel(config().HitModuleLabel()),
    fViews(config().Views()),
    fInputMVALabel(config().InputMVALabel()),
    fInputMVACut(config().InputMVACut()),
    fVertexMVACut(config().VertexMVACut()),
    fDriftLimit(config().DriftLimit()),
    fMatch3DLimit(config().Match3DLimit()),
    fIgnoreCosmics(config().IgnoreCosmics()),
    fIgnoreLongTracks(config().IgnoreLongTracks()),
    fLongTrackCut(config().LongTrackCut()),
    fCosmicModuleLabel(config().CosmicModuleLabel()),
    fCalcPiZeroEff(config().CalcPiZeroEfficiency()),
    fTruthMatchDist(config().TruthMatchDist()),
    fSaveHitMaps(config().SaveHitMaps())
    {
      produces< std::vector<recob::Vertex> >();
      fNTrueConvs = 0;
      fNRecoConvs = 0;
    }
  // ------------------------------------------------------

  // Access the information from the MVA method. Very similar to the code in PMAlgTrackMaker_module
  template<size_t N> void ShowerVertexFinder::ReadMVA(std::vector<float> &weights, art::Event& evt){
    // Access the MVA. It could have 4 or 3 outputs so try to see which.
    auto cluResults = anab::MVAReader<recob::Cluster,N>::create(evt,fInputMVALabel);
    if (!cluResults){
      // We don't seem to have an MVA.
      return;
    }

    // This is the best way to get hold of the MVA clusters and tag used to associate the clusters
    // back to the underlying hit objects.
    fMVAClusters = cluResults->dataHandle();
    fMVAClusterLabel = cluResults->dataTag();

    // Now we have the MVA, we need to access the outputs
    int trkLikeIdx = cluResults->getIndex("track");
    int emLikeIdx = cluResults->getIndex("em");
    if ((trkLikeIdx < 0) || (emLikeIdx < 0)) {
      // If the variables we want don't exist then we have a problem.
      return;
    }

    // Actually extract the weights
    const auto & cnnOuts = cluResults->outputs();
    for (size_t i = 0; i < cnnOuts.size(); ++i){

      double trkOrEm = cnnOuts[i][trkLikeIdx] + cnnOuts[i][emLikeIdx];
      double val = 0;
      if (trkOrEm > 0){
        // Make sure output is normalised to fall between 0 and 1. 
        val = cnnOuts[i][trkLikeIdx] / trkOrEm;
      }
      weights.push_back(val);
    }
  }

  void ShowerVertexFinder::produce(art::Event & evt)
  {
    // Make sure we don't have any left overs from the previous event
    fHitMap.clear();

    mf::LogVerbatim("ShowerVertexFinder") << "next event: " << evt.run() << " / " << evt.id().event();

    // We may wish to make some plots
    art::ServiceHandle<art::TFileService> tfs;

    // Read the results from the track / shower CNN
    // These are the MVA weights
    std::vector<float> mvaWeights;

    if(fInputMVACut != 0){
      // Access the MVA. It could have 4 or 3 outputs so try to see which.
      ReadMVA<4>(mvaWeights,evt);
      if(mvaWeights.size() == 0){
        // If we didn't get any weights then try the version with three outputs instead
        ReadMVA<3>(mvaWeights,evt);
      }
    }

    std::vector<unsigned int> presentViews;
    std::vector<unsigned int> hitsToIgnore;

    // Get the hits to ignore according to the MVA
    if(mvaWeights.size() != 0){

      const art::FindManyP<recob::Hit> hitsFromClusters(fMVAClusters,evt,fMVAClusterLabel);

      // At this stage we have the hits from the clusters and the tag for each cluster.
      // Loop over the clusters and run the shower vertex CNN for each cluster
      for (size_t c = 0; c != fMVAClusters->size(); ++c){

        // If this is a shower cluster then don't do anything
        if(mvaWeights[c] <= fInputMVACut){
          continue;
        }

        // Get the hits from this cluster
        auto const& hits = hitsFromClusters.at(c);

        // Add all hits from this cluster to a vector
        for (auto const & h : hits){
          unsigned int view = h->WireID().Plane;
          if (!isViewSelected(view)){
            continue;
          }
          hitsToIgnore.push_back(h.key());
        } // End loop over hits
      } // End for loop over clusters
    }
    unsigned int nMVAIgnore = hitsToIgnore.size();

    // Get the hits to ignore from cosmics
    if(fIgnoreCosmics || fIgnoreLongTracks){
      // Get the reconstructed tracks
      auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(fCosmicModuleLabel);
      // We need the association between the tracks and the hits
      const art::FindManyP<recob::Hit> findTrackHits(recoTracks, evt, fCosmicModuleLabel);
      // Also want the association between tracks and anab::CosmicTag objects
      const art::FindManyP<anab::CosmicTag> findCosmicTags(recoTracks,evt,fCosmicModuleLabel);

      // Loop over the tracks
      for(unsigned int t = 0; t < recoTracks->size(); ++t){
        // Look for a cosmic tag object
        auto const& trackCosmicTag = findCosmicTags.at(t);

        float thisLength = ((*recoTracks)[t].Vertex() - (*recoTracks)[t].End()).Mag();

        bool isCosmic = (trackCosmicTag.size() != 0);
        bool isLong = thisLength > 100.*fLongTrackCut;

        if(!fIgnoreCosmics) isCosmic = false;
        if(!fIgnoreLongTracks) isLong = false;

        if(!isCosmic && !isLong) continue;

        if(fIgnoreLongTracks) mf::LogDebug("ShowerVertexFinder") << "Ignoring hits from track " << t << ", with length " << thisLength << " longer than " << 100.*fLongTrackCut << std::endl;

        // Get the cosmic tagged track hits
        auto const& trackHits = findTrackHits.at(t);

        // Add these hit to the vector to ignore
        for (auto const & h : trackHits){
          unsigned int view = h->WireID().Plane;
          if (!isViewSelected(view)){
            continue;
          }
          hitsToIgnore.push_back(h.key());
        }
      }
    }
    unsigned int nCosmicIgnore = hitsToIgnore.size() - nMVAIgnore;

    mf::LogDebug("ShowerVertexFinder") << "Hits ignored: " << std::endl;
    mf::LogDebug("ShowerVertexFinder") << " - From MVA     = " << nMVAIgnore << std::endl;
    mf::LogDebug("ShowerVertexFinder") << " - From Cosmics = " << nCosmicIgnore << std::endl;

    // Next job is to fill a hit map to form the basis of the images
    unsigned int cryo, tpc, view;

    auto hitListHandle = evt.getValidHandle< std::vector<recob::Hit> >(fHitModuleLabel);
    std::vector< art::Ptr<recob::Hit> > hitPtrList;
    art::fill_ptr_vector(hitPtrList, hitListHandle);
    size_t nHits = 0;
    for (auto const& h : hitPtrList)
    {
      view = h->WireID().Plane;
      if (!isViewSelected(view)) continue;
          
      // Store the view if we haven't already noted it.
      if(std::find(presentViews.begin(),presentViews.end(),view)==presentViews.end()){
        presentViews.push_back(view);
      }

      cryo = h->WireID().Cryostat;
      tpc = h->WireID().TPC;

      // Should we ignore this hit due to information from MVA or cosmic tag?
      if(std::find(hitsToIgnore.begin(),hitsToIgnore.end(),h.key())!=hitsToIgnore.end()){
        continue;
      }

      fHitMap[cryo][tpc][view].push_back(h.key());
      ++nHits;
    }

    std::map<unsigned int, TH2D*> histMap;
    if(fSaveHitMaps){
      // Make hit maps for all views.
      // Values filled as:
      // 0 = no hit (could be no deposit or because it was ignored)
      // 1 = hit but no vertex
      // 2 = hit with 1-view vertex
      // 3 = hit with 2-view vertex
      // 4 = hit with 3-view vertex
      for(unsigned int p = 0; p < presentViews.size(); ++p){
        std::stringstream histName;
        histName << "hitMap_" << evt.event() << "_" << presentViews[p];
        TH2D* temp = tfs->make<TH2D>(histName.str().c_str(),";Wire;Drift (cm)",1200,0,1200,720,-360,360);
        histMap.insert(std::make_pair(presentViews[p],temp));
      }
    }

    mf::LogDebug("ShowerVertexFinder") << "Hits remaining = " << nHits << std::endl;

    // Now we need to classify the shower hits with the shower vertex CNN
    auto wireHandle = evt.getValidHandle< std::vector<recob::Wire> >(fWireProducerLabel);

    if(!wireHandle.isValid()){
      mf::LogWarning("ShowerVertexFinder") << "Wire handle not valid... skipping event." << std::endl;
      return;
    }

    // We want to make a vector of recob::Vertex objects

    auto const* detProp = lar::providerFrom<detinfo::DetectorPropertiesService>();

    std::vector<VertexCandidate> finalMatchedVertices;

    // Store our vertices here
    auto allVertices = std::make_unique< std::vector< recob::Vertex > >(); 
    
    std::vector< char > hitInFA(hitPtrList.size(), 0); // tag hits in fid. area as 1, use 0 for hits close to the projectrion edges
    for (auto const & pcryo : fHitMap)
    {
      cryo = pcryo.first;
      for (auto const & ptpc : pcryo.second)
      {
        tpc = ptpc.first;

        // Store the potential verticies we find in each view
        // Do this more generically soon... some kind of map?
        std::map<unsigned int, std::vector<std::pair<geo::WireID,float> > > vertices;

        for (auto const & pview : ptpc.second)
        {
          view = pview.first;
          if (!isViewSelected(view)) continue; // should not happen, hits were selected

          fPointIdAlg.setWireDriftData(*wireHandle, view, tpc, cryo);

          std::vector<HitCNNOutput> allHitCNNVals;
          auto const* geom = lar::providerFrom<geo::Geometry>();

          // How many calls do we make to the CNN out of interest?
          unsigned int nCalls = 0;

          // (1) do all hits in this plane ------------------------------------------------
          for (size_t h : pview.second) // h is the Ptr< recob::Hit >::key()
          {
            const recob::Hit & hit = *(hitPtrList[h]);
            unsigned int thisWire = hit.WireID().Wire;
            float thisTime = hit.PeakTime();

            // Create a HitCNNOutput object to store all local values for this hit
            // This is because we can't expect the conversion point to be in one
            // exact pixel. For this reason we allow the network to check neighbouring
            // points, giving a 3x3 region centred on the hit position
            HitCNNOutput hitCNNVals(hit.WireID(),thisTime,hit.RMS());
            allHitCNNVals.push_back(hitCNNVals);
            // The current object will always be the one on the end of the vector
            unsigned int endPos = allHitCNNVals.size() - 1;
            //std::cout << " - Hit " << cryo << ", " << tpc << ", " << view << ", " << endPos << " :: " << thisWire << ", " << thisTime << std::endl;
            // Loop over +/- 1 wires and time +/- timeRMS
            for(unsigned int w = thisWire - 1; w <= thisWire + 1; ++w){
              if(w < 0 || w >= geom->Nwires(view,tpc,cryo)){
                // This wire doesn't exist, so we should move on. The three values will be filled
                // with a default value of -999.
                continue;
              }

              for(short t = -1; t <= 1; ++t){
                float time = thisTime + t*hit.RMS();
                float thisCNNVal = -999.;
                // Before checking the CNN, check this object and the other objects to see if we have
                // the same CNN value stored. Iterate backwards as most likely to find a match within
                // this object
                for(unsigned int hc = endPos; hc > 0; --hc ){
                  // Check all 9 bins
                  for(unsigned int b = 0; b < 9; ++b){
                    // Must have the same wire
                    if(w != allHitCNNVals[hc].GetWireFromIndex(b)) continue;
                    // If this output value isn't set then move on
                    if(allHitCNNVals[hc].GetCNNValue(b) < -998.) continue;
                    // If these two points are from the same patch then they
                    // must have the same output
                    if(fPointIdAlg.isSamePatch(w,time,w,allHitCNNVals[hc].GetTimeFromIndex(b))){
                      if(w == thisWire && t == 0){
                        thisCNNVal = allHitCNNVals[hc].GetCNNValue(b);
                        allHitCNNVals[hc].SetCNNValue(b,0.0);
                      }
                      else{
                        // Set this value to zero since it has already been accounted for
                        thisCNNVal = 0.0;
                      }
                      break; 
                    }
                  }
                } // End loop over other HitCNNOutput objects

                // Finally, access the CNN if we need to
                if(thisCNNVal < -998.){
                  thisCNNVal = 999; // Use this as a flag to identify where we need calls later
                }

                allHitCNNVals[endPos].AddOutput(w,time,thisCNNVal);
              }

            }
          } // Finished creating the HitCNNOutput objects

          // It is much more efficient to query the CNN with a vector of points as opposed to one by one.
          // Loop over the HitCNNOutput objects and make a vector of <idx,bin> calls to pass to the network.
          std::vector<std::pair<unsigned int,unsigned short>> prepForCNN;
          std::vector<std::pair<unsigned int,float>> valuesForCNN;

          for(unsigned int v = 0; v < allHitCNNVals.size(); ++v){
            for(unsigned short b = 0; b < 9; ++b){
              if(allHitCNNVals[v].GetCNNValue(b) > 998){
                prepForCNN.push_back(std::make_pair(v,b));
                valuesForCNN.push_back(std::make_pair(allHitCNNVals[v].GetWireFromIndex(b),allHitCNNVals[v].GetTimeFromIndex(b)));
              }
            }
          }

          // Now lets call the CNN once with this vector
          std::vector<std::vector<float>> cnnOutputs = fPointIdAlg.predictIdVectors(valuesForCNN);   
          nCalls = cnnOutputs.size();
          // Associate these points back to the HitCNNOutput objects
          for(unsigned int v = 0; v < cnnOutputs.size(); ++v){
            allHitCNNVals[prepForCNN[v].first].SetCNNValue(prepForCNN[v].second,cnnOutputs[v][0]);
          }
       
          // Loop over the HitCNNOutput objects
          for(auto const hc : allHitCNNVals){

            // Want to convert the peak time into drift position.
            double x = detProp->ConvertTicksToX(hc.GetTime(),hc.GetWireID().planeID());

            if(hc.GetMaxCNNValue() > fVertexMVACut){
//              std::cout << " We dun found a vertex with scores " << cnnOutputMin[0] << ", " << cnnOutput[0] << ", " << cnnOutputMax[0] << std::endl;
              std::pair<geo::WireID,float> vtx = std::make_pair(hc.GetWireID(),x);

              if(vertices.find(view) == vertices.end()){
                std::vector<std::pair<geo::WireID,float> > temp;
                temp.push_back(vtx);
                vertices.insert(std::make_pair(view,temp));
              }
              else{
                vertices[view].push_back(vtx);
              }
              if(fSaveHitMaps) histMap[view]->SetBinContent(histMap[view]->FindBin(hc.GetWire(),x),2);
            }
            else{
              if(fSaveHitMaps) histMap[view]->SetBinContent(histMap[view]->FindBin(hc.GetWire(),x),1);
            }

          } // End loop over HitCNNObjects
          mf::LogDebug("ShowerVertexFinder") << "Finished finding vertices for view " << view << " after " << nCalls << " calls to the CNN for " << allHitCNNVals.size() << " hits " << std::endl;
        } // End loop over views

        // We can only match vertices if we find them in more than one view
        unsigned int n0 = 0;
        unsigned int n1 = 0;
        unsigned int n2 = 0;
        if(vertices.find(0) != vertices.end()) n0 = vertices[0].size();
        if(vertices.find(1) != vertices.end()) n1 = vertices[1].size();
        if(vertices.find(2) != vertices.end()) n2 = vertices[2].size();
        mf::LogDebug("ShowerVertexFinder") << "Vertices in TPC " << tpc << ": " << n0 << ", " << n1 << ", " << n2 << std::endl;

        if(vertices.size() > 1){
      
          // Now we want to try to match the vertices between views
          std::vector<VertexCandidate> matchedVertices;
          GetVertexCandidates(evt,vertices,matchedVertices);
    
          // Using the score as the key is helpful as they are then sorted by default
          // and a low score is best!
          std::map<float,VertexCandidate> scoreMap;
          GetFinalCandidates(matchedVertices,scoreMap);
    
          // Loop through the map once more to print out the final matches and store in the vector
          if(scoreMap.size()){
            mf::LogDebug("ShowerVertexFinder") << "- Matches for TPC " << tpc << std::endl;
          }
          for (auto const &f : scoreMap){
            if(f.second.IsFailed()) continue;
            finalMatchedVertices.push_back(f.second);
            f.second.Print();
            if(fSaveHitMaps){
              for(unsigned int i = 0; i < f.second.GetSize(); ++i){
                histMap[f.second.GetView(i)]->SetBinContent(histMap[f.second.GetView(i)]->FindBin(f.second.GetWire(i).Wire,f.second.GetDrift(i)),f.second.GetSize()+1);
              }
            }
          }
        }

      } // End loop over tpcs
    } // End loop over cryostats

    // Get the matching efficiency and fill the reco - true plots if we have MC and request it
    if(fCalcPiZeroEff && !evt.isRealData()){
      auto particleHandle = evt.getValidHandle< std::vector<simb::MCParticle> >("largeant");
      GetEfficiency(particleHandle,finalMatchedVertices);
    }

    // Make the recob::Vertex objects
    for(auto const v : finalMatchedVertices){
      double xyz[3];
      v.Get3DPos().GetXYZ(xyz);
      size_t vidx = allVertices->size();
      allVertices->push_back(recob::Vertex(xyz, vidx));
    }
 
    // Write the outputs to the event
    evt.put(std::move(allVertices));

  }

  void ShowerVertexFinder::GetVertexCandidates(const art::Event& evt, const std::map<unsigned int,std::vector<std::pair<geo::WireID,float> > > &inputMap, std::vector<VertexCandidate> &vtxCands){
    
    mf::LogDebug("ShowerVertexFinder") << " - Preparing candidates " << std::endl;

    if(vtxCands.size()){
      vtxCands.clear();
    }

    // Let's save some characters.
    typedef std::map<unsigned int,std::vector<std::pair<geo::WireID,float> > > thisMap; 
    typedef thisMap::const_iterator thisIt; 

    thisIt m1 = inputMap.begin();
    thisIt m2 = m1; ++m2;
    thisIt m3 = m2; ++m3;
    
    std::vector<VertexCandidate> finalCands;
    // Look at the vector for the first view
    for(auto const &v1 : (*m1).second){

      // Use this candidate as something to compare to.
      VertexCandidate thisCand(v1.first,v1.second,(*m1).first);
//      std::cout << "View 0 Cand " << v1.first << ", " << v1.second << std::endl;

      // Look at the second view vector
      // Look at the vector
      for(auto const &v2 : (*m2).second){
        if(thisCand.IsCompatible(v2.first,v2.second,fDriftLimit)){
          // Need a copy of our candidate so that we don't change the original
          // since we might need it later
          VertexCandidate candCopy1(thisCand);
          candCopy1.AddView(v2.first,v2.second,(*m2).first);
          // Can we find a match in view 3 too? 
          //bool match3D = false;
          if(inputMap.size() == 3){
            for(auto const &v3 : (*m3).second){
              // Make another copy incase we have multiple matches
              VertexCandidate candCopy2(candCopy1);
              if(candCopy2.IsCompatible(v3.first,v3.second,fDriftLimit)){
                candCopy2.AddView(v3.first,v3.second,(*m3).first);
                finalCands.push_back(candCopy2);
                //match3D = true;
                //any3DMatch = true;
//                std::cout << "  - 3D match made " << std::endl;
              }
            }
          }
          // Save the 2D vertex too
          finalCands.push_back(candCopy1);
//          std::cout << "  - 2D match made " << std::endl;
        }
      } // End loop over the second view vector

      // How about all those matches between just views 1 and 3?
      if(inputMap.size() == 3){
        // Make a copy
        for(auto const &v3 : (*m3).second){
          VertexCandidate candCopy1(thisCand);
//          std::cout << "View 2 Cand " << v3.first << ", " << v3.second << std::endl;
          if(candCopy1.IsCompatible(v3.first,v3.second,fDriftLimit)){
            candCopy1.AddView(v3.first,v3.second,(*m3).first);
            finalCands.push_back(candCopy1);
//            std::cout << "View 2 matched" << std::endl;
          }
        } // End loop over the third view vector
      }

    } // End loop over first view vector

    // Also need to try matching between just planes two and three. 
    // This could make some repeats, but we can deal with them later.
    if(inputMap.size() == 3){
      for(auto const &v1 : (*m2).second){

        // Use this candidate as something to compare to.
        VertexCandidate thisCand(v1.first,v1.second,(*m2).first);
//        std::cout << "View 1 Cand " << v1.first << ", " << v1.second << std::endl;

        // Look at the second view vector
        for(auto const &v2 : (*m3).second){
          // Make a copy
          VertexCandidate candCopy1(thisCand);
          if(candCopy1.IsCompatible(v2.first,v2.second,fDriftLimit)){
            candCopy1.AddView(v2.first,v2.second,(*m3).first);
            finalCands.push_back(candCopy1);
//            std::cout << "View 2 Cand " << v2.first << ", " << v2.second << std::endl;
          }
        } // End loop over the second view vector
    
      } // End loop to find matches between views 2 and 3
    }

    // We can now push back all of our vertices into the main list
    for(auto const &thisVtx : finalCands){
      if(thisVtx.GetSize() == 3){
        vtxCands.push_back(thisVtx);
      }
      else{
        // Make sure that any 2D vertices have a corresponding hit in the third view
        unsigned int missingView = 999;
        if(thisVtx.HasView(0) && thisVtx.HasView(1)) missingView = 2;
        if(thisVtx.HasView(0) && thisVtx.HasView(2)) missingView = 1;
        if(thisVtx.HasView(1) && thisVtx.HasView(2)) missingView = 0;
        // Returns true if good
        if(Check3DPosition(evt,thisVtx.Get3DPos(),thisVtx.GetDrift(0),missingView)){
          vtxCands.push_back(thisVtx);
        }
      } 
    }

  }

  bool ShowerVertexFinder::Check3DPosition(const art::Event& evt, TVector3 pos, float drift, unsigned int missingView){

    // For now just return true, the code below needs some work.
    return true;    

    // Get a geometry instance to convert a 3D position to a wire coordinate
    auto const* geom = lar::providerFrom<geo::Geometry>();
    // Get the TPC first
    geo::TPCID thisTPC = geom->FindTPCAtPosition(pos);
    // Now get the WireCoordinate
    double wire = geom->WireCoordinate(pos.Y(),pos.Z(),missingView,thisTPC.TPC,thisTPC.Cryostat);

    // Now we have our wire and time, look for a match in the hit vector
    auto hitListHandle = evt.getValidHandle< std::vector<recob::Hit> >(fHitModuleLabel);
    std::vector< art::Ptr<recob::Hit> > hitPtrList;
    art::fill_ptr_vector(hitPtrList, hitListHandle);
    for (auto const& h : hitPtrList){
      std::vector<size_t> keysOfInterest = fHitMap[thisTPC.Cryostat][thisTPC.TPC][missingView];
      if(std::find(keysOfInterest.begin(),keysOfInterest.end(),h.key()) != keysOfInterest.end()){
        // This hit is interesting, so we should see if it matches our position
        const recob::Hit & hit = *(hitPtrList[h.key()]);
        unsigned int thisWire = hit.WireID().Wire;
        float thisTime = hit.PeakTime();
        if(fabs(thisWire - wire) < 2 && fabs(thisTime - drift) < 2){
          return true;
        }
      }
    }

    mf::LogDebug("ShowerVertexFinder") << "2D vertex did not correspond to a hit in 3D, so was rejected. " << std::endl;
    return false;
  }

  // Function takes the candidates and weeds out obvious ones before adding those real candidates to a map.
  // It then flags a number of these as failed, leaving just the best matches.
  void ShowerVertexFinder::GetFinalCandidates(std::vector<VertexCandidate> &input, std::map<float,VertexCandidate> &output){

    mf::LogDebug("ShowerVertexFinder") << " - Finalising candidates " << std::endl;

    if(output.size()){
      output.clear();
    }

    for(auto &mv : input){
      // If the score is >9998 then it is a failed attempt so don't store.
      float score = mv.GetScore();
      if(score > 9998.){
        mv.SetFailed();
        continue;
      }

      // 3D vertices should predict the same 3D position using two pairs.
      // This should agree very nicely if it is a good match.
      if(mv.GetSize()==3){
        if(score > fMatch3DLimit){
          mv.SetFailed();
          continue;
        }
      }

      // Punish 2D vertices compared to 3D
      if(mv.GetSize()==2){
        score += 100;
      }
      output.insert(std::make_pair(score,mv));
    }

    // Now we have a map sorted by score, with obvious junk removed. Lowest will be the best 3D vertex, and highest (probably) the worst 2D one.
    // We now need to make sure to accept the best scores, and then fail other matches using the same WireID / time combinations.
    for(auto &m1 : output){
      if(m1.second.IsFailed()){
        continue;
      }

      // Look to see if any vertices in this candidate match another candidate
      for(unsigned int i = 0; i < m1.second.GetSize(); ++i){
        // Loop over the other entries in the map
        for(auto &m2 : output){
          if(m2.second.IsFailed()){
            continue;
          }
          if(m1.first == m2.first) continue;
          if(m2.second.ContainsVertex(m1.second.GetWire(i),m1.second.GetDrift(i),m1.second.GetView(i))){
            // m2, which has a worse score than m1, contains a wire that is the same. Set it to false and move on.       
            m2.second.SetFailed();
            continue;
          }
        }
      }
    }

  }

  bool ShowerVertexFinder::IsFiducial(float x, float y, float z, float dx, float dy, float dz){
  
    if(fDimensionsMin.size() == 0){
      this->GetDetectorLimits();
    }

    // If the position is literally outside the detector then it is obviously a failure.
    if(x < fDimensionsMin[0] || x > fDimensionsMax[0]) return false;
    if(y < fDimensionsMin[1] || y > fDimensionsMax[1]) return false;
    if(z < fDimensionsMin[2] || z > fDimensionsMax[2]) return false;

    // If it is inside but close to the edge, will the direction save us?
    float fid = 50.;
    if(x < (fDimensionsMin[0] + fid) && dx < 0.0) return false;
    if(x > (fDimensionsMax[0] - fid) && dx > 0.0) return false;
    if(y < (fDimensionsMin[1] + fid) && dy < 0.0) return false;
    if(y > (fDimensionsMax[1] - fid) && dy > 0.0) return false;
    if(z < (fDimensionsMin[2] + fid) && dz < 0.0) return false;
    if(z > (fDimensionsMax[2] - fid) && dz > 0.0) return false;
    return true; 
  }

  void ShowerVertexFinder::GetDetectorLimits(){
  
    // Use the geometry to find the detector limits.
    // Need to find the minimum and maximum height values from the geometry.
    double minX = 1.e6;
    double maxX = -1.e6;
    double minY = 1.e6;
    double maxY = -1.e6;
    double minZ = 1.e6;
    double maxZ = -1.e6;
  
    auto const* geom = lar::providerFrom<geo::Geometry>();
  
    // Since we can stack TPCs, we can't just use geom::TPCGeom::Height()
    for (geo::TPCID const& tID: geom->IterateTPCIDs()) {
      geo::TPCGeo const& TPC = geom->TPC(tID);
  
      // get center in world coordinates
      double origin[3] = {0.};
      double center[3] = {0.};
      TPC.LocalToWorld(origin, center);
      double tpcDim[3] = {TPC.HalfWidth(), TPC.HalfHeight(), 0.5*TPC.Length() };
  
      // If the TPC is less than 1m wide, it is a dummy so we should ignore it.
      if (tpcDim[0] < 25.0){
        continue;
      }

      if( center[0] - tpcDim[0] < minX ) minX = center[0] - tpcDim[0];
      if( center[0] + tpcDim[0] > maxX ) maxX = center[0] + tpcDim[0];
      if( center[1] - tpcDim[1] < minY ) minY = center[1] - tpcDim[1];
      if( center[1] + tpcDim[1] > maxY ) maxY = center[1] + tpcDim[1];
      if( center[2] - tpcDim[2] < minZ ) minZ = center[2] - tpcDim[2];
      if( center[2] + tpcDim[2] > maxZ ) maxZ = center[2] + tpcDim[2];
    } // for all TPC
  
    fDimensionsMin.clear();
    fDimensionsMin.push_back(minX);
    fDimensionsMin.push_back(minY);
    fDimensionsMin.push_back(minZ);
    fDimensionsMax.clear();
    fDimensionsMax.push_back(maxX);
    fDimensionsMax.push_back(maxY);
    fDimensionsMax.push_back(maxZ);
  
  }

  void ShowerVertexFinder::GetEfficiency(art::ValidHandle< std::vector<simb::MCParticle> > particleHandle, std::vector<VertexCandidate> const &finalMatchedVertices){
    std::vector<TVector3> trueDecayPhotons;
    std::vector<TVector3> trueDecayDirs;
    std::vector<bool> shouldReco;
    std::vector<float> energy;

    // First things first, loop over the MCParticles to find all of the pizeros
    std::map<unsigned int,ShowerVertexFinder::PiZero> truePiZeros;
    for(unsigned int p = 0; p < particleHandle->size(); ++p){
      if((*particleHandle)[p].PdgCode() == 111){
        ShowerVertexFinder::PiZero newPiZero(p);
        truePiZeros.insert(std::make_pair(p,newPiZero));
      }
    }

    // Now loop again looking for photons and match them to a pizero
    for(unsigned int p = 0; p < particleHandle->size(); ++p){
      simb::MCParticle thisP = (*particleHandle)[p];
      if(thisP.PdgCode() == 22){
        // Check if the mother is a PiZero
        for(unsigned int m = 0; m < particleHandle->size(); ++m){
          simb::MCParticle thisM = (*particleHandle)[m];
          // Find the mother and see if it is a PiZero
          if(thisP.Mother() == thisM.TrackId() && thisM.PdgCode() == 111){
            truePiZeros[m].fPhotonIndex.push_back(p);
            truePiZeros[m].fMatchedPhoton.push_back(false);
            truePiZeros[m].fShouldReco.push_back(false);
          }
        }
      }
    }

    // We now have a complete map of piZeros and their decay photons.
    // Now iterate through the map and see which photons we expect to reconstruct
    for(auto &pizero : truePiZeros){
      simb::MCParticle thisPiZero = (*particleHandle)[pizero.second.fPiZeroIndex];
      mf::LogDebug("ShowerVertexFinder") << "Found true pi-zero decay at point (" << thisPiZero.EndY() << ", " << thisPiZero.EndZ() << ")" << std::endl;
      for(unsigned int p = 0; p < pizero.second.fPhotonIndex.size(); ++p){
        simb::MCParticle photon = (*particleHandle)[pizero.second.fPhotonIndex[p]];
        float energy = 1000.*(photon.E()-photon.Mass());
        mf::LogDebug("ShowerVertexFinder") << "Found true pi-zero decay photon at end point (" << photon.EndY() << ", " << photon.EndZ() << ") with energy " << energy << std::endl;
        // We should double check that these photons convert inside the detector
        bool isFid = IsFiducial(photon.EndX(),photon.EndY(),photon.EndZ(),photon.EndPx(),photon.EndPy(),photon.EndPz());
//        if(energy > 40.0 && isFid){
        if(isFid){
          pizero.second.fShouldReco[p] = true;
        
          // Fill some photon level histograms
          fDirXAll->Fill(photon.Momentum().Vect().Unit().X());
          fDirYAll->Fill(photon.Momentum().Vect().Unit().Y());
          fEnAll->Fill(energy);
        }
      } 
      // Fill some pion level histograms
      if(pizero.second.fShouldReco[0] && pizero.second.fShouldReco[1]){
        simb::MCParticle photon1 = (*particleHandle)[pizero.second.fPhotonIndex[0]];
        simb::MCParticle photon2 = (*particleHandle)[pizero.second.fPhotonIndex[1]];
        fPhotAngAll->Fill(photon1.EndPosition().Vect().Angle(photon2.EndPosition().Vect()));
        fPhotDistAll->Fill((photon1.EndPosition().Vect()-photon2.EndPosition().Vect()).Mag());
      }
    }

    // Now match to the reconstructed vertices
    std::vector<unsigned int> recoVtxUsed;
    for(auto &pizero : truePiZeros){
      for(unsigned int p = 0; p < pizero.second.fPhotonIndex.size(); ++p){
        // Get the photon
        simb::MCParticle photon = (*particleHandle)[pizero.second.fPhotonIndex[p]];

        // The photon conversion point
        TVector2 photonEnd(photon.EndY(),photon.EndZ());

        // Check if we should reconstruction this photon
        if(!pizero.second.fShouldReco[p]){
          mf::LogDebug("ShowerVertexFinder") << " Not expected to reconstruct gamma vertex (" << photonEnd.X() << ", " << photonEnd.Y() << ") " << std::endl;
          continue;
        }
        ++fNTrueConvs;

        double recoTrueScore = 9999;
        TVector2 bestMatch;
        unsigned int bestIndex = 9999;
        bool got3ViewVtx = false;
        // Loop over the 3-view vertices first
        for(unsigned int v = 0; v < finalMatchedVertices.size(); ++v){
          // Try to match 3D vertices first
          if(finalMatchedVertices[v].GetSize() != 3) continue;

          // If we used this already, then don't consider it again
          if(!(std::find(recoVtxUsed.begin(),recoVtxUsed.end(),v)==recoVtxUsed.end())) continue;

          TVector3 recoVtx = finalMatchedVertices[v].Get3DPos();
          TVector2 vtxPos(recoVtx.Y(),recoVtx.Z());
          float dist = (photonEnd-vtxPos).Mod();
          if(dist < recoTrueScore){
            recoTrueScore = dist;
            bestMatch = vtxPos;
            bestIndex = v;
          }
        }

        if(recoTrueScore < fTruthMatchDist){
          got3ViewVtx = true;
          recoVtxUsed.push_back(bestIndex);
        }

        // If no 3 view vertex, look for a 2 view one
        bool got2ViewVtx = false;
        if(!got3ViewVtx){
          for(unsigned int v = 0; v < finalMatchedVertices.size(); ++v){
            // Try to match 3D vertices first
            if(finalMatchedVertices[v].GetSize() != 2) continue;
  
            // If we used this already, then don't consider it again
            if(!(std::find(recoVtxUsed.begin(),recoVtxUsed.end(),v)==recoVtxUsed.end())) continue;
  
            TVector3 recoVtx = finalMatchedVertices[v].Get3DPos();
            TVector2 vtxPos(recoVtx.Y(),recoVtx.Z());
            float dist = (photonEnd-vtxPos).Mod();
            if(dist < recoTrueScore){
              recoTrueScore = dist;
              bestMatch = vtxPos;
              bestIndex = v;
            }
          }

          if(recoTrueScore < fTruthMatchDist){
            got2ViewVtx = true;
            recoVtxUsed.push_back(bestIndex);
          }

        } // Check for a 2 view vertex

        // If we matched a vertex, fill some plots
        if(got3ViewVtx || got2ViewVtx){
          ++fNRecoConvs;
          pizero.second.fMatchedPhoton[p] = true;
          fVtxDist->Fill(recoTrueScore);
          fDirXVtx->Fill(photon.Momentum().Vect().Unit().X());
          fDirYVtx->Fill(photon.Momentum().Vect().Unit().Y());
          fEnVtx->Fill(1000.*(photon.E()-photon.Mass()));
          mf::LogDebug("ShowerVertexFinder") << "Matched a true vertex (" << photonEnd.X() << ", " << photonEnd.Y() << ") to a reco vertex (" << bestMatch.X() << ", " << bestMatch.Y() << ")" << std::endl;
        }
        else{
          mf::LogDebug("ShowerVertexFinder") << "Didn't find a match for true vertex (" << photonEnd.X() << ", " << photonEnd.Y() << ")" << std::endl;

          //Get the true conversion point of this photon and look at where it should be in wire, time coordinates
          auto const* geom = lar::providerFrom<geo::Geometry>();

          geo::TPCID thisTPC = geom->FindTPCAtPosition(photon.EndPosition().Vect());
          if(thisTPC.TPC < 10000){
            double wireCoord[3] = {-999.,-999.,-999.};
            wireCoord[0] = geom->WireCoordinate(photonEnd.X(),photonEnd.Y(),0,thisTPC.TPC,thisTPC.Cryostat);
            wireCoord[1] = geom->WireCoordinate(photonEnd.X(),photonEnd.Y(),1,thisTPC.TPC,thisTPC.Cryostat);
            wireCoord[2] = geom->WireCoordinate(photonEnd.X(),photonEnd.Y(),2,thisTPC.TPC,thisTPC.Cryostat);
            mf::LogDebug("ShowerVertexFinder") << " - True position corresponds to wires " << wireCoord[0] << ", " << wireCoord[1] << ", " << wireCoord[2] << std::endl; 
          }
        }

      } // End loop over piZero photons

      // If we matched both photons
      if(pizero.second.fMatchedPhoton[0] && pizero.second.fMatchedPhoton[1]){
        simb::MCParticle photon1 = (*particleHandle)[pizero.second.fPhotonIndex[0]];
        simb::MCParticle photon2 = (*particleHandle)[pizero.second.fPhotonIndex[1]];
        fPhotDistVtx->Fill((photon1.EndPosition().Vect()-photon2.EndPosition().Vect()).Mag());
        fPhotAngVtx->Fill(photon1.EndPosition().Vect().Angle(photon2.EndPosition().Vect()));
      }
    } // End loop over piZeros

  }

  // ------------------------------------------------------

  bool ShowerVertexFinder::isViewSelected(int view) const
  {
    if (fViews.empty()) return true;
    else
    {
      bool selected = false;
      for (auto k : fViews) if (k == view) { selected = true; break; }
      return selected;
    }
  }
  // ------------------------------------------------------

  void ShowerVertexFinder::beginJob() {

    art::ServiceHandle<art::TFileService> tfs;
  
    if(fCalcPiZeroEff){  
      fVtxDist = tfs->make<TH1D>("rmtVtxDist",";Distance between reconstructed and true photon conversion point (cm)",50,0,1);
      fDirXAll = tfs->make<TH1D>("dirXAll",";Direction cosine x",50,-1,1); // For all true events
      fDirXVtx = tfs->make<TH1D>("dirXVtx",";Direction cosine x",50,-1,1); // For those matched to reco
      fDirYAll = tfs->make<TH1D>("dirYAll",";Direction cosine y",50,-1,1); // For all true events
      fDirYVtx = tfs->make<TH1D>("dirYVtx",";Direction cosine y",50,-1,1); // For those matched to reco
      fPhotDistAll = tfs->make<TH1D>("photDistAll","Distance to other photon (cm)",50,0,50);
      fPhotDistVtx = tfs->make<TH1D>("photDistVtx","Distance to other photon (cm)",50,0,50);
      fPhotAngAll = tfs->make<TH1D>("photAngAll","Angle between photons (degrees)",50,0,0.5*TMath::Pi());
      fPhotAngVtx = tfs->make<TH1D>("photAngVtx","Angle between photons (degrees)",50,0,0.5*TMath::Pi());
      fEnAll   = tfs->make<TH1D>("enAll",";Energy (MeV)",50,0,1000);
      fEnVtx   = tfs->make<TH1D>("enVtx",";Energy (MeV)",50,0,1000);
    }
  }

  void ShowerVertexFinder::endJob() {

    if(fNTrueConvs != 0){
      std::cout << "Pi-zero decay photon vertex efficiency = " << fNRecoConvs / (float)fNTrueConvs << "(" << fNRecoConvs << "/" << fNTrueConvs << ")" << std::endl;
    }
  }

  DEFINE_ART_MODULE(ShowerVertexFinder)

} // End nnet namespace

