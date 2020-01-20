// Chris Backhouse - bckhouse@fnal.gov

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"

#include "larcore/Geometry/Geometry.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h" // Uncompress()

#include <iostream>

#include "TGraph.h"
#include "TH1.h"

namespace hf
{
// ----------------------------------------------------------------------------
class HitPlotter: public art::EDAnalyzer
{
public:
  explicit HitPlotter(const fhicl::ParameterSet& pset);

  void analyze(const art::Event& evt) override;
  void endJob() override;

protected:
  art::ServiceHandle<geo::Geometry> fGeom;

  std::string fWireLabel;
  std::string fRawDigitLabel;
  std::vector<std::string> fHitLabels;
};

DEFINE_ART_MODULE(HitPlotter)

// ---------------------------------------------------------------------------
HitPlotter::HitPlotter(const fhicl::ParameterSet& pset)
: EDAnalyzer(pset),
  fWireLabel(pset.get<std::string>("WireLabel")),
  fRawDigitLabel(pset.get<std::string>("RawDigitLabel")),
  fHitLabels(pset.get<std::vector<std::string>>("HitLabels"))
{
}

// ----------------------------------------------------------------------------
void HitPlotter::endJob()
{
}

std::string safe_name(geo::WireID id)
{
  std::string ret(id);
  for(unsigned int i = 0; i < ret.size(); ++i) if(ret[i] == ' ' || ret[i] == ':') ret[i] = '_';
  return ret;
}

// ----------------------------------------------------------------------------
void HitPlotter::analyze(const art::Event& evt)
{
  art::ServiceHandle<art::TFileService> tfs;

  art::TFileDirectory evtdir = tfs->mkdir(Form("evt_%d", evt.event()));


  std::map<std::string, std::map<geo::WireID, std::vector<recob::Hit>>> hitmap;

  for(const std::string& label: fHitLabels){
    art::Handle<std::vector<recob::Hit>> hits;
    evt.getByLabel(label, hits);

    for(const recob::Hit& hit: *hits){
      const geo::WireID id = fGeom->ChannelToWire(hit.Channel())[0];
      hitmap[label][id].push_back(hit);
    }
  }


  std::map<geo::WireID, const raw::RawDigit*> digmap;

  art::Handle<std::vector<raw::RawDigit>> digs;
  evt.getByLabel(fRawDigitLabel, digs);

  if(!digs.failedToGet()){
    for(const raw::RawDigit& dig: *digs){
      digmap[fGeom->ChannelToWire(dig.Channel())[0]] = &dig;
    }
  }
  else if(!fRawDigitLabel.empty()){
    std::cout << "HitPlotter: Warning, RawDigits not found under label '"
              << fRawDigitLabel << "'" << std::endl;
  }


  art::Handle<std::vector<recob::Wire>> wires;
  evt.getByLabel(fWireLabel, wires);

  for(const recob::Wire& wire: *wires){
    const geo::WireID id = fGeom->ChannelToWire(wire.Channel())[0];
    art::TFileDirectory wiredir = evtdir.mkdir(safe_name(id));

    std::vector<short> adcs;
    const raw::RawDigit* dig = digmap[id];
    if(dig) raw::Uncompress(dig->ADCs(), adcs, dig->Compression());

    for(const lar::sparse_vector<float>::datarange_t& range: wire.SignalROI().get_ranges()){

      art::TFileDirectory rangedir = wiredir.mkdir(TString::Format("range_%lu_%lu", range.begin_index(), range.end_index()).Data());

      TH1* hwire = rangedir.make<TH1F>("wire", "", range.end_index()-range.begin_index()+1, range.begin_index()-.5, range.end_index()+.5);

      int tick = range.begin_index();
      for(float y: range.data()) hwire->Fill((tick++), y);

      TH1* hdig = rangedir.make<TH1F>("dig", "", range.end_index()-range.begin_index()+1, range.begin_index()-.5, range.end_index()+.5);

      for(unsigned int tick = range.begin_index();
          tick < std::min(range.end_index()+1, adcs.size());
          ++tick){
        const int adc = adcs[tick] ? int(adcs[tick])-dig->GetPedestal() : 0;
        hdig->Fill(tick, adc);
      }

      for(const std::string& label: fHitLabels){
        // Figure out whether to draw wires or raw digits based on what was fit
        art::Handle<art::Assns<recob::Hit, recob::Wire>> passn;
        evt.getByLabel(label, passn);
        const bool plotRawDigs = passn.failedToGet();

        art::TFileDirectory hitdir = rangedir.mkdir(label);

        TH1* hsum = hitdir.make<TH1F>("sum", "", range.end_index()-range.begin_index()+1, range.begin_index()-.5, range.end_index()+.5);

        int hitIdx = 0;
        for(const recob::Hit& hit: hitmap[label][id]){
          if(hit.PeakTime() > range.begin_index() && hit.PeakTime() < range.end_index()){
            TGraph* g = hitdir.make<TGraph>();
            g->SetName(TString::Format("hit_%d", hitIdx++).Data());
            for(double x = -3; x < +3.05; x += .1){

              double y;
              if(hit.SignalType() == geo::kCollection || !plotRawDigs)
                y = hit.PeakAmplitude() * exp(-x*x/2);
              else
                y = hit.PeakAmplitude() * -x*exp(-x*x/2);

              g->SetPoint(g->GetN(), hit.PeakTime()+x*hit.RMS(), y);
            }
            g->Write();
          }
        }

        for(int i = 1; i <= hsum->GetNbinsX(); ++i){
          const double xi = hsum->GetBinCenter(i);

          for(const recob::Hit& hit: hitmap[label][id]){
            const double x = (xi-hit.PeakTime())/hit.RMS();

            double y;
            if(hit.SignalType() == geo::kCollection || !plotRawDigs)
              y = hit.PeakAmplitude() * exp(-x*x/2);
            else
              y = hit.PeakAmplitude() * -x*exp(-x*x/2);

            hsum->Fill(xi, y);
          }
        }
      }
    } // end for range

  } // end for wire
}

} // namespace
