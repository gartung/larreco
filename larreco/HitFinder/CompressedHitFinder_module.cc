// Chris Backhouse - bckhouse@fnal.gov

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"

#include "larcore/Geometry/Geometry.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"


#include <iostream>

#include "TMatrixD.h"
#include "TVectorD.h"

const double L1_strength = .1; // TODO

const double kRMS = 2.75; // TODO

// ----------------------------------------------------------------------------
/// TMatrixD::operator() does various sanity checks and shows up in profiles
class TMatrixDFast: public TMatrixD
{
public:
  TMatrixDFast(int N, int M) : TMatrixD(N, M) {}

  inline double operator()(int i, int j) const
  {
    return fElements[fNrows*i+j];
  }
  inline double& operator()(int i, int j)
  {
    return fElements[fNrows*i+j];
  }
};

// ----------------------------------------------------------------------------
class TVectorDFast: public TVectorD
{
public:
  TVectorDFast(int N) : TVectorD(N) {}
  TVectorDFast(const TVectorD& v) : TVectorD(v) {}

  // Skip bounds-checking overhead
  inline double  operator[](unsigned int i) const {return fElements[i];}
  inline double& operator[](unsigned int i)       {return fElements[i];}
};

// ----------------------------------------------------------------------------
class GausKernel
{
public:
  GausKernel(double rms, int range) : fRange(range), fData(2*range+1)
  {
    const double norm = 1/(rms*sqrt(2*M_PI));
    for(int i = 0; i < 2*range+1; ++i){
      const double x = (i-range)/rms;
      fData[i] = norm*exp(-x*x/2);
    }
  }

  inline double y(int dx) const {return fData[fRange+dx];}

  inline int Range() const {return fRange;}
protected:
  int fRange;
  std::vector<double> fData;
};

// ----------------------------------------------------------------------------
inline double chisq(const TVectorDFast& y,
                    const TVectorDFast& pred,
                    const TVectorDFast& alpha)
{
  return (pred-y).Norm2Sqr() + L1_strength * alpha.Norm1();
}

// ----------------------------------------------------------------------------
// w/ pts at 0, 0.5, 1
double find_min(double y0, double y1, double y2, double* pred)
{
  const double c = y0;
  //  a*.25 + b*.5 + c = y1; (1)
  //  a + b + c = y2; (2)
  // 4(1)-(2) -> b + 3c = 4y1-y2
  const double b = 4*y1-y2-3*c;
  const double a = y2-b-c;

  // flat - either +ve or -ve infinity...
  if(a == 0){
    if(pred) *pred = y0;
    return 0;
  }

  // 2*a*x+b = 0 -> x = -b/2a

  const double x = -b/(2*a);

  if(pred) *pred = a*x*x+b*x+c;

  return x;
}

// ----------------------------------------------------------------------------
struct LiteHit
{
  LiteHit(double t, double a, double r) : time(t), amp(a), rms(r) {}

  double time;
  double amp;
  double rms;
};

// ----------------------------------------------------------------------------
// Sum all consecutive runs into their average position
std::vector<LiteHit> consolidate(TVectorDFast& alpha)
{
  std::vector<LiteHit> ret;

  double tot = 0;
  double tavg = 0;

  const int N = alpha.GetNrows();

  for(int i = 0; i <= N; ++i){
    // end of a run
    if((i == N || alpha[i] == 0) && tot > 0){
      ret.emplace_back(tavg/tot, tot, kRMS);
      tot = 0;
      tavg = 0;
    }

    if(i == N) break;

    tot += alpha[i];
    tavg += alpha[i]*i;
  }

  return ret;
}

// ----------------------------------------------------------------------------
std::vector<LiteHit> Fit(const std::vector<float>& yvec)
{
  const int N = yvec.size();
  TVectorDFast y(N);
  for(int i = 0; i < N; ++i) y[i] = yvec[i];

  const int range = 3*kRMS;
  const GausKernel kern(kRMS, range);

  TVectorDFast alpha(N);

  TVectorDFast pred(N);

  double chisq0 = chisq(y, pred, alpha);

  TVectorDFast diff_static(N);
  for(int k = 0; k < N; ++k){
    for(int i = std::max(0, k-range); i < std::min(N, k+range+1); ++i){
      diff_static[i] += -2*y[k]*kern.y(i-k);
    }
    diff_static[k] += L1_strength;
  }

  for(int trial = 0; ; ++trial){
    TVectorDFast diff = diff_static;

    // This is the hot loop
    for(int k = 0; k < N; ++k){
      // Suprisingly this is faster without checking pred[k] != 0 first
      for(int i = std::max(0, k-range); i < std::min(N, k+range+1); ++i){
        diff[i] += 2*pred[k]*kern.y(i-k);
      }
    }

    // Step downhill, but not out of the space
    TVectorDFast step = -1.*diff;
    for(int i = 0; i < N; ++i){
      if(alpha[i] == 0 && step[i] < 0){step[i] = 0;}
    }

    // Now we have a candidate direction

    double stepsize = std::numeric_limits<double>::infinity();;
    int closestIdx = -1;
    for(int i = 0; i < N; ++i){
      if(step[i] >= 0) continue; // stepping away from the wall
      if(alpha[i] + stepsize*step[i] < 0){
        const double distToWall = -alpha[i]/step[i];
        stepsize = distToWall;
        closestIdx = i;
      }
    }

    if(closestIdx == -1) stepsize = 1;

    TVectorDFast predStep(N);
    for(int i = 0; i < N; ++i){
      // Whereas this check does seem to be worthwhile
      if(step[i] != 0){
        for(int j = std::max(0, i-range); j < std::min(N, i+range+1); ++j){
          predStep[j] += step[i]*kern.y(i-j);
        }
      }
    }


    const TVectorDFast alpha1 = alpha + .5*stepsize*step;
    const TVectorDFast alpha2 = alpha +    stepsize*step;

    const TVectorDFast pred1 = pred + .5*stepsize*predStep;
    const TVectorDFast pred2 = pred +    stepsize*predStep;

    const double chisq1 = chisq(y, pred1, alpha1);
    const double chisq2 = chisq(y, pred2, alpha2);

    double predict;
    const double opt_step = find_min(chisq0, chisq1, chisq2, &predict)*stepsize;

    if(opt_step < stepsize){
      stepsize = opt_step;
      closestIdx = -1;

      if(opt_step == 0){
        std::cout << "Zero step. Done" << std::endl;
        break;
      }
      if(chisq0-predict < 1e-4){//1e-6){
        //        std::cout << "Tiny improvement expected. Done" << std::endl;
        break;
      }

      alpha += stepsize*step;
      pred += stepsize*predStep;
      chisq0 = chisq(y, pred, alpha);
    }
    else{
      alpha = alpha2;
      pred = pred2;
      chisq0 = chisq2;

      if(closestIdx >= 0) alpha[closestIdx] = 0; // make sure
    }
  } // end for trial

  return consolidate(alpha);
}



namespace hf
{
// ----------------------------------------------------------------------------
class CompressedHitFinder: public art::EDProducer
{
public:
  explicit CompressedHitFinder(const fhicl::ParameterSet& pset);

  void produce(art::Event& evt) override;
  void endJob() override;

protected:
  recob::Hit LiteHitToHit(const LiteHit& hit, const recob::Wire& wire, int t0) const;

  art::ServiceHandle<geo::Geometry> fGeom;

  std::string fWireLabel;
  double fThreshold;
};

DEFINE_ART_MODULE(CompressedHitFinder)

// ---------------------------------------------------------------------------
CompressedHitFinder::CompressedHitFinder(const fhicl::ParameterSet& pset)
: EDProducer(pset),
  fWireLabel(pset.get<std::string>("WireLabel")),
  fThreshold(pset.get<double>("Threshold"))
{
  produces<std::vector<recob::Hit>>();
}

// ----------------------------------------------------------------------------
void CompressedHitFinder::endJob()
{
}

// ----------------------------------------------------------------------------
void CompressedHitFinder::produce(art::Event& evt)
{
  auto hitcol = std::make_unique<std::vector<recob::Hit>>();

  art::Handle<std::vector<recob::Wire>> wires;
  evt.getByLabel(fWireLabel, wires);

  int wireIdx = 0;
  for(const recob::Wire& wire: *wires){
    for(const lar::sparse_vector<float>::datarange_t& range: wire.SignalROI().get_ranges()){
      for(const LiteHit& hit: Fit(range.data())){
        if(hit.amp < fThreshold) continue;

        hitcol->push_back(LiteHitToHit(hit, wire, range.begin_index()));
      } // end for hit
    } // end for range

    ++wireIdx;
    if(wireIdx%100 == 0) std::cout << wireIdx << " / " << wires->size() << std::endl;
  } // end for wire

  evt.put(std::move(hitcol));
}

// ----------------------------------------------------------------------------
recob::Hit CompressedHitFinder::LiteHitToHit(const LiteHit& hit,
                                             const recob::Wire& wire,
                                             int t0) const
{
  const double peaktime = hit.time + t0;
  const double integral = hit.amp;
  const double summedADC = hit.amp; // TODO - how is this supposed to be defined?
  const double peakamp =  hit.amp/(hit.rms*sqrt(2*M_PI));

  const geo::SigType_t sigType = fGeom->SignalType(wire.Channel());

  // Convention to always put it on the first one?
  const geo::WireID id = fGeom->ChannelToWire(wire.Channel())[0];

  return recob::Hit(wire.Channel(),
                    0/*starttick*/, 0/*endtick*/, peaktime,
                    0/*sigmaT*/, hit.rms, peakamp, 0/*sigmaA*/,
                    summedADC, integral, 0/*sigmaI*/,
                    1/*multiplicity*/, -1/*localIdx*/,
                    0/*GoF*/, 0/*DoF*/,
                    wire.View(), sigType, id);
}

} // namespace
