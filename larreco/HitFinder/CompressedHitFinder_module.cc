// Chris Backhouse - bckhouse@fnal.gov

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"

#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"

#include "larcore/Geometry/Geometry.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h" // Uncompress()

#include <iostream>

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompLU.h"

// Larger is faster
const double effective_threshold = 0.1;//1; // effectively this much is subtracted from each alpha. TODO seems to be more than that?
const double L1_strength = 2*effective_threshold;

const double kRMS = 2.75; // TODO
const double kRMSBiPolar = 5.5; // TODO

// ----------------------------------------------------------------------------
/// TMatrixD::operator() does various sanity checks and shows up in profiles
class TMatrixDFast: public TMatrixD
{
public:
  TMatrixDFast(int N, int M) : TMatrixD(N, M) {}

  // Skip bounds-checking overhead
  inline double  operator()(int i, int j) const {return fElements[fNrows*i+j];}
  inline double& operator()(int i, int j)       {return fElements[fNrows*i+j];}
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

TVectorDFast RowReductionSolver(TMatrixDFast M, TVectorDFast v)
{
  const int N = v.GetNrows();

  // Make unit diagonal. Not mathematically necessary, but seems to help
  // stability.
  for(int row = 0; row < N; ++row){
    const double factor = 1/M(row, row);
    for(int col = 0; col < N; ++col){
      M(row, col) *= factor;
    }
    v(row) *= factor;
  }

  // Clear out the region below the diagonal
  for(int row = 0; row < N; ++row){
    // We know our matrix in practice can have a lot of zeros before the
    // diagonal. The initialization of col0 can be moved outside of the loop
    // for a small speedup if we're additionally willing to rely on the fact
    // that in our matrices the length of that string of zeros is monotonically
    // increasing. Likewise for col1 (though the logic will have to change a
    // little).
    int col0 = 0;
    while(M(row, col0) == 0 && col0 < row) ++col0;
    int col1 = N-1;
    while(M(row, col1) == 0 && col1 >= 0) --col1;

    for(int col = col0; col < row; ++col){
      const double factor = -M(row, col);
      const int origin_row = col;
      M(row, col) = 0;
      // All columns outside this range are already zero in the origin row
      for(int i = col+1; i <= col1; ++i){
        M(row, i) += factor * M(origin_row, i);
      }
      v(row) += factor * v(origin_row);
    }

    // Make sure we have 1 on the diagonal
    const double div = 1/M(row, row);
    M(row, row) = 1;
    for(int i = row+1; i < N; ++i) M(row, i) *= div;
    v(row) *= div;
  }

  // Now clear out the upper-right region
  /*
  for(int row = N-1; row >= 0; --row){
    for(int col = row+1; col < N; ++col){
      const int origin_row = col;
      const double factor = -M(row, col);
      for(int i = 0; i < N; ++i) M(row, i) += factor * M(origin_row, i);
      v(row) += factor * v(origin_row);
    }
  }
  */

  // Which is equivalent to just this
  for(int i = N; i >= 0; --i)
    for(int j = i+1; j < N; ++j)
      v(i) -= M(i, j) * v(j);
  // Could also have made use of the trailing zeros in the rows here for
  // another small speedup.

  return v;
}

// ----------------------------------------------------------------------------
class Kernel
{
public:
  static Kernel Gaus(double rms, int range)
  {
    Kernel ret(rms, range);

    const double norm = 1/(rms*sqrt(2*M_PI));
    for(int i = 0; i < 2*range+1; ++i){
      const double x = (i-range)/rms;
      ret.fData[i] = norm*exp(-x*x/2);
    }

    return ret;
  }

  static Kernel BiPolar(double rms, int range)
  {
    Kernel ret(rms, range);

    const double norm = 1/(2*rms); // makes the total absolute area = 1
    for(int i = 0; i < 2*range+1; ++i){
      const double x = (i-range)/rms;
      ret.fData[i] = -norm*x*exp(-x*x/2);
    }

    return ret;
  }

  inline double operator()(int dx) const {return fData[fRange+dx];}
  inline double operator[](int dx) const {return fData[fRange+dx];}

  inline int Range() const {return fRange;}
  double RMS() const {return fRMS;}

protected:
  Kernel(double rms, int range) : fRMS(rms), fRange(range), fData(2*range+1) {}

  double fRMS;
  int fRange;
  std::vector<double> fData;
};

// Find the exact minimum of the chisq function for any subset of variables
class ExactSolver
{
public:
  ExactSolver(const Kernel& kern, const std::vector<float>& y)
    : A(y.size(), y.size()), B(y.size())
  {
    const int N = y.size();
    const int range = kern.Range();

    // This matrix and vector give the solution to the full problem with all
    // variables free (which probably falls outside of the physical region).

    // TODO this thing has such a simple structure we probably don't need to
    // store it like this but can eg just keep one slice through it.
    for(int k = 0; k < N; ++k){
      for(int i = std::max(0, k-range); i < std::min(N, k+range+1); ++i){
        for(int j = std::max(0, i-range); j < std::min(N, i+range+1); ++j){
          A(k, j) += kern(i-j)*kern(i-k);
        }
        B(k) += kern(i-k)*y[i];
      }
      B(k) -= L1_strength/2;
    }
  }

  // Solve the reduced problem with only indices 'idxs' free, and all other
  // variables held fixed at zero.
  void Solve(const std::vector<int>& idxs, TVectorDFast& ret) const
  {
    const int N_sub = idxs.size();

    if(N_sub == 0) return;

    TMatrixDFast A_sub(N_sub, N_sub);
    TVectorDFast B_sub(N_sub);

    for(int i = 0; i < N_sub; ++i){
      for(int j = 0; j < N_sub; ++j){
        A_sub(i, j) = A(idxs[i], idxs[j]);
      }
      B_sub[i] = B[idxs[i]];
    }

    //    TDecompLU d(A_sub);
    //    bool ok;
    //    const TVectorDFast proposed = d.Solve(B_sub, ok);

    const TVectorDFast proposed = RowReductionSolver(A_sub, B_sub);

    for(int i = 0; i < N_sub; ++i) ret[idxs[i]] = proposed[i];
  }

protected:
  TMatrixDFast A;
  TVectorDFast B;
};

// ----------------------------------------------------------------------------
struct LiteHit
{
  LiteHit(double t, double a, double r) : time(t), amp(a), rms(r) {}

  double time;
  double amp;
  double rms;
};

// ----------------------------------------------------------------------------
// Search through the array of amplitudes, and consolidate each consecutive
// runs of non-zero amplitudes into its mean and total magnitude
std::vector<LiteHit> consolidate(const TVectorDFast& alpha, const Kernel& kern)
{
  std::vector<LiteHit> ret;

  double tot = 0;
  double tavg = 0;

  const int N = alpha.GetNrows();

  for(int i = 0; i <= N; ++i){
    // end of a run
    if((i == N || alpha[i] == 0) && tot > 0){
      ret.emplace_back(tavg/tot, tot, kern.RMS());
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
void take_step(const TVectorDFast& step,
               TVectorDFast& alpha,
               TVectorDFast& pred,
               const Kernel& kern)
{
  alpha += step;

  // Update the predictions
  const int N = step.GetNrows();
  const int range = kern.Range();

  for(int j = 0; j < N; ++j){
    // This check does seem to be worthwhile
    if(step[j] != 0){
      for(int i = std::max(0, j-range); i < std::min(N, j+range+1); ++i){
        pred[i] += step[j]*kern(i-j);
      }
    }
  }
}

// ----------------------------------------------------------------------------
// Updates 'alpha' by 'step', unless that position would have a negative
// coordinate. In that case , update to the first position in direction 'step'
// where one of the coordinates becomes zero. Adjusts 'pred'
// correspondingly. Returns whether the step was indeed to a wall.
bool project_to_wall(const TVectorDFast& step,
                     TVectorDFast& alpha,
                     TVectorDFast& pred,
                     const Kernel& kern)
{
  const int N = step.GetNrows();

  double stepsize = 1; // try to take the full step
  int closestIdx = -1;
  for(int i = 0; i < N; ++i){
    // Stepping towards the wall, and our current smallest step would take us
    // through it.
    if(step[i] < 0 && alpha[i] + stepsize*step[i] < 0){
      stepsize = -alpha[i]/step[i];
      closestIdx = i;
    }
  }

  take_step(stepsize*step, alpha, pred, kern);

  if(closestIdx >= 0) alpha[closestIdx] = 0; // make sure to touch wall exactly

  return closestIdx >= 0;
}

// ----------------------------------------------------------------------------
std::vector<LiteHit> Fit(const std::vector<float>& yvec, const Kernel& kern)
{
  // We follow the algorithm from "Efficient sparse coding algorithms", H. Lee
  // et al.
  const int N = yvec.size();
  const int range = kern.Range();

  // Convert the input data to TVector format
  TVectorDFast y(N);
  for(int i = 0; i < N; ++i) y[i] = yvec[i];

  // Solver that we will be using throughout
  const ExactSolver exact(kern, yvec);

  // Amplitudes - initialize all zeros
  TVectorDFast alpha(N);

  // Prediction implied by those amplitudes - ditto
  TVectorDFast pred(N);

  // Those parts of the derivative dchisq/dalpha that are independent of alpha
  TVectorDFast diff_static(N);
  for(int k = 0; k < N; ++k){
    for(int i = std::max(0, k-range); i < std::min(N, k+range+1); ++i){
      diff_static[k] += -2*y[i]*kern(i-k);
    }
    diff_static[k] += L1_strength;
  }

  for(int trial = 0; ; ++trial){
    // Update the derivative by adding the alpha-depending terms
    TVectorDFast diff = diff_static;

    // This is the hottest loop in the fit
    for(int k = 0; k < N; ++k){
      // Suprisingly this is faster without checking pred[i] != 0 first
      for(int i = std::max(0, k-range); i < std::min(N, k+range+1); ++i){
        diff[k] += 2*pred[i]*kern(i-k);
      }
    }

    // Search for the variable at a physical boundary that has the largest
    // derivative downwards into the space.
    double bestDiff = 0;
    int bestIdx = -1;
    for(int i = 0; i < N; ++i){
      if(alpha[i] == 0 && diff[i] < bestDiff){
        bestDiff = diff[i];
        bestIdx = i;
      }
    }

    // There are no variables left that it would be advantageous to add, so
    // we're done
    if(bestIdx == -1) break;

    while(true){
      // The active set is all those variables not on boundaries, plus, in the
      // first loop, the one most powerful variable from above.
      std::vector<int> idxs;
      for(int i = 0; i < N; ++i) if(alpha[i] > 0) idxs.push_back(i);
      if(bestIdx >= 0){idxs.push_back(bestIdx); bestIdx = -1;}

      // Where the exact solution would want us to move
      TVectorDFast proposed(N);
      exact.Solve(idxs, proposed);

      // Do the subtraction this way to be certain we get exact zeros where we
      // expect.
      TVectorDFast step(N);
      for(int i: idxs) step[i] = proposed[i]-alpha[i];

      // Step as far as we can towards the exact minimum. If we hit a wall we
      // remove that variable and go round again. If we don't hit any walls
      // we'll break out and try to add a new variable.
      if(!project_to_wall(step, alpha, pred, kern)) break;
    }
  } // end for trial

  // Turn the alphas into actual hits
  return consolidate(alpha, kern);
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

  bool fFitRawDigits;
  std::string fWireLabel;
  std::string fRawDigitLabel;
  double fThreshold;
};

DEFINE_ART_MODULE(CompressedHitFinder)

// ---------------------------------------------------------------------------
CompressedHitFinder::CompressedHitFinder(const fhicl::ParameterSet& pset)
: EDProducer(pset),
  fFitRawDigits(pset.get<bool>("FitRawDigits")),
  fWireLabel(pset.get<std::string>("WireLabel")),
  fRawDigitLabel(pset.get<std::string>("RawDigitLabel")),
  fThreshold(pset.get<double>("Threshold"))
{
  produces<std::vector<recob::Hit>>();

  if(fFitRawDigits)
    produces<art::Assns<recob::Hit, raw::RawDigit>>();
  else
    produces<art::Assns<recob::Hit, recob::Wire>>();
}

// ----------------------------------------------------------------------------
void CompressedHitFinder::endJob()
{
}

// ----------------------------------------------------------------------------
void CompressedHitFinder::produce(art::Event& evt)
{
  auto hitcol = std::make_unique<std::vector<recob::Hit>>();

  auto assns_dig  = std::make_unique<art::Assns<recob::Hit, raw::RawDigit>>();
  auto assns_wire = std::make_unique<art::Assns<recob::Hit, recob::Wire>>();

  art::PtrMaker<recob::Hit> pmaker(evt);


  art::Handle<std::vector<recob::Wire>> wires;
  evt.getByLabel(fWireLabel, wires);


  std::map<geo::WireID, const art::Ptr<raw::RawDigit>> digmap;

  art::Handle<std::vector<raw::RawDigit>> digs;

  if(fFitRawDigits){
    evt.getByLabel(fRawDigitLabel, digs);

    for(unsigned int digIdx = 0; digIdx < digs->size(); ++digIdx){
      const art::Ptr<raw::RawDigit> dig(digs, digIdx);
      digmap.emplace(fGeom->ChannelToWire(dig->Channel())[0], dig);
    }
  }


  for(unsigned int wireIdx = 0; wireIdx < wires->size(); ++wireIdx){
    const recob::Wire& wire = (*wires)[wireIdx];

    const geo::SigType_t sigType = fGeom->SignalType(wire.Channel());

    const Kernel kern = (fFitRawDigits && sigType == geo::kInduction) ? Kernel::BiPolar(kRMSBiPolar, 3*kRMSBiPolar) : Kernel::Gaus(kRMS, 3*kRMS);

    const geo::WireID id = fGeom->ChannelToWire(wire.Channel())[0];
    const raw::RawDigit* dig = digmap[id].get();

    for(const lar::sparse_vector<float>::datarange_t& range: wire.SignalROI().get_ranges()){
      std::vector<LiteHit> hits;
      if(!fFitRawDigits){
        hits = Fit(range.data(), kern);
      }
      else{
        std::vector<short> adcs;
        if(dig) raw::Uncompress(dig->ADCs(), adcs, dig->Compression());

        std::vector<float> tofit(std::min(adcs.end(), adcs.begin()+range.begin_index()),
                                std::min(adcs.end(), adcs.begin()+range.end_index()+1));

        for(float& y: tofit) if(y != 0) y -= dig->GetPedestal();

        hits = Fit(tofit, kern);
      }

      for(const LiteHit& hit: hits){
        if(hit.amp < fThreshold) continue;

        hitcol->push_back(LiteHitToHit(hit, wire, range.begin_index()));

        // Synthesize a Ptr for the hit we just pushed
        art::Ptr<recob::Hit> phit = pmaker(hitcol->size()-1);

        if(fFitRawDigits){
          assns_dig->addSingle(phit, digmap[id]);
        }
        else{
          assns_wire->addSingle(phit, art::Ptr<recob::Wire>(wires, wireIdx));
        }
      } // end for hit
    } // end for range

    if(wireIdx%100 == 0) std::cout << wireIdx << " / " << wires->size() << std::endl;
  } // end for wire

  evt.put(std::move(hitcol));

  if(fFitRawDigits) evt.put(std::move(assns_dig)); else evt.put(std::move(assns_wire));
}

// ----------------------------------------------------------------------------
recob::Hit CompressedHitFinder::LiteHitToHit(const LiteHit& hit,
                                             const recob::Wire& wire,
                                             int t0) const
{
  const double peaktime = hit.time + t0;
  const double integral = hit.amp;
  const double summedADC = hit.amp; // TODO - how is this supposed to be defined?
  // TODO what formula to use for induction hits?
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
