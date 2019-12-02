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

//#include "TDecompLU.h"
//#include "TDecompSVD.h" // much slower
#include "TDecompSparse.h" // faster
#include "TMatrixD.h"
#include "TVectorD.h"

#include "TMatrixDSparse.h"

// Larger is faster
const double effective_threshold = 0.1;//1; // effectively this much is subtracted from each alpha. TODO seems to be more than that?
const double L1_strength = 2*effective_threshold;

const double kRMS = 2.75; // TODO

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

  inline double operator()(int dx) const {return fData[fRange+dx];}
  inline double operator[](int dx) const {return fData[fRange+dx];}

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
// Find the minimum of a parabola, given three points sampled from it.
//
// y0, y1, y2 should be evenly spaced, and the result is given as a fraction of
// the distance from y0 to y2.
double find_min(double y0, double y1, double y2, double* pred = 0)
{
  const double c = y0;
  //  a/4 + b/2 + c = y1 (1)
  //  a   + b   + c = y2 (2)
  // 4(1)-(2) -> b + 3c = 4y1 - y2
  const double b =  -y2 + 4*y1 - 3*y0;
  //  const double a = y2 - b - c;
  const double a = 2*y2 - 4*y1 + 2*y0;

  // flat - either +ve or -ve infinity...
  if(a == 0){
    if(pred) *pred = y0;
    return 0;
  }

  // 2ax + b = 0 -> x = -b/2a

  const double x = -b/(2*a);

  if(pred) *pred = a*x*x + b*x + c;

  return x;
}

// Find the exact minimum of the chisq function for any subset of variables
class ExactSolver
{
public:
  ExactSolver(const GausKernel& kern, const std::vector<float>& y)
    : A(y.size(), y.size()), B(y.size())
  {
    const int N = y.size();
    const int range = kern.Range();

    // This matrix and vector give the solution to the full problem with all
    // variables free (which probably falls outside of the physical region).

    for(int i = 0; i < N; ++i){
      for(int j = std::max(0, i-range); j < std::min(N, i+range+1); ++j){
        for(int k = std::max(0, i-range); k < std::min(N, i+range+1); ++k){
          A(j, k) += kern(i-j)*kern(i-k);
        }
        B(i) += kern(i-j)*y[j];
      }
      B(i) -= L1_strength/2;
    }
  }

  // Solve the reduced problem with only indices 'idxs' free, and all other
  // variables held fixed at zero. Returns false if something went wrong, in
  // that case 'ret' can't be trusted.
  bool Solve(const std::vector<int>& idxs, TVectorDFast& ret) const
  {
    const int N_sub = idxs.size();

    if(N_sub == 0) return false;

    // Surprisingly it's faster to accumulate this densely and sparsen after
    TMatrixDFast A_sub(N_sub, N_sub);
    TVectorDFast B_sub(N_sub);

    for(int i = 0; i < N_sub; ++i){
      for(int j = 0; j < N_sub; ++j){
        A_sub(i, j) = A(idxs[i], idxs[j]);
      }
      B_sub[i] = B[idxs[i]];
    }

    const TMatrixDSparse sparse(A_sub);
    TDecompSparse decomp(sparse, 0/*verbosity*/);
    bool ok;
    const TVectorDFast proposed = decomp.Solve(B_sub, ok);
    if(!ok) return false;

    for(int i = 0; i < N_sub; ++i) ret[idxs[i]] = proposed[i];

    return true;
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
std::vector<LiteHit> consolidate(const TVectorDFast& alpha)
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
// Updates 'alpha' to the first position in direction 'step' where one of the
// coordinates becomes zero. Adjusts 'pred' correspondingly. If there is no
// such intercept, steps an arbitrary distance. Returns whether the step was
// indeed to a wall.
bool project_to_wall(const TVectorDFast& step,
                     TVectorDFast& alpha,
                     TVectorDFast& pred,
                     const GausKernel& kern)
{
  const int N = step.GetNrows();

  double stepsize = std::numeric_limits<double>::infinity();;
  int closestIdx = -1;
  for(int i = 0; i < N; ++i){
    // Stepping towards the wall, and our current smallest step would take us
    // through it.
    if(step[i] < 0 && alpha[i] + stepsize*step[i] < 0){
      stepsize = -alpha[i]/step[i];
      closestIdx = i;
    }
  }

  // Motion away from all walls, just step by the original step size
  if(closestIdx == -1) stepsize = 1;

  alpha += stepsize*step;

  if(closestIdx >= 0) alpha[closestIdx] = 0; // make sure to touch wall exactly

  const int range = kern.Range();

  for(int i = 0; i < N; ++i){
    // This check does seem to be worthwhile
    if(step[i] != 0){
      for(int j = std::max(0, i-range); j < std::min(N, i+range+1); ++j){
        pred[j] += stepsize*step[i]*kern(i-j);
      }
    }
  }

  return closestIdx >= 0;
}

// ----------------------------------------------------------------------------
std::vector<LiteHit> Fit(const std::vector<float>& yvec)
{
  const int N = yvec.size();
  // Convert the input data to TVector format
  TVectorDFast y(N);
  for(int i = 0; i < N; ++i) y[i] = yvec[i];

  const int range = 3*kRMS;
  const GausKernel kern(kRMS, range);

  // Solver that we will be using throughout
  const ExactSolver exact(kern, yvec);

  // Amplitudes - initialize all zeros
  TVectorDFast alpha(N);

  // Prediction implied by those amplitudes - ditto
  TVectorDFast pred(N);

  // The chisq at the current position
  double chisq0 = chisq(y, pred, alpha);

  // Those parts of the derivative dchisq/dalpha that are independent of alpha
  TVectorDFast diff_static(N);
  for(int k = 0; k < N; ++k){
    for(int i = std::max(0, k-range); i < std::min(N, k+range+1); ++i){
      diff_static[i] += -2*y[k]*kern(i-k);
    }
    diff_static[k] += L1_strength;
  }

  std::vector<int> prevIdxs;
  prevIdxs.reserve(N);

  for(int trial = 0; ; ++trial){
    // Update the derivative by adding the alpha-depending terms
    TVectorDFast diff = diff_static;

    // This is the hottest loop in the fit
    for(int k = 0; k < N; ++k){
      // Suprisingly this is faster without checking pred[k] != 0 first
      for(int i = std::max(0, k-range); i < std::min(N, k+range+1); ++i){
        diff[i] += 2*pred[k]*kern(i-k);
      }
    }

    std::vector<int> idxs;
    idxs.reserve(N);

    // Proposed step for gradient descent. Step downhill, but not if we're
    // against a wall and the step is out of the space.
    TVectorDFast step = -1.*diff;
    for(int i = 0; i < N; ++i){
      if(alpha[i] == 0 && step[i] < 0){
        step[i] = 0;
      }
      else{
        // All good, so this variable will be free in the fit
        idxs.push_back(i);
      }
    }

    // If we have a different set of variables fixed to zero than before, then
    // it is worth trying the exact fit again.
    if(idxs != prevIdxs){
      prevIdxs = idxs;

      const int N_sub = idxs.size();

      if(N_sub == 0) return {}; // no hits at all

      TVectorDFast proposed(N);
      bool ok = exact.Solve(idxs, proposed);
      if(!ok) continue; // something went wrong, don't step

      // If the exact solution is inside the space then this will be the
      // ultimate solution.
      for(int i = 0; i < N; ++i){
        if(proposed[i] < 0) ok = false;
      }

      if(ok){
        std::cout << "Exact solution after " << trial
                  << " with " << N_sub << std::endl;

        alpha = proposed;
        break;
      }

      // Do the subtraction this way to be certain we get exact zeros where we
      // expect.
      TVectorDFast exactStep(N);
      for(int i: idxs) exactStep[i] = proposed[i]-alpha[i];

      // Step as far as we can towards the exact minimum (which is outside the
      // space), and go around the loop again (now with a new set of idxs).
      project_to_wall(exactStep, alpha, pred, kern);
      chisq0 = chisq(y, pred, alpha);
      continue;
    }

    // We didn't do the exact solution thing, so we should try to take our best
    // step downhill instead.

    // How much would the prediction change if we took the full step?
    TVectorDFast predStep(N);
    for(int i = 0; i < N; ++i){
      // Whereas this check does seem to be worthwhile
      if(step[i] != 0){
        for(int j = std::max(0, i-range); j < std::min(N, i+range+1); ++j){
          predStep[j] += step[i]*kern(i-j);
        }
      }
    }

    // Step all the way to the wall (or random distance otherwise)
    TVectorDFast alpha2 = alpha;
    TVectorDFast pred2 = pred;
    const bool toWall = project_to_wall(step, alpha2, pred2, kern);

    // Also step half that distance
    const TVectorDFast alpha1 = .5*(alpha + alpha2);
    const TVectorDFast pred1 = .5*(pred + pred2);

    // And then we'll interpolate the chisq and figure out the optimum step
    const double chisq1 = chisq(y, pred1, alpha1);
    const double chisq2 = chisq(y, pred2, alpha2);
    const double opt_step = find_min(chisq0, chisq1, chisq2);

    // If the optimum step is into open space, or not as far as the wall, take
    // it
    if(!toWall || opt_step < 1){
      if(opt_step <= 0){ // TODO where do slightly negative values come from?
        std::cout << "Zero step. Done" << std::endl;
        break;
      }

      alpha = opt_step*alpha2 + (1-opt_step)*alpha;
      pred  = opt_step*pred2  + (1-opt_step)*pred;
      chisq0 = chisq(y, pred, alpha);
    }
    else{
      // Otherwise take the full step to the wall
      alpha = alpha2;
      pred = pred2;
      chisq0 = chisq2;
    }
  } // end for trial

  // Turn the alphas into actual hits
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
