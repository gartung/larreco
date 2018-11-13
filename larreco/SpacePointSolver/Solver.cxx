#include "Solver.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <iostream>

#include "TMatrixD.h"
#include "TMatrixDSparse.h"
#include "TMatrixDSym.h"
#include "TMatrixF.h"
#include "TMatrixFSym.h"
#include "TVectorD.h"
#include "TVectorF.h"
#include "TDecompSparse.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
#include "TDecompChol.h" // symmetric
#include "TDecompBK.h" // symmetric
#include "TDecompQRH.h"
//#define TDecompSparse TDecompLU
//#define TMatrixDSparse TMatrixD

template<class T> T sqr(T x){return x*x;}

// ---------------------------------------------------------------------------
InductionWireHit::InductionWireHit(int chan, double q)
  : fChannel(chan), fCharge(q), fPred(0)
{
}

// ---------------------------------------------------------------------------
Neighbour::Neighbour(SpaceCharge* sc, double coupling)
  : fSC(sc), fCoupling(coupling)
{
}

// ---------------------------------------------------------------------------
SpaceCharge::SpaceCharge(double x, double y, double z,
                         CollectionWireHit* cwire,
                         InductionWireHit* wire1, InductionWireHit* wire2)
  : fX(x), fY(y), fZ(z),
    fCWire(cwire), fWire1(wire1), fWire2(wire2),
    fPred(0),
    fNeiPotential(0)
{
}

// ---------------------------------------------------------------------------
void SpaceCharge::AddCharge(double dq)
{
  fPred += dq;

  for(Neighbour& nei: fNeighbours)
    nei.fSC->fNeiPotential += dq * nei.fCoupling;

  if(fWire1) fWire1->fPred += dq;
  if(fWire2) fWire2->fPred += dq;
}

// ---------------------------------------------------------------------------
CollectionWireHit::CollectionWireHit(int chan, double q,
                                     const std::vector<SpaceCharge*>& cross)
  : fChannel(chan), fCharge(q), fCrossings(cross)
{
  if(q < 0){
    std::cout << "Trying to construct collection wire with negative charge " << q << " this should never happen." << std::endl;
    abort();
  }

  const double p = q/cross.size();

  for(SpaceCharge* sc: cross) sc->AddCharge(p);
}

// ---------------------------------------------------------------------------
CollectionWireHit::~CollectionWireHit()
{
  // Each SpaceCharge can only be in one collection wire, so we own them. But
  // not SpaceCharge does not clean up its induction wires since they're
  // shared.
  for(SpaceCharge* sc: fCrossings) delete sc;
}

// ---------------------------------------------------------------------------
double Metric(double q, double p)
{
  return sqr(q-p);
}

// ---------------------------------------------------------------------------
QuadExpr Metric(double q, QuadExpr p)
{
  return sqr(q-p);
}

// ---------------------------------------------------------------------------
double Metric(const std::vector<SpaceCharge*>& scs, double alpha)
{
  double ret = 0;

  std::set<InductionWireHit*> iwires;
  for(const SpaceCharge* sc: scs){
    if(sc->fWire1) iwires.insert(sc->fWire1);
    if(sc->fWire2) iwires.insert(sc->fWire2);

    if(alpha != 0){
      ret -= alpha*sqr(sc->fPred);
      // "Double-counting" of the two ends of the connection is
      // intentional. Otherwise we'd have a half in the line above.
      ret -= alpha * sc->fPred * sc->fNeiPotential;
    }
  }

  for(const InductionWireHit* iwire: iwires){
    ret += Metric(iwire->fCharge, iwire->fPred);
  }

  return ret;
}

// ---------------------------------------------------------------------------
double Metric(const std::vector<CollectionWireHit*>& cwires, double alpha)
{
  std::vector<SpaceCharge*> scs;
  for(CollectionWireHit* cwire: cwires)
    scs.insert(scs.end(), cwire->fCrossings.begin(), cwire->fCrossings.end());
  return Metric(scs, alpha);
}

// ---------------------------------------------------------------------------
QuadExpr Metric(const SpaceCharge* sci, const SpaceCharge* scj, double alpha)
{
  QuadExpr ret = 0;

  // How much charge moves from scj to sci
  QuadExpr x = QuadExpr::X();

  if(alpha != 0){
    const double scip = sci->fPred;
    const double scjp = scj->fPred;

    // Self energy. SpaceCharges are never the same object
    ret -= alpha*sqr(scip + x);
    ret -= alpha*sqr(scjp - x);

    // Interaction. We're only seeing one end of the double-ended connection
    // here, so multiply by two.
    ret -= 2 * alpha * (scip + x) * sci->fNeiPotential;
    ret -= 2 * alpha * (scjp - x) * scj->fNeiPotential;

    // This miscounts if i and j are neighbours of each other
    for(const Neighbour& n: sci->fNeighbours){
      if(n.fSC == scj){
        // If we detect that case, remove the erroneous terms
        ret += 2 * alpha * (scip + x) * scjp * n.fCoupling;
        ret += 2 * alpha * (scjp - x) * scip * n.fCoupling;

        // And replace with the correct interaction terms
        ret -= 2 * alpha * (scip + x) * (scjp - x) * n.fCoupling;
        break;
      }
    }
  }


  const InductionWireHit* iwire1 = sci->fWire1;
  const InductionWireHit* jwire1 = scj->fWire1;

  const double qi1 = iwire1 ? iwire1->fCharge : 0;
  const double pi1 = iwire1 ? iwire1->fPred : 0;

  const double qj1 = jwire1 ? jwire1->fCharge : 0;
  const double pj1 = jwire1 ? jwire1->fPred : 0;

  if(iwire1 == jwire1){
    // Same wire means movement of charge cancels itself out
    if(iwire1) ret += Metric(qi1, pi1);
  }
  else{
    if(iwire1) ret += Metric(qi1, pi1 + x);
    if(jwire1) ret += Metric(qj1, pj1 - x);
  }

  const InductionWireHit* iwire2 = sci->fWire2;
  const InductionWireHit* jwire2 = scj->fWire2;

  const double qi2 = iwire2 ? iwire2->fCharge : 0;
  const double pi2 = iwire2 ? iwire2->fPred : 0;

  const double qj2 = jwire2 ? jwire2->fCharge : 0;
  const double pj2 = jwire2 ? jwire2->fPred : 0;

  if(iwire2 == jwire2){
    if(iwire2) ret += Metric(qi2, pi2);
  }
  else{
    if(iwire2) ret += Metric(qi2, pi2 + x);
    if(jwire2) ret += Metric(qj2, pj2 - x);
  }

  return ret;
}

// ---------------------------------------------------------------------------
QuadExpr Metric(const SpaceCharge* sc, double alpha)
{
  QuadExpr ret = 0;

  // How much charge is added to sc
  QuadExpr x = QuadExpr::X();

  if(alpha != 0){
    const double scp = sc->fPred;

    // Self energy
    ret -= alpha*sqr(scp + x);

    // Interaction. We're only seeing one end of the double-ended connection
    // here, so multiply by two.
    ret -= 2 * alpha * (scp + x) * sc->fNeiPotential;
  }

  // Prediction of the induction wires
  ret += Metric(sc->fWire1->fCharge, sc->fWire1->fPred + x);
  ret += Metric(sc->fWire2->fCharge, sc->fWire2->fPred + x);

  return ret;
}

// ---------------------------------------------------------------------------
double SolvePair(CollectionWireHit* cwire,
                 SpaceCharge* sci, SpaceCharge* scj,
                 double alpha)
{
  const QuadExpr chisq = Metric(sci, scj, alpha);
  const double chisq0 = chisq.Eval(0);

  // Find the minimum of a quadratic expression
  double x = -chisq.Linear()/(2*chisq.Quadratic());

  // Don't allow either SpaceCharge to go negative
  const double xmin = -sci->fPred;
  const double xmax =  scj->fPred;

  // Clamp to allowed range
  x = std::min(xmax, x);
  x = std::max(xmin, x);

  const double chisq_new = chisq.Eval(x);

  // Should try these too, because the function might be convex not concave, so
  // d/dx=0 gives the max not the min, and the true min is at one extreme of
  // the range.
  const double chisq_p = chisq.Eval(xmax);
  const double chisq_n = chisq.Eval(xmin);

  if(std::min(std::min(chisq_p, chisq_n), chisq_new) > chisq0+1){
    std::cout << "Solution at " << x << " is worse than current state! Scan from " << xmin << " to " << xmax << std::endl;
    for(double x = xmin; x < xmax; x += .01*(xmax-xmin)){
      std::cout << x << " " << chisq.Eval(x) << std::endl;
    }

    std::cout << "Soln, original, up edge, low edge:" << std::endl;
    std::cout << chisq_new << " " << chisq0 << " " << chisq_p << " " << chisq_n << std::endl;
    abort();
  }

  if(std::min(chisq_n, chisq_p) < chisq_new){
    if(chisq_n < chisq_p) return xmin;
    return xmax;
  }

  return x;
}

// ---------------------------------------------------------------------------
void Iterate(CollectionWireHit* cwire, double alpha)
{
  // Consider all pairs of crossings
  const unsigned int N = cwire->fCrossings.size();

  //  std::cout << N << std::endl;

  for(unsigned int i = 0; i+1 < N; ++i){
    SpaceCharge* sci = cwire->fCrossings[i];

    for(unsigned int j = i+1; j < N; ++j){
      SpaceCharge* scj = cwire->fCrossings[j];

      const double x = SolvePair(cwire, sci, scj, alpha);

      if(x == 0) continue;

      // Actually make the update
      sci->AddCharge(+x);
      scj->AddCharge(-x);
    } // end for j
  } // end for i
}

// ---------------------------------------------------------------------------
template <class MatrixT, class DecompT, class VectorT>
class IterateQuadProg
{
public:
  static void Iterate(CollectionWireHit* cwire, double alpha)
  {
    //    std::unordered_set<SpaceCharge*> mask;

    bool printing = false;

    DecompT dc;

    // Clear out all the indices so that >= 0 means we're using it
    for(SpaceCharge* sc: cwire->fCrossings){
      //      if(sc->fWire1) sc->fWire1->fScratch = -1;
      //      if(sc->fWire2) sc->fWire2->fScratch = -1;
      for(Neighbour& n: sc->fNeighbours){
        n.fSC->fScratch = -1;
      }
    }

    for(SpaceCharge* sc: cwire->fCrossings){
      sc->fScratch = 0; // activate them, will be numbered later
      // Zero out all the space charges to get blank slate induction wire
      // predictions from the charges we're not optimizing
      sc->AddCharge(-sc->fPred);
    }

    while(true){
      std::vector<SpaceCharge*> scsPruned;
      scsPruned.reserve(cwire->fCrossings.size());
      for(SpaceCharge* sc: cwire->fCrossings){
        //        if(mask.count(sc) == 0)
        if(sc->fScratch >= 0){
          sc->fScratch = scsPruned.size();
          scsPruned.push_back(sc);
        }
      }

      const unsigned int Nscs = scsPruned.size();
      if(Nscs < 2) break;

      // Zero out all the space charges to get blank slate induction wire
      // predictions from the charges we're not optimizing
      //      for(SpaceCharge* sc: scsPruned) sc->AddCharge(-sc->fPred);

      //      std::unordered_map<InductionWireHit*, int> iwIdx;
      //      std::unordered_map<SpaceCharge*, int> neiIdx;

      // Need to re-clear the induction IDs so we can re-assign them
      // for(SpaceCharge* sc: scsPruned){
      //   if(sc->fWire1) sc->fWire1->fScratch = -1;
      //   if(sc->fWire2) sc->fWire2->fScratch = -1;
      // }

      // unsigned int nextIndID = 0;
      // //      for(unsigned int scIdx = 0; scIdx < Nscs; ++scIdx){
      // //        SpaceCharge* sc = scsPruned[scIdx];
      // for(SpaceCharge* sc: scsPruned){

      //   //        sc->fScratch = scIdx;
      //   //        neiIdx[scsPruned[scIdx]] = scIdx;

      //   if(sc->fWire1 && sc->fWire1->fScratch < 0) sc->fWire1->fScratch = nextIndID++;
      //   if(sc->fWire2 && sc->fWire2->fScratch < 0) sc->fWire2->fScratch = nextIndID++;

      //   // if(sc->fWire1 && iwIdx.find(sc->fWire1) == iwIdx.end()){
      //   //   iwIdx[sc->fWire1] = nextIndID;
      //   //   ++nextIndID;
      //   // }
      //   // if(sc->fWire2 && iwIdx.find(sc->fWire2) == iwIdx.end()){
      //   //   iwIdx[sc->fWire2] = nextIndID;
      //   //   ++nextIndID;
      //   // }
      // }

      //  std::cout << Nscs << " SCs" << std::endl;
      //  std::cout << iwIdx.size() << " inductions" << std::endl;

      //      const unsigned int Nind = nextIndID;//iwIdx.size();

      //      MatrixT T(Nind, Nscs+1); // Induction/spacepoint transfer matrix
      //      VectorT D(Nind); // induction wire residuals

      VectorT mc(Nscs+1); // minus c vector
      for(SpaceCharge* sc: scsPruned){
        //        if(sc->fWire1 && sc->fWire1->fScratch >= 0) D(sc->fWire1->fScratch) = sc->fWire1->fCharge - sc->fWire1->fPred;
        //        if(sc->fWire2 && sc->fWire2->fScratch >= 0) D(sc->fWire2->fScratch) = sc->fWire2->fCharge - sc->fWire2->fPred;

        if(sc->fWire1) mc(sc->fScratch) += sc->fWire1->fCharge - sc->fWire1->fPred;
        if(sc->fWire2) mc(sc->fScratch) += sc->fWire2->fCharge - sc->fWire2->fPred;
      }

      //      for(auto it: iwIdx) D(it.second) = it.first->fCharge - it.first->fPred;

      //      for(unsigned int scIdx = 0; scIdx < Nscs; ++scIdx){
      //        const SpaceCharge* sc = scsPruned[scIdx];
      //      for(const SpaceCharge* sc: scsPruned){

        //        if(sc->fWire1) T(sc->fWire1->fScratch, sc->fScratch) = 1;
        //        if(sc->fWire2) T(sc->fWire2->fScratch, sc->fScratch) = 1;

        //        if(sc->fWire1) T(iwIdx[sc->fWire1], scIdx) = 1;
        //        if(sc->fWire2) T(iwIdx[sc->fWire2], scIdx) = 1;
      //      }

      //  T.Print();
      //  D.Print();

      //      MatrixT TT(MatrixT::kTransposed, T);
      //      MatrixT Q = TT*T;

      TMatrixDSym Q(Nscs+1);
      //      MatrixT Q(Nscs+1, Nscs+1);
      // TODO - would be faster as a loop over induction wires?
      // for(const SpaceCharge* sc1: scsPruned){
      //   for(const SpaceCharge* sc2: scsPruned){
      //     if(sc2 < sc1) continue; // Only visit each pair once

      //     // TODO - we're visiting every pair twice, hence the half
      //     if(sc1->fWire1 && sc1->fWire1 == sc2->fWire1){
      //       Q(sc1->fScratch, sc2->fScratch) += 1;
      //       if(sc1 != sc2) Q(sc2->fScratch, sc1->fScratch) += 1;
      //     }
      //     // both can be true
      //     if(sc1->fWire2 && sc1->fWire2 == sc2->fWire2){
      //       Q(sc1->fScratch, sc2->fScratch) += 1;
      //       if(sc1 != sc2) Q(sc2->fScratch, sc1->fScratch) += 1;
      //     }
      //   }
      // }

      // Only visit each pair once
      for(unsigned int i = 0; i+1 < Nscs; ++i){
        const SpaceCharge* sci = scsPruned[i];
        for(unsigned int j = i+1; j < Nscs; ++j){
          const SpaceCharge* scj = scsPruned[j];
          if(sci->fWire1 && sci->fWire1 == scj->fWire1){
            Q(i, j) += 1;
            Q(j, i) += 1;
          }
          if(sci->fWire2 && sci->fWire2 == scj->fWire2){
            Q(i, j) += 1;
            Q(j, i) += 1;
          }
        }
      }

      // Fill in the diagonal
      for(unsigned int i = 0; i < Nscs; ++i){
        const SpaceCharge* sci = scsPruned[i];
        if(sci->fWire1) Q(i, i) += 1;
        if(sci->fWire2) Q(i, i) += 1;
      }

      // if(Nscs == 4){
      //   T.Print();
      //   Q.Print();
      //   abort();
      // }

      //      VectorT mc = TT*D; // minus c vector

      // Interactions
      if(alpha != 0){
        for(const SpaceCharge* sc: scsPruned){
          //        for(unsigned int scIdx = 0; scIdx < Nscs; ++scIdx){
          //          const SpaceCharge* sc = scsPruned[scIdx];

          const int scIdx = sc->fScratch;

          // Self interaction
          Q(scIdx, scIdx) -= 2*alpha;

          // Fixed interaction with all the points that haven't been zeroed out
          mc(scIdx) += alpha*sc->fNeiPotential;

          for(const Neighbour& n: sc->fNeighbours){
            if(n.fSC->fScratch >= 0){
            //            if(neiIdx.find(n.fSC) == neiIdx.end()){
              // Interaction with an alien point (charge effectively held
              // fixed)
              //              mc(scIdx) += alpha*n.fCoupling*n.fSC->fPred;
              //            }
              //            else{
              // Interaction within this cwire
              Q(scIdx, n.fSC->fScratch) -= alpha*n.fCoupling;
              Q(n.fSC->fScratch, scIdx) -= alpha*n.fCoupling;
              //              Q(scIdx, neiIdx[n.fSC]) -= alpha*n.fCoupling;
              //              Q(neiIdx[n.fSC], scIdx) -= alpha*n.fCoupling;
            }
          }
        }
      }

      // Conservation of charge
      //      Q.ResizeTo(Nscs+1, Nscs+1);
      for(unsigned int i = 0; i < Nscs; ++i){
        Q(Nscs, i) = 1;
        Q(i, Nscs) = 1;
      }

      //      VectorT mcd(Nscs+1);
      //      for(unsigned int i = 0; i < Nscs; ++i) mcd(i) = mc(i);
      //      VectorT& mcd = mc;
      mc(Nscs) = cwire->fCharge; // todo factor of 2?

      //  mc.Print();
      //  abort();

      //      TMatrixFSym Qs(Nscs+1, Q.GetMatrixArray());

      //  std::cout << "Decomp!" << std::endl;
      dc.SetMatrix(Q);
      //  std::cout << "A" << std::endl;
      bool ok;
      VectorT x = dc.Solve(mc, ok);
      //    std::cout << "ok? " << ok << std::endl;

      if(!ok) break;

      if(false){//printing || Q.GetNcols() == 5){
        printing = true;
        std::cout << std::endl << std::endl;
        //        for(auto it: iwIdx) std::cout << it.first << " " << it.second << std::endl;
        //        std::cout << "T = ";
        //        T.Print();
        //        std::cout << "D = ";
        //        D.Print();
        std::cout << "Q = ";
        Q.Print();
        std::cout << "mc = ";
        mc.Print();
        std::cout << "x = ";
        x.Print();
      }

      // Apply findings
      bool any = false;

      // for(unsigned int scIdx = 0; scIdx < Nscs; ++scIdx){
      //   SpaceCharge* sc = scsPruned[scIdx];

      //   if(x[scIdx] <= 0){
      //     //          mask.insert(sc);
      //     sc->fScratch = -1;
      //     //          if(sc->fWire1) sc->fWire1->fScratch = -1;
      //     //          if(sc->fWire2) sc->fWire2->fScratch = -1;
      //     any = true;
      //     // The space charge is already set to zero from the prologue
      //   }
      //   else{
      //     sc->AddCharge(x[scIdx]);
      //   }

      for(SpaceCharge* sc: scsPruned){
        if(x[sc->fScratch] <= 0){
          sc->fScratch = -1;
          any = true;
        }
      }

      if(!any){
        for(SpaceCharge* sc: scsPruned) sc->AddCharge(x[sc->fScratch]);
        break;
      }
      //  std::cout << "Masked " << mask.size() << " space charges" << std::endl;
    } // again!

    if(printing) abort();
  }
};

// ---------------------------------------------------------------------------
void Iterate(SpaceCharge* sc, double alpha)
{
  const QuadExpr chisq = Metric(sc, alpha);

  // Find the minimum of a quadratic expression
  double x = -chisq.Linear()/(2*chisq.Quadratic());

  // Don't allow the SpaceCharge to go negative
  const double xmin = -sc->fPred;

  // Clamp to allowed range
  x = std::max(xmin, x);

  const double chisq_new = chisq.Eval(x);

  // Should try here too, because the function might be convex not concave, so
  // d/dx=0 gives the max not the min, and the true min is at one extreme of
  // the range.
  const double chisq_n = chisq.Eval(xmin);

  if(chisq_n < chisq_new)
    sc->AddCharge(xmin);
  else
    sc->AddCharge(x);
}

// ---------------------------------------------------------------------------
void Iterate(const std::vector<CollectionWireHit*>& cwires,
             const std::vector<SpaceCharge*>& orphanSCs,
             double alpha)
{
  //  SolveQuadProg(cwires, alpha);

  for(CollectionWireHit* cwire: cwires){
    // Appears to never be quicker to be sparse
    if(false){//cwire->fCrossings.size() > 500){
      //      IterateQuadProg<TMatrixDSparse, TDecompSparse, TVectorD>::Iterate(cwire, alpha);
    }
    else{
      // Floats have no effect on speed but also don't seem to harm result
      //      IterateQuadProg<TMatrixD, TDecompLU>::Iterate(cwire, alpha);
      // best?  IterateQuadProg<TMatrixF, TDecompLU, TVectorF>::Iterate(cwire, alpha);
      // slow!      IterateQuadProg<TMatrixF, TDecompSVD, TVectorF>::Iterate(cwire, alpha);
      //      IterateQuadProg<TMatrixF, TDecompQRH, TVectorF>::Iterate(cwire, alpha);

      // Symmetric only
      IterateQuadProg<TMatrixF, TDecompLU, TVectorF>::Iterate(cwire, alpha);
      //      IterateQuadProg<TMatrixD, TDecompBK, TVectorD>::Iterate(cwire, alpha);
      //      IterateQuadProg<TMatrixF, TDecompChol, TVectorF>::Iterate(cwire, alpha);
    }
  }

  //  for(CollectionWireHit* cwire: cwires) Iterate(cwire, alpha);
  //  for(SpaceCharge* sc: orphanSCs) Iterate(sc, alpha);
}

// ---------------------------------------------------------------------------
void SolveQuadProg(const std::vector<CollectionWireHit*>& cwires,
                   double alpha)
{
  // Zero out so we can see which part we work with
  for(CollectionWireHit* c: cwires){
    for(SpaceCharge* sc: c->fCrossings){
      sc->fPred = 0;
    }
  }

  std::set<SpaceCharge*> mask;

  for(int again = 0; again < 100; ++again){

  unsigned int Nscs = 0;
  std::map<InductionWireHit*, int> iwIdx;
  std::map<SpaceCharge*, int> neiIdx;

  unsigned int Ncol = 0;
  for(CollectionWireHit* c: cwires){
    if(c->fCrossings.size() < 2) continue;
    ++Ncol;

    for(SpaceCharge* sc: c->fCrossings){
      if(mask.count(sc) > 0) continue;

      neiIdx[sc] = Nscs;

      ++Nscs;

      if(sc->fWire1 && iwIdx.find(sc->fWire1) == iwIdx.end()) iwIdx[sc->fWire1] = iwIdx.size();
      if(sc->fWire2 && iwIdx.find(sc->fWire2) == iwIdx.end()) iwIdx[sc->fWire2] = iwIdx.size();
    }

    if(Ncol == 300) break; // approx max memory we can handle
  }

  std::cout << Ncol << " cwires" << std::endl;
  std::cout << Nscs << " SCs" << std::endl;
  std::cout << iwIdx.size() << " inductions" << std::endl;
  //  Nscs = std::min(Nscs, 3000u); // TODO: memory. Can we initialize the sparse matrix smarter?
  //  Nscs = 10;

  const unsigned int Nind = iwIdx.size();

  TMatrixDSparse T(Nind, Nscs);
  TVectorD D(Nind);

  TVectorD d(Ncol);
  TMatrixDSparse E(Ncol, Nscs);

  unsigned int colIdx = 0;
  for(CollectionWireHit* c: cwires){
    if(c->fCrossings.size() < 2) continue;

    d(colIdx) = cwires[colIdx]->fCharge;
    ++colIdx;
    if(colIdx >= Ncol) break;
  }

  for(auto it: iwIdx) D(it.second) = it.first->fCharge;

  unsigned int scIdx = 0;
  colIdx = 0;
  for(CollectionWireHit* c: cwires){
    if(c->fCrossings.size() < 2) continue;

    for(SpaceCharge* sc: c->fCrossings){
      if(mask.count(sc) > 0) continue;

      E(colIdx, scIdx) = 1;
      if(sc->fWire1) T(iwIdx[sc->fWire1], scIdx) = 1;
      if(sc->fWire2) T(iwIdx[sc->fWire2], scIdx) = 1;
      ++scIdx;
      if(scIdx >= Nscs) break;
    }
    if(scIdx >= Nscs) break;
    ++colIdx;
    if(colIdx >= Ncol) break;
  }

  //  T.Print();
  //  D.Print();

  TMatrixDSparse TT(TMatrixDSparse::kTransposed, T);
  TMatrixDSparse Q = TT*T;

  // interactions
  colIdx = 0;
  for(CollectionWireHit* c: cwires){
    if(c->fCrossings.size() < 2) continue;

    for(SpaceCharge* sc: c->fCrossings){
      if(mask.count(sc) > 0) continue;

      for(const Neighbour& n: sc->fNeighbours){
        const unsigned int ni = neiIdx[n.fSC];
        if(ni >= Nscs) continue;
        Q(colIdx, ni) += alpha*n.fCoupling;
        Q(ni, colIdx) += alpha*n.fCoupling;
      }
    }
    if(scIdx >= Nscs) break;
    ++colIdx;
    if(colIdx >= Ncol) break;
  }

  //  Q.Print();

  TVectorD mc = TT*D; // minus c vector

  T.ResizeTo(0, 0); TT.ResizeTo(0, 0); // clear these out

  Q.ResizeTo(Nscs+Ncol, Nscs+Ncol);
  Q.SetSub(Nscs, 0, E);
  Q.SetSub(0, Nscs, TMatrixDSparse(TMatrixDSparse::kTransposed, E));

  //  Q.Print();

  TVectorD mcd(Nscs+Ncol);
  for(unsigned int i = 0; i < Nscs; ++i) mcd(i) = mc(i);
  for(unsigned int i = Nscs; i < Nscs+Ncol; ++i) mcd(i) = d(i-Nscs);

  //  mcd.Print();
  //  abort();

  std::cout << "Decomp!" << std::endl;
  TDecompSparse dc;
  dc.SetMatrix(Q);
  std::cout << "A" << std::endl;
  bool ok;
  TVectorD x = dc.Solve(mcd, ok);
  std::cout << "ok? " << ok << std::endl;

  //  x.Print();

  scIdx = 0;
  for(CollectionWireHit* c: cwires){
    if(c->fCrossings.size() < 2) continue;

    for(SpaceCharge* sc: c->fCrossings){
      if(mask.count(sc) > 0) continue;

      if(x[scIdx] < 0) mask.insert(sc);
      sc->fPred = std::max(0., x[scIdx]);

      ++scIdx;
      if(scIdx >= Nscs) break;
    }
    if(scIdx >= Nscs) break;
  }

  std::cout << "Masked " << mask.size() << " space charges" << std::endl;
} // again!
}
