#ifndef GEOMETRIC3DVERTEXFITTER_H
#define GEOMETRIC3DVERTEXFITTER_H

////////////////////////////////////////////////////////////////////////
// Class:       Geometric3DVertexFitter
// File:        Geometric3DVertexFitter.h
//
// Author: Giuseppe Cerati, cerati@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardata/RecoObjects/TrackStatePropagator.h"

namespace trkf {
  //
  using SMatrixSym22 = ROOT::Math::SMatrix< double, 2, 2, ROOT::Math::MatRepSym< double, 2 > >;
  using SVector2     = ROOT::Math::SVector<double,2>;
  using SMatrixSym33 = ROOT::Math::SMatrix< double, 3, 3, ROOT::Math::MatRepSym< double, 3 > >;
  using SVector3     = ROOT::Math::SVector<double,3>;
  //
  struct FittedVertex {
  public:
    FittedVertex() : valid_(false) {}
    FittedVertex(const Point_t& pos, const SMatrixSym33& cov, double chi2, int ndof)
    : pos_(pos), cov_(cov), chi2_(chi2), ndof_(ndof), valid_(true) {}
    //
    void addTrack(art::Ptr<recob::Track> tk, double dist) { vtxtracks_.push_back(tk); propdists_.push_back(dist); }
    //
    void addTrackAndUpdateVertex(const Point_t& pos, const SMatrixSym33& cov, double chi2, int ndof, art::Ptr<recob::Track> tk, double dist) { 
      pos_ = pos;
      cov_ = cov;
      chi2_+=chi2;
      ndof_+=ndof;
      addTrack(tk, dist);
    }
    //
    const Point_t& position() const { return pos_; }
    const SMatrixSym33& covariance() const { return cov_; }
    const std::vector< art::Ptr<recob::Track> >& tracks() const { return vtxtracks_; }
    std::vector< art::Ptr<recob::Track> > tracksCopy() const { return vtxtracks_; }
    const std::vector< double >& distances() const { return propdists_; }
    //
    double chi2() const { return chi2_; }
    double ndof() const { return ndof_; }
    //
    bool isValid() const { return valid_; }
    size_t findTrack(art::Ptr<recob::Track> tk) const {
      for (size_t it = 0; it!=vtxtracks_.size(); ++it) {
	if (tk==vtxtracks_[it]) return it;
      }
      return vtxtracks_.size();
    }
  private:
    Point_t pos_;
    SMatrixSym33 cov_;
    double chi2_;
    int ndof_;
    std::vector< art::Ptr<recob::Track> > vtxtracks_;
    std::vector< double > propdists_;
    bool valid_;
  };
  //
  //
  class Geometric3DVertexFitter {
    // TODO:
    // test chi2, ip, ipErr, sip functions
    // find a cut value to accept tracks in vertex fit
    // find a way to have core functionalities rely only on plain recob::Tracks and not art::Ptr so that they can be used outside full art
  public:
    //
    struct Options {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int> debugLevel {
	Name("debugLevel"),
	Comment("Debugging level: 0 for no printouts, 1 for minimal, 2 for full.")
      };
    };

    // Constructor
    Geometric3DVertexFitter(const fhicl::Table<Options>& o, const fhicl::Table<TrackStatePropagator::Config>& p)
      : debugLevel(o().debugLevel())
      {
	prop = std::make_unique<TrackStatePropagator>(p);
      }

    FittedVertex fitPFP(size_t iPF, const art::ValidHandle<std::vector<recob::PFParticle> >& inputPFParticle, 
			const std::unique_ptr<art::FindManyP<recob::Track> >& assocTracks) const;
    FittedVertex fitTracks(std::vector<art::Ptr<recob::Track> >& tracks) const;
    FittedVertex fitTwoTracks(const art::Ptr<recob::Track> tk1, const art::Ptr<recob::Track> tk2) const;
    void addTrackToVertex(FittedVertex& vtx, const art::Ptr<recob::Track> tk) const;
    double chi2(const FittedVertex& vtx, const art::Ptr<recob::Track> tk) const;
    double ip(const FittedVertex& vtx, const art::Ptr<recob::Track> tk) const;
    double ipErr(const FittedVertex& vtx, const art::Ptr<recob::Track> tk) const;
    double sip(const FittedVertex& vtx, const art::Ptr<recob::Track> tk) const;
    FittedVertex unbiasedVertex(const FittedVertex& vtx, const art::Ptr<recob::Track> tk) const;
    double chi2Unbiased(const FittedVertex& vtx, const art::Ptr<recob::Track> tk) const;
    double ipUnbiased(const FittedVertex& vtx, const art::Ptr<recob::Track> tk) const;
    double ipErrUnbiased(const FittedVertex& vtx, const art::Ptr<recob::Track> tk) const;
    double sipUnbiased(const FittedVertex& vtx, const art::Ptr<recob::Track> tk) const;
  private:
    std::unique_ptr<TrackStatePropagator> prop;
    int debugLevel;
    //
    struct ParsCovsPlaneDist {
    ParsCovsPlaneDist(SVector2 p1, SVector2 p2, SMatrixSym22 c1, SMatrixSym22 c2, recob::tracking::Plane p, double d)
    : par1(p1), par2(p2), cov1(c1), cov2(c2), plane(p), dist(d) {}
      SVector2& par1, par2;
      SMatrixSym22& cov1, cov2;
      recob::tracking::Plane& plane;
      double dist;
    };
    ParsCovsPlaneDist getParsCovsPlaneDist(const trkf::FittedVertex& vtx, const art::Ptr<recob::Track> tk) const;
    double chi2(const ParsCovsPlaneDist& pcp) const;
    double ip(const ParsCovsPlaneDist& pcp) const;
    double ipErr(const ParsCovsPlaneDist& pcp) const;
    double sip(const ParsCovsPlaneDist& pcp) const;
  };
  //
}

#endif
