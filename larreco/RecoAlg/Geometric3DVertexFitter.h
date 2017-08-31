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
    void addTrack(const recob::Track* tk, int pid, double dist) { vtxtracks_.push_back(tk); trackpids_.push_back(pid); propdists_.push_back(dist); }
    //
    void addTrackAndUpdateVertex(const Point_t& pos, const SMatrixSym33& cov, double chi2, int ndof, const recob::Track* tk, int pid, double dist) { 
      pos_ = pos;
      cov_ = cov;
      chi2_+=chi2;
      ndof_+=ndof;
      addTrack(tk, pid, dist);
    }
    //
    const Point_t& position() const { return pos_; }
    const SMatrixSym33& covariance() const { return cov_; }
    const std::vector< const recob::Track* >& tracks() { return vtxtracks_; }
    const std::vector< int >& pids() { return trackpids_; }
    const std::vector< double >& distances() { return propdists_; }
    //
    double chi2() const { return chi2_; }
    double ndof() const { return ndof_; }
    //
    bool isValid() const { return valid_; }
  private:
    Point_t pos_;
    SMatrixSym33 cov_;
    double chi2_;
    int ndof_;
    std::vector< const recob::Track* > vtxtracks_;
    std::vector< int > trackpids_;
    std::vector< double > propdists_;
    bool valid_;
  };
  //
  struct TrackWithPid { 
  public:
    TrackWithPid(const recob::Track* t, int p) : track(t), pid(p) {}
    const recob::Track* track; 
    int pid; 
  };
  //
  //
  //
  class Geometric3DVertexFitter {
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
    FittedVertex fitTracks(std::vector<TrackWithPid>& tracks) const;
    FittedVertex fitTwoTracks(TrackWithPid tk1, TrackWithPid tk2) const;
    void addTrackToVertex(FittedVertex& vtx, TrackWithPid tk) const;
  private:
    std::unique_ptr<TrackStatePropagator> prop;
    int debugLevel;
  };
  //
}

#endif
