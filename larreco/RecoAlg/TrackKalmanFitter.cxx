#include "TrackKalmanFitter.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"

#include "lardata/RecoObjects/SurfWireX.h"
#include "lardata/RecoObjects/KHitWireX.h"
#include "lardata/RecoObjects/Surface.h"
#include "lardata/RecoObjects/KGTrack.h"
#include "lardata/RecoObjects/SurfXYZPlane.h"
#include "lardata/RecoObjects/SurfYZPlane.h"
#include "lardata/RecoObjects/PropYZPlane.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

trkf::KTrack trkf::TrackKalmanFitter::convertRecobTrackIntoKTrack(const TVector3& pos, const TVector3&  dir,  const double pval, const int pdgid) {
  //build a KTrack on the surface orthogonal to the dir at the vertex
  const std::shared_ptr<const trkf::SurfXYZPlane> vtxsurf(new trkf::SurfXYZPlane(pos.X(),  pos.Y(),  pos.Z(), dir.X(), dir.Y(), dir.Z()));
  const double vtxsintheta = std::sin(vtxsurf->theta());
  const double vtxcostheta = std::cos(vtxsurf->theta());
  const double vtxsinphi = std::sin(vtxsurf->phi());
  const double vtxcosphi = std::cos(vtxsurf->phi());
  const double vtxpu = dir.X()*vtxcostheta + dir.Y()*vtxsintheta*vtxsinphi - dir.Z()*vtxsintheta*vtxcosphi;
  const double vtxpv = dir.Y()*vtxcosphi + dir.Z()*vtxsinphi;
  const double vtxpw = dir.X()*vtxsintheta - dir.Y()*vtxcostheta*vtxsinphi + dir.Z()*vtxcostheta*vtxcosphi;
  trkf::TrackVector vtxvec(5);
  vtxvec(0) = (pos.X()-vtxsurf->x0())*vtxcostheta + (pos.Y()-vtxsurf->y0())*vtxsintheta*vtxsinphi - (pos.Z()-vtxsurf->z0())*vtxsintheta*vtxcosphi;
  vtxvec(1) = (pos.Y()-vtxsurf->y0())*vtxcosphi + (pos.Z()-vtxsurf->z0())*vtxsinphi;
  vtxvec(2) = vtxpu/vtxpw;
  vtxvec(3) = vtxpv/vtxpw;
  vtxvec(4) = (pval>0 ? 1./pval : 1.);
  return trkf::KTrack(vtxsurf,vtxvec,trkf::Surface::FORWARD, pdgid);
}

bool trkf::TrackKalmanFitter::fitTrack(const recob::Track& track, const std::vector<art::Ptr<recob::Hit> >& hits,
				       const double pval, const int pdgid, const bool flipDirection,
				       recob::Track& outTrack,    art::PtrVector<recob::Hit>& outHits) {

  auto position = track.Vertex();
  auto direction = track.VertexDirection();

  if (flipDirection) {
    position = track.End();
    direction = -track.EndDirection();
  }

  const trkf::KTrack vtxTrack = convertRecobTrackIntoKTrack(position, direction, pval, pdgid);

  if (vtxTrack.isValid()==0) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__;
    return false;
  }

  //now add an error to this track: KETrack
  trkf::TrackError vtxerr;
  vtxerr.resize(5, false);
  vtxerr.clear();
  //if covariance exists use it (scaled by 100 to remove bias from previous fits), otherwise "dummy" initialization
  if (track.NumberCovariance()) {
    auto covariance = flipDirection ? track.VertexCovariance() : track.EndCovariance();
    for(int i=0; i<5; ++i) {
      for(int j=0; j<5; ++j) {
	vtxerr(i,j) = 100.*covariance(i,j);
      }
    }
  } else {
    vtxTrack.getSurface()->getStartingError(vtxerr);
  }
  const trkf::KETrack tre(vtxTrack, vtxerr);

  //setup the KFitTrack we'll use throughout the fit, it's initial status in INVALID
  trkf::KFitTrack trf = trkf::KFitTrack(tre, 0., 0., trkf::KFitTrack::INVALID);

  //std::cout << "INITIAL TRACK" << std::endl;
  //std::cout << trf.Print(std::cout) << std::endl;
  
  //figure out if hit vector is sorted along or opposite to track direction (ideally the track should be aware of it...)
  trkf::KFitTrack trfVtxF = trkf::KFitTrack(trf, 0., 0.);
  trkf::KFitTrack trfVtxB = trkf::KFitTrack(trf, 0., 0.);
  std::shared_ptr<const trkf::SurfWireX> vtxsurfF(new trkf::SurfWireX(hits.front()->WireID()));
  std::shared_ptr<const trkf::SurfWireX> vtxsurfB(new trkf::SurfWireX(hits.back()->WireID()));
  const KHitWireX hitF(hits.front(),vtxsurfF);
  const KHitWireX hitB(hits.back(),vtxsurfB);
  boost::optional<double> distF = prop_->err_prop(trfVtxF,hitF.getMeasSurface(),trkf::Propagator::UNKNOWN,false);
  boost::optional<double> distB = prop_->err_prop(trfVtxB,hitB.getMeasSurface(),trkf::Propagator::UNKNOWN,false);

  unsigned int nplanes = 0;
  std::vector<KHitWireX> hitsv;
  if (distB.get_value_or(DBL_MAX)<distF.get_value_or(DBL_MAX)) {
    for (auto hit = hits.rbegin(); hit<hits.rend(); ++hit) {
      std::shared_ptr<const trkf::SurfWireX> vtxsurf(new trkf::SurfWireX((*hit)->WireID()));
      hitsv.push_back(std::move(KHitWireX(*hit,vtxsurf)));
    }
  } else {
    for (auto hit = hits.begin(); hit<hits.end(); ++hit) {
      std::shared_ptr<const trkf::SurfWireX> vtxsurf(new trkf::SurfWireX((*hit)->WireID()));
      hitsv.push_back(std::move(KHitWireX(*hit,vtxsurf)));
    }
  }
  for (auto khit = hitsv.begin(); khit<hitsv.end(); ++khit) {
    if ((khit->getHit()->WireID().Plane+1)>nplanes) nplanes = khit->getHit()->WireID().Plane+1;
    if (useRMS_) {
      //0.0833333333 is 1/12, see KHitWireX.cxx line 69
      khit->setMeasError(khit->getMeasError()*khit->getHit()->RMS()*khit->getHit()->RMS()/(std::max(khit->getHit()->SigmaPeakTime()*khit->getHit()->SigmaPeakTime(),0.08333333333f)));
    }
    khit->setMeasError(hitErrScaleFact_*khit->getMeasError());
  }

  //setup the track vector we use to store the fit results
  std::vector<trkf::KFitTrack> fwdPrdTracks;
  std::vector<trkf::KHitTrack> fwdUpdTracks;

  if (sortHitsByPlane_) {
    std::vector<std::vector<unsigned int> > hitsInPlanes(nplanes);
    unsigned int ihit = 0;
    for (auto khit : hitsv) {
      hitsInPlanes[khit.getHit()->WireID().Plane].push_back(ihit++);
    }
    std::vector<unsigned int> iterHitsInPlanes;
    for (auto it : hitsInPlanes) iterHitsInPlanes.push_back(0);
    for (unsigned int p = 0; p<hitsv.size(); ++p) {
      int min_plane = -1;
      double min_dist = DBL_MAX;
      for (unsigned int iplane = 0; iplane<iterHitsInPlanes.size(); ++iplane) {
	for (unsigned int& ih = iterHitsInPlanes[iplane]; ih<hitsInPlanes[iplane].size(); ++ih) {
	  auto& khit = hitsv[hitsInPlanes[iplane][iterHitsInPlanes[iplane]]];
	  //propagate to measurement surface
	  auto trftmp = trf;
	  boost::optional<double> pdist = prop_->noise_prop(trftmp,khit.getMeasSurface(),trkf::Propagator::FORWARD,true);
	  if (!pdist) pdist = prop_->noise_prop(trftmp,khit.getMeasSurface(),trkf::Propagator::BACKWARD,true);
	  if (!pdist) {
	    //mf::LogWarning("TrackKalmanFitter") << "WARNING: both forward and backward propagation failed. Skip this hit...";
	    continue;
	  }
	  const double dist = *pdist;
	  if (skipNegProp_ && dist<0.) continue;
	  if (dist<min_dist) {
	    min_plane = iplane;
	    min_dist = dist;
	  }
	  break;
	}
      }
      //now we know what is the closest plane
      if (min_plane<0) continue;
      auto& khit = hitsv[hitsInPlanes[min_plane][iterHitsInPlanes[min_plane]]];
      iterHitsInPlanes[min_plane]++;
      //propagate to measurement surface
      boost::optional<double> pdist = prop_->noise_prop(trf,khit.getMeasSurface(),trkf::Propagator::FORWARD,true);
      if (!pdist) pdist = prop_->noise_prop(trf,khit.getMeasSurface(),trkf::Propagator::BACKWARD,true);
      if (!pdist) {
	//mf::LogWarning("TrackKalmanFitter") << "WARNING: both forward and backward propagation failed. Skip this hit...";
	continue;
      }
      bool okpred = khit.predict(trf, prop_);
      if (khit.getPredSurface()!=khit.getMeasSurface()) {
	//mf::LogWarning("TrackKalmanFitter") << "WARNING: khit.getPredSurface()!=khit.getMeasSurface(). Skip this hit...";
	continue;
      }
      //now update the forward fitted track
      if (okpred) {
	trf.setPath(trf.getPath()+(*pdist));
	trf.setChisq(trf.getChisq()+khit.getChisq());
	trf.setStat(trkf::KFitTrack::FORWARD_PREDICTED);
	fwdPrdTracks.push_back(trf);
	khit.update(trf);
	trf.setStat(trkf::KFitTrack::FORWARD);
	//std::cout << "UPDATED FWD TRACK" << std::endl;
	//std::cout << trf.Print(std::cout) << std::endl;
	//store this track for the backward fit+smooth
	const std::shared_ptr< const KHitBase > strp(new trkf::KHitWireX(khit));
	trkf::KHitTrack khitTrack(trf, strp);
	fwdUpdTracks.push_back(khitTrack);
      } else {
	mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " " << trf.getStat();
	return false;
      }
      if (trf.isValid()==0) {
	mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " " << trf.getStat();
	return false;
      }
    }
  } else {
    for (auto khit : hitsv) {
      if (skipNegProp_) {
	auto trftmp = trf;//is there a way to avoid creating this copy (i.e. to check propagation distance without modifying the track)?
	boost::optional<double> pdisttest = prop_->noise_prop(trftmp,khit.getMeasSurface(),trkf::Propagator::FORWARD,true);
	if (pdisttest.get_value_or(-1.)<0.) {
	  //mf::LogWarning("TrackKalmanFitter") << "WARNING: negative propagation distance. Skip this hit...";
	  continue;
	}
      }
      //propagate to measurement surface
      boost::optional<double> pdist = prop_->noise_prop(trf,khit.getMeasSurface(),trkf::Propagator::FORWARD,true);
      if (!pdist) pdist = prop_->noise_prop(trf,khit.getMeasSurface(),trkf::Propagator::BACKWARD,true);
      if (!pdist) {
	//mf::LogWarning("TrackKalmanFitter") << "WARNING: both forward and backward propagation failed. Skip this hit...";
	continue;
      }
      bool okpred = khit.predict(trf, prop_);
      if (khit.getPredSurface()!=khit.getMeasSurface()) {
	//mf::LogWarning("TrackKalmanFitter") << "WARNING: khit.getPredSurface()!=khit.getMeasSurface(). Skip this hit...";
	continue;
      }
      //now update the forward fitted track
      if (okpred) {
	trf.setPath(trf.getPath()+(*pdist));
	trf.setChisq(trf.getChisq()+khit.getChisq());
	trf.setStat(trkf::KFitTrack::FORWARD_PREDICTED);
	fwdPrdTracks.push_back(trf);
	khit.update(trf);
	trf.setStat(trkf::KFitTrack::FORWARD);
	//store this track for the backward fit+smooth
	const std::shared_ptr< const KHitBase > strp(new trkf::KHitWireX(khit));
	trkf::KHitTrack khitTrack(trf, strp);
	fwdUpdTracks.push_back(khitTrack);
      } else {
	mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " " << trf.getStat();
	return false;
      }
      if (trf.isValid()==0) {
	mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " " << trf.getStat();
	return false;
      }
    } //for (auto khit : hitsv)
  }

  //std::cout << "TRACK AFTER FWD" << std::endl;
  //std::cout << trf.Print(std::cout) << std::endl;

  //reinitialize trf for backward fit, scale the error to avoid biasing the backward fit
  trf.setError(100.*trf.getError());
  trf.setStat(trkf::KFitTrack::BACKWARD_PREDICTED);
  trf.setChisq(0.);
  
  //backward loop over track states and hits in fwdUpdTracks: use hits for backward fit and fwd track states for smoothing
  float totChi2 = 0.;
  auto fwdPrdTrackIt = fwdPrdTracks.rbegin();
  for (auto fwdUpdTrackIt = fwdUpdTracks.rbegin(); fwdUpdTrackIt != fwdUpdTracks.rend(); ++fwdUpdTrackIt, ++fwdPrdTrackIt) {
    auto& fwdPrdTrack = *fwdPrdTrackIt;
    auto& fwdUpdTrack = *fwdUpdTrackIt;
    trkf::KHitWireX khit(dynamic_cast<const trkf::KHitWireX&>(*fwdUpdTrack.getHit().get()));//need a non const copy in case we want to modify the error
    boost::optional<double> pdist = prop_->noise_prop(trf,khit.getMeasSurface(),trkf::Propagator::BACKWARD,true);
    if (!pdist) pdist = prop_->noise_prop(trf,khit.getMeasSurface(),trkf::Propagator::FORWARD,true);
    if (!pdist) {
      //mf::LogWarning("TrackKalmanFitter") << "WARNING: both forward and backward propagation failed. Skip this hit...";
      continue;
    }
    bool okpred = khit.predict(trf, prop_);
    if (okpred) {
      trf.setPath(trf.getPath()+(*pdist));
      trf.setChisq(trf.getChisq()+khit.getChisq());
      trf.setStat(trkf::KFitTrack::BACKWARD_PREDICTED);
      //combine forward updated and backward, add this to the output track
      fwdPrdTrack.combineFit(trf);
      fwdUpdTrack.combineFit(trf);
      fwdUpdTrack.setPath(trf.getPath()+(*pdist));
      //compute the chi2 between the combined predicted and the hit
      totChi2+=((khit.getMeasVector()[0]-fwdPrdTrack.getVector()[0])*(khit.getMeasVector()[0]-fwdPrdTrack.getVector()[0])/(fwdPrdTrack.getError()(0,0)+khit.getMeasError()(0,0)));
      //now update the backward fitted track
      khit.update(trf);
      trf.setStat(trkf::KFitTrack::BACKWARD);
      //std::cout << "UPDATED BWD TRACK" << std::endl;
      //std::cout << trf.Print(std::cout) << std::endl;
      //std::cout << "COMBINED TRACK" << std::endl;
      //std::cout << fwdUpdTrack.Print(std::cout) << std::endl;
    } else {
      mf::LogWarning("TrackKalmanFitter") << "WARNING invalid predicted surface";
    }
    if (trf.isValid()==0) {
      mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " " << trf.getStat();
      return false;
    }
  }//for (auto fwdUpdTrackIt = fwdUpdTracks.rbegin(); fwdUpdTrackIt != fwdUpdTracks.rend(); ++fwdUpdTrackIt) {

  //std::cout << "TRACK AFTER BWD" << std::endl;
  //std::cout << trf.Print(std::cout) << std::endl;
  
  if (fwdUpdTracks.size()<2) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " ";
    return false;
  }

  //prepare output track
  trkf::KGTrack fittedTrack(0);
  if (sortOutputHitsMinLength_) {
    //try to sort fixing wires order on planes and picking the closest next plane
    std::vector<std::vector<unsigned int> > tracksInPlanes(nplanes);
    unsigned int itrk = 0;
    for (auto fwdUpdTrack : fwdUpdTracks) {
      trkf::KHitWireX khit(dynamic_cast<const trkf::KHitWireX&>(*fwdUpdTrack.getHit().get()));
      tracksInPlanes[khit.getHit()->WireID().Plane].push_back(itrk++);
    }
    //this assumes that the first hit/state is a good one, may want to check if that's the case
    std::vector<unsigned int> iterTracksInPlanes;
    for (auto it : tracksInPlanes) iterTracksInPlanes.push_back(0);
    assert(fwdUpdTracks.front().isValid());
    double pos[3], dir[3];
    fwdUpdTracks.front().getPosition(pos);
    fwdUpdTracks.front().getMomentum(dir);
    for (unsigned int p = 0; p<fwdUpdTracks.size(); ++p) {
      int min_plane = -1;
      double min_dotp = DBL_MAX;
      double tmppos[3], tmpdir[3];
      for (unsigned int iplane = 0; iplane<iterTracksInPlanes.size(); ++iplane) {
	for (unsigned int& itk = iterTracksInPlanes[iplane]; itk<tracksInPlanes[iplane].size(); ++itk) {
	  auto& track = fwdUpdTracks[tracksInPlanes[iplane][iterTracksInPlanes[iplane]]];
	  assert(track.isValid());
	  track.getPosition(tmppos);
	  track.getMomentum(tmpdir);
	  const double dotp = (tmppos[0]-pos[0])*dir[0]+(tmppos[1]-pos[1])*dir[1]+(tmppos[2]-pos[2])*dir[2];
	  if (dotp<min_dotp) {
	    min_plane = iplane;
	    min_dotp = dotp;
	  }
	  break;
	}
      }
      if (min_plane<0) continue;
      auto& track = fwdUpdTracks[tracksInPlanes[min_plane][iterTracksInPlanes[min_plane]]];
      track.setPath(1000.*p);//this preserves the order
      fittedTrack.addTrack(track);
      track.getPosition(pos);
      track.getMomentum(dir);
      iterTracksInPlanes[min_plane]++;
    }
  } else {
    int p = 0;
    for (auto fwdUpdTrackIt = fwdUpdTracks.begin(); fwdUpdTrackIt != fwdUpdTracks.end(); ++fwdUpdTrackIt, ++p) {
      auto& fwdUpdTrack = *fwdUpdTrackIt;
      fwdUpdTrack.setPath(1000.*p);//this preserves the order
      fittedTrack.addTrack(fwdUpdTrack);
    }
  }

  if (fittedTrack.getTrackMap().size()<4) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " ";
    return false;
  }

  bool zeromom = false;
  for (auto itertrack = fittedTrack.getTrackMap().rbegin(); itertrack != fittedTrack.getTrackMap().rend(); ++itertrack) {
    trkf::KHitTrack& trh = itertrack->second;
    double mom[3];
    trh.getMomentum(mom);
    double p = std::sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
    if (p == 0.) zeromom = true;
  }
  if (zeromom) {
    mf::LogWarning("TrackKalmanFitter") << "Fit failure at " << __FILE__ << " " << __LINE__ << " " << trf.getStat();
    return false;
  }

  //fill return objects with smoothed track and its hits
  recob::Track outTrackTmp;
  fittedTrack.fillTrack(outTrackTmp,track.ID());
  std::vector<unsigned int> hittpindex;
  fittedTrack.fillHits(outHits, hittpindex);

  auto covs = outTrackTmp.Covariances();
  int ndof = fittedTrack.getTrackMap().size()-4;//hits are 1D measurement, i.e. each hit is one d.o.f.; no B field: 4 fitted parameters
  outTrack = recob::Track(std::move(outTrackTmp.Trajectory()),pdgid,totChi2,ndof,std::move(covs.first),std::move(covs.second),track.ID());

  return true;

}
