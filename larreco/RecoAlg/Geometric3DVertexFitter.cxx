#include "larreco/RecoAlg/Geometric3DVertexFitter.h"

trkf::FittedVertex trkf::Geometric3DVertexFitter::fitPFP(size_t iPF, const art::ValidHandle<std::vector<recob::PFParticle> >& inputPFParticle, 
							 const std::unique_ptr<art::FindManyP<recob::Track> >& assocTracks) const
{
  using namespace std;
  //
  art::Ptr<recob::PFParticle> pfp(inputPFParticle, iPF);
  //
  if (debugLevel>1) std::cout << "pfp#" << iPF << " PdgCode=" << pfp->PdgCode() 
			      << " IsPrimary=" << pfp->IsPrimary()
			      << " NumDaughters=" << pfp->NumDaughters()
			      << std::endl;
  if (pfp->IsPrimary()==false || pfp->NumDaughters()<2) FittedVertex();
  
  vector< art::Ptr<recob::Track> > tracks;
  
  auto& pfd = pfp->Daughters();
  for (auto ipfd : pfd) {
    vector< art::Ptr<recob::Track> > pftracks = assocTracks->at(ipfd);
    art::Ptr<recob::PFParticle> pfpd(inputPFParticle, ipfd);
    for (auto t : pftracks) {
      tracks.push_back(t);
    }
  }
  if (tracks.size()<2) return FittedVertex();
  //
  return fitTracks(tracks);
}

trkf::FittedVertex trkf::Geometric3DVertexFitter::fitTracks(std::vector<art::Ptr<recob::Track> >& tracks) const
{
  if (debugLevel>0) std::cout << "fitting vertex with ntracks=" << tracks.size() << std::endl;
  if (tracks.size()==0) return FittedVertex();
  if (tracks.size()==1) {
    FittedVertex vtx;
    vtx.addTrackAndUpdateVertex(tracks[0]->Start(), tracks[0]->VertexCovarianceGlobal6D().Sub<SMatrixSym33>(0,0), 0, 0, tracks[0], 0);
    return vtx;
  }
  // sort tracks by number of hits
  std::sort(tracks.begin(), tracks.end(), [](art::Ptr<recob::Track> a, art::Ptr<recob::Track> b) {return a->CountValidPoints() > b->CountValidPoints();});
  //
  // find vertex between the two tracks with most hits
  FittedVertex vtx = fitTwoTracks(tracks[0], tracks[1]);
  if (vtx.isValid()==false) return vtx;
  //
  // then add other tracks and update vertex measurement
  for (auto tk = tracks.begin()+2; tk<tracks.end(); ++tk) {
    addTrackToVertex(vtx, *tk);
  }
  return vtx;
}

trkf::FittedVertex trkf::Geometric3DVertexFitter::fitTwoTracks(const art::Ptr<recob::Track> tk1, const art::Ptr<recob::Track> tk2) const
{
  // find the closest approach points
  auto start1 = tk1->Trajectory().Start();
  auto start2 = tk2->Trajectory().Start();

  auto dir1 = tk1->Trajectory().StartDirection();
  auto dir2 = tk2->Trajectory().StartDirection();

  if (debugLevel>0) {
    std::cout << "track1 with start=" << start1 << " dir=" << dir1
	      << " length=" << tk1->Length() << " points=" << tk1->CountValidPoints()
	      << std::endl;
    std::cout << "covariance=\n" << tk1->VertexCovarianceGlobal6D() << std::endl;
    std::cout << "track2 with start=" << start2 << " dir=" << dir2
              << " length=" << tk2->Length() << " points=" << tk2->CountValidPoints()
  	      << std::endl;
    std::cout << "covariance=\n" << tk2->VertexCovarianceGlobal6D() << std::endl;
  }

  auto dpos = start1-start2;
  auto dotd1d2 = dir1.Dot(dir2);

  auto dotdpd1 = dpos.Dot(dir1);
  auto dotdpd2 = dpos.Dot(dir2);

  auto dist2 = ( dotd1d2*dotdpd1 - dotdpd2 )/( dotd1d2*dotd1d2 - 1 );
  auto dist1 = ( dotd1d2*dist2 - dotdpd1 );

  //by construction both point of closest approach on the two lines lie on this plane
  recob::tracking::Plane target(start1+dist1*dir1, dir1);

  recob::tracking::Plane plane1(start1, dir1);
  trkf::TrackState state1(tk1->VertexParametersLocal5D(), tk1->VertexCovarianceLocal5D(), plane1, true, tk1->ParticleId());
  bool propok1 = true;
  state1 = prop->propagateToPlane(propok1, state1, target, true, true, trkf::TrackStatePropagator::UNKNOWN);
  if (!propok1) return FittedVertex();

  recob::tracking::Plane plane2(start2, dir2);
  trkf::TrackState state2(tk2->VertexParametersLocal5D(), tk2->VertexCovarianceLocal5D(), plane2, true, tk2->ParticleId());
  bool propok2 = true;
  state2 = prop->propagateToPlane(propok2, state2, target, true, true, trkf::TrackStatePropagator::UNKNOWN);
  if (!propok2) return FittedVertex();

  if (debugLevel>0) {
    std::cout << "track1 pca=" << start1+dist1*dir1 << " dist=" << dist1 << std::endl;
    std::cout << "track2 pca=" << start2+dist2*dir2 << " dist=" << dist2 << std::endl;
  }

  if (debugLevel>1) {
    std::cout << "target pos=" << target.position() << " dir=" << target.direction() << std::endl;
    //test if orthogonal
    auto dcp = state1.position()-state2.position();
    std::cout << "dot dcp-dir1=" << dcp.Dot(tk1->Trajectory().StartDirection()) << std::endl;
    std::cout << "dot dcp-dir2=" << dcp.Dot(tk2->Trajectory().StartDirection()) << std::endl;
  }

  //here we neglect that to find dist1 and dist2 one track depends on the other...
  SMatrixSym22 cov1 = state1.covariance().Sub<SMatrixSym22>(0,0);
  SMatrixSym22 cov2 = state2.covariance().Sub<SMatrixSym22>(0,0);

  SVector2 par1(state1.parameters()[0],state1.parameters()[1]);
  SVector2 par2(state2.parameters()[0],state2.parameters()[1]);
  SVector2 deltapar = par2 - par1;

  SMatrixSym22 covsum = (cov2+cov1);
  //
  if (debugLevel>1) {
    std::cout << "par1=" << par1 << std::endl;
    std::cout << "par2=" << par2 << std::endl;
    std::cout << "deltapar=" << deltapar << std::endl;
    //
    std::cout << "cov1=\n" << cov1 << std::endl;
    std::cout << "cov2=\n" << cov2 << std::endl;
    std::cout << "covsum=\n" << covsum << std::endl;
  }
  //
  if (debugLevel>1) {
    double det1;
    bool d1ok = cov1.Det2(det1);
    std::cout << "cov1 det=" << det1 << " ok=" << d1ok << std::endl;
    double det2;
    bool d2ok = cov2.Det2(det2);
    std::cout << "cov2 det=" << det2 << " ok=" << d2ok << std::endl;
    double detsum;
    bool dsok = covsum.Det2(detsum);
    std::cout << "covsum det=" << detsum << " ok=" << dsok << std::endl;
  }
  //
  bool invertok = covsum.Invert();
  if (!invertok) return FittedVertex();
  auto k = cov1*covsum;

  if (debugLevel>1) {
    std::cout << "inverted covsum=\n" << covsum << std::endl;
    std::cout << "k=\n" << k << std::endl;
    std::cout << "dist1=" << dist1 << " dist2=" << dist2 << " propok1=" << propok1 << " propok2=" << propok2 << " invertok=" << invertok << std::endl;
  }

  SVector2 vtxpar2 = par1 + k*deltapar;
  SMatrixSym22 vtxcov22 = cov1 - ROOT::Math::SimilarityT(cov1,covsum);

  if (debugLevel>1) {
    std::cout << "vtxpar2=" << vtxpar2 << std::endl;
    std::cout << "vtxcov22=\n" << vtxcov22 << std::endl;
  }
  
  auto chi2 = ROOT::Math::Similarity(deltapar,covsum);

  SVector5 vtxpar5(vtxpar2[0],vtxpar2[1],0,0,0);
  SMatrixSym55 vtxcov55;vtxcov55.Place_at(vtxcov22,0,0);
  TrackState vtxstate(vtxpar5, vtxcov55, target, true, 0);

  const int ndof = 1; //1=2*2-3; each measurement is 2D because it is defined on a plane
  trkf::FittedVertex vtx(vtxstate.position(), vtxstate.covariance6D().Sub<SMatrixSym33>(0,0), chi2, ndof);
  vtx.addTrack(tk1, dist1);
  vtx.addTrack(tk2, dist2);

  if (debugLevel>0) {
    std::cout << "vtxpos=" << vtx.position() << std::endl;
    std::cout << "vtxcov=\n" << vtx.covariance() << std::endl;
    std::cout << "chi2=" << chi2 << std::endl;
  }

  return vtx;
}

trkf::Geometric3DVertexFitter::ParsCovsPlaneDist trkf::Geometric3DVertexFitter::getParsCovsPlaneDist(const trkf::FittedVertex& vtx, const art::Ptr<recob::Track> tk) const {
  auto start = tk->Trajectory().Start();
  auto dir = tk->Trajectory().StartDirection();

  if (debugLevel>0) {
    std::cout << "adding track with start=" << start << " dir=" << dir
	      << " length=" << tk->Length() << " points=" << tk->CountValidPoints()
	      << std::endl;
    std::cout << "covariance=\n" << tk->VertexCovarianceGlobal6D() << std::endl;
  }

  auto vtxpos = vtx.position();
  auto vtxcov = vtx.covariance();

  auto dist = dir.Dot(vtxpos-start);

  //by construction vtxpos also lies on this plane
  recob::tracking::Plane target(start+dist*dir, dir);

  recob::tracking::Plane plane(start, dir);
  trkf::TrackState state(tk->VertexParametersLocal5D(), tk->VertexCovarianceLocal5D(), plane, true, tk->ParticleId());
  bool propok = true;
  state = prop->propagateToPlane(propok, state, target, true, true, trkf::TrackStatePropagator::UNKNOWN);

  if (debugLevel>0) {
    std::cout << "input vtx=" << vtxpos << std::endl;
    std::cout << "track pca=" << start+dist*dir << " dist=" << dist << std::endl;
  }

  //rotate the vertex position and covariance to the target plane
  SVector6 vtxpar6(vtxpos.X(),vtxpos.Y(),vtxpos.Z(),0,0,0);
  SMatrixSym66 vtxcov66;vtxcov66.Place_at(vtxcov,0,0);

  auto vtxpar5 = target.Global6DToLocal5DParameters(vtxpar6);
  SVector2 par1(vtxpar5[0],vtxpar5[1]);
  SVector2 par2(state.parameters()[0],state.parameters()[1]);

  //here we neglect that to find dist, the vertex is used...
  SMatrixSym22 cov1 = target.Global6DToLocal5DCovariance(vtxcov66,false,Vector_t()).Sub<SMatrixSym22>(0,0);
  SMatrixSym22 cov2 = state.covariance().Sub<SMatrixSym22>(0,0);

  return ParsCovsPlaneDist(par1,par2,cov1,cov2,target,dist);
}

void trkf::Geometric3DVertexFitter::addTrackToVertex(trkf::FittedVertex& vtx, const art::Ptr<recob::Track> tk) const
{

  ParsCovsPlaneDist pcp = getParsCovsPlaneDist(vtx, tk);
  SVector2& par1 = pcp.par1;
  SVector2& par2 = pcp.par2;
  SMatrixSym22& cov1 = pcp.cov1;
  SMatrixSym22& cov2 = pcp.cov2;
  recob::tracking::Plane& target = pcp.plane;
  double dist = pcp.dist;

  SVector2 deltapar = par2 - par1;

  SMatrixSym22 covsum = (cov2+cov1);
  //
  if (debugLevel>1) {
    std::cout << "par1=" << par1 << std::endl;
    std::cout << "par2=" << par2 << std::endl;
    std::cout << "deltapar=" << deltapar << std::endl;
    //
    std::cout << "cov1=\n" << cov1 << std::endl;
    std::cout << "cov2=\n" << cov2 << std::endl;
    std::cout << "covsum=\n" << covsum << std::endl;
  }
  //
  if (debugLevel>1) {
    double det1;
    bool d1ok = cov1.Det2(det1);
    std::cout << "cov1 det=" << det1 << " ok=" << d1ok << std::endl;
    double det2;
    bool d2ok = cov2.Det2(det2);
    std::cout << "cov2 det=" << det2 << " ok=" << d2ok << std::endl;
    double detsum;
    bool dsok = covsum.Det2(detsum);
    std::cout << "covsum det=" << detsum << " ok=" << dsok << std::endl;
  }
  //
  bool invertok = covsum.Invert();
  if (!invertok) return;
  auto chi2 = ROOT::Math::Similarity(deltapar,covsum);

  auto k = cov1*covsum;

  if (debugLevel>1) {
    std::cout << "inverted covsum=\n" << covsum << std::endl;
    std::cout << "k=\n" << k << std::endl;
    // std::cout << "dist=" << dist << " propok=" << propok << " invertok=" << invertok << std::endl;
  }

  SVector2 updvtxpar2 = par1 + k*deltapar;
  SMatrixSym22 updvtxcov22 = cov1 - ROOT::Math::SimilarityT(cov1,covsum);

  if (debugLevel>1) {
    std::cout << "updvtxpar2=" << updvtxpar2 << std::endl;
    std::cout << "updvtxcov22=\n" << updvtxcov22 << std::endl;
  }
  
  SVector5 updvtxpar5(updvtxpar2[0],updvtxpar2[1],0,0,0);
  SMatrixSym55 updvtxcov55;updvtxcov55.Place_at(updvtxcov22,0,0);
  TrackState updvtxstate(updvtxpar5, updvtxcov55, target, true, 0);

  const int ndof = 2;// Each measurement is 2D because it is defined on a plane
  vtx.addTrackAndUpdateVertex(updvtxstate.position(), updvtxstate.covariance6D().Sub<SMatrixSym33>(0,0), chi2, ndof, tk, dist);

  if (debugLevel>0) {
    std::cout << "updvtxpos=" << vtx.position() << std::endl;
    std::cout << "updvtxcov=\n" << vtx.covariance() << std::endl;
    std::cout << "add chi2=" << chi2 << std::endl;
  }

  return;
}

double trkf::Geometric3DVertexFitter::chi2(const trkf::Geometric3DVertexFitter::ParsCovsPlaneDist& pcp) const
{
  const SVector2 deltapar = pcp.par2 - pcp.par1;
  SMatrixSym22 covsum = (pcp.cov2+pcp.cov1);
  //
  bool invertok = covsum.Invert();
  if (!invertok) return -1.;
  //
  return ROOT::Math::Similarity(deltapar,covsum);
}

double trkf::Geometric3DVertexFitter::chi2(const FittedVertex& vtx, const art::Ptr<recob::Track> tk) const 
{
  return chi2(getParsCovsPlaneDist(vtx, tk));
}

double trkf::Geometric3DVertexFitter::ip(const trkf::Geometric3DVertexFitter::ParsCovsPlaneDist& pcp) const
{
  const SVector2 deltapar = pcp.par2 - pcp.par1;
  return std::sqrt( deltapar[0]*deltapar[0] + deltapar[1]*deltapar[1] );
}

double trkf::Geometric3DVertexFitter::ip(const FittedVertex& vtx, const art::Ptr<recob::Track> tk) const
{
  return ip(getParsCovsPlaneDist(vtx, tk));
}

double trkf::Geometric3DVertexFitter::ipErr(const trkf::Geometric3DVertexFitter::ParsCovsPlaneDist& pcp) const
{
  SVector2 deltapar = pcp.par2 - pcp.par1;
  deltapar/=std::sqrt( deltapar[0]*deltapar[0] + deltapar[1]*deltapar[1] );
  SMatrixSym22 covsum = (pcp.cov2+pcp.cov1);
  return ROOT::Math::Similarity(deltapar,covsum);
}

double trkf::Geometric3DVertexFitter::ipErr(const FittedVertex& vtx, const art::Ptr<recob::Track> tk) const
{
  return ipErr(getParsCovsPlaneDist(vtx, tk));
}

double trkf::Geometric3DVertexFitter::sip(const trkf::Geometric3DVertexFitter::ParsCovsPlaneDist& pcp) const
{
  return ip(pcp)/ipErr(pcp);
}

double trkf::Geometric3DVertexFitter::sip(const FittedVertex& vtx, const art::Ptr<recob::Track> tk) const
{
  return sip(getParsCovsPlaneDist(vtx, tk));
}

trkf::FittedVertex trkf::Geometric3DVertexFitter::unbiasedVertex(const trkf::FittedVertex& vtx, const art::Ptr<recob::Track> tk) const 
{
  auto ittoerase = vtx.findTrack(tk);
  if (ittoerase == vtx.tracks().size()) {
    return vtx;
  } else {
    std::vector< art::Ptr<recob::Track> > tks = vtx.tracksCopy();
    tks.erase(tks.begin()+ittoerase);
    return fitTracks(tks);
  }
}

double trkf::Geometric3DVertexFitter::chi2Unbiased(const trkf::FittedVertex& vtx, const art::Ptr<recob::Track> tk) const {
  auto ittoerase = vtx.findTrack(tk);
  if (ittoerase == vtx.tracks().size()) {
    return chi2(vtx,tk);
  } else {
    std::vector< art::Ptr<recob::Track> > tks = vtx.tracksCopy();
    tks.erase(tks.begin()+ittoerase);
    return chi2(fitTracks(tks),tk);
  }
}

double trkf::Geometric3DVertexFitter::ipUnbiased(const trkf::FittedVertex& vtx, const art::Ptr<recob::Track> tk) const {
  auto ittoerase = vtx.findTrack(tk);
  if (ittoerase == vtx.tracks().size()) {
    return ip(vtx,tk);
  } else {
    std::vector< art::Ptr<recob::Track> > tks = vtx.tracksCopy();
    tks.erase(tks.begin()+ittoerase);
    return ip(fitTracks(tks),tk);
  }
}

double trkf::Geometric3DVertexFitter::ipErrUnbiased(const trkf::FittedVertex& vtx, const art::Ptr<recob::Track> tk) const {
  auto ittoerase = vtx.findTrack(tk);
  if (ittoerase == vtx.tracks().size()) {
    return ipErr(vtx,tk);
  } else {
    std::vector< art::Ptr<recob::Track> > tks = vtx.tracksCopy();
    tks.erase(tks.begin()+ittoerase);
    return ipErr(fitTracks(tks),tk);
  }
}

double trkf::Geometric3DVertexFitter::sipUnbiased(const trkf::FittedVertex& vtx, const art::Ptr<recob::Track> tk) const {
  auto ittoerase = vtx.findTrack(tk);
  if (ittoerase == vtx.tracks().size()) {
    return sip(vtx,tk);
  } else {
    std::vector< art::Ptr<recob::Track> > tks = vtx.tracksCopy();
    tks.erase(tks.begin()+ittoerase);
    return sip(fitTracks(tks),tk);
  }
}
