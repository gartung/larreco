#include "larreco/RecoAlg/ImagePatternAlgs/ShowerVertex/VertexCandidate.h"

#include "TVector3.h"

nnet::VertexCandidate::VertexCandidate(){
  fConfirmed = true;
}

nnet::VertexCandidate::VertexCandidate(geo::WireID thisWire, float thisDrift, unsigned int thisView){
  this->AddView(thisWire,thisDrift,thisView);
  fConfirmed = true;
}

nnet::VertexCandidate::VertexCandidate(geo::WireID w1, float d1, unsigned int v1, geo::WireID w2, float d2, unsigned int v2){
  this->AddView(w1,d1,v1);
  this->AddView(w2,d2,v2);
  fConfirmed = true;
}

// Add a potential vertex from one view
void nnet::VertexCandidate::AddView(geo::WireID thisWire, float thisDrift, unsigned int thisView){
  fWireIDs.push_back(thisWire);
  fDriftX.push_back(thisDrift);
  fViews.push_back(thisView);
}

// Copy constructor
nnet::VertexCandidate::VertexCandidate(const VertexCandidate &rhs){
  fConfirmed = rhs.fConfirmed;
  fWireIDs = rhs.fWireIDs;
  fDriftX = rhs.fDriftX;
  fViews  = rhs.fViews;
}

// If there are three elements, we can see how similar the predicted
// position is between two pairs and return the difference between
// them as a score.
float nnet::VertexCandidate::GetScore() const{

  if(this->GetSize() < 2) return 99999.;

  auto const *geo = lar::providerFrom<geo::Geometry>();

  // For 2D make sure the intersection makes sense but then use
  // the drift agreement for the score 
  if(this->GetSize() == 2){
    TVector3 p1(0,0,0);
    bool s1 = geo->WireIDsIntersect(fWireIDs[0],fWireIDs[1],p1);
    if(s1){
      return this->GetDriftScore();
    }
    else return 9999.;
  }

  // Try seeing if these coordinates make sense.
  TVector3 p1(0,0,0);
  TVector3 p2(0,0,0);

  // If the position isn't inside the TPC it returns false
  bool s1 = geo->WireIDsIntersect(fWireIDs[0],fWireIDs[1],p1);
  bool s2 = geo->WireIDsIntersect(fWireIDs[0],fWireIDs[2],p2);

  if(!s1 || !s2){
    return 9999.;
  }
  else{
    return (p1-p2).Mag();
  }
}

TVector3 nnet::VertexCandidate::Get3DPos() const{
  if(GetSize() < 2){
    return TVector3();
  }
  auto const *geo = lar::providerFrom<geo::Geometry>();
  TVector3 p1(0,0,0);
  geo->WireIDsIntersect(fWireIDs[0],fWireIDs[1],p1); 
  // For the drift coordinate, use the average
  float avDrift = 0;
  for(unsigned int i = 0; i < this->GetSize(); ++i){
    avDrift += fDriftX[i];
  }
  avDrift = avDrift / (1.0*this->GetSize());
  p1.SetX(avDrift);
  return p1;
}

bool nnet::VertexCandidate::IsCompatible(geo::WireID newWire, float drift,float limit) const{

  if(this->GetSize() == 0) return true;
  else{

    // Make sure the WireIDs are in the same TPC
    if(fWireIDs[0].asTPCID().TPC != newWire.asTPCID().TPC){
      std::cout << " Trying to add wires from different TPCs!!! Not compatible!" << std::endl;
      return false;
    }

    // Check drift is compatible
    float dDrift = fDriftX[0] - drift;
    return (fabs(dDrift) < limit);
  }

}

float nnet::VertexCandidate::GetDriftScore() const{

  // For 3D, doesn't really make sense
  if(this->GetSize() == 3) return 999.;
  else{
    return abs(fDriftX[0]-fDriftX[1]);
  }
}

unsigned int nnet::VertexCandidate::GetSize() const{
  return fWireIDs.size();
}

bool nnet::VertexCandidate::HasView(unsigned int view) const{
  return !(std::find(fViews.begin(),fViews.end(),view)==fViews.end());
}

geo::WireID nnet::VertexCandidate::GetWire(unsigned int i) const{
  return fWireIDs[i];
}

float nnet::VertexCandidate::GetDrift(unsigned int i) const{
  return fDriftX[i];
}

unsigned int nnet::VertexCandidate::GetView(unsigned int i) const{
  return fViews[i];
}

void nnet::VertexCandidate::SetFailed(){
  fConfirmed = false;
}

bool nnet::VertexCandidate::IsFailed() const{
  return !fConfirmed;
}

bool nnet::VertexCandidate::ContainsVertex(geo::WireID wire, float x, unsigned int view){
  for(unsigned int i = 0; i < this->GetSize(); ++i){
    if(fWireIDs[i].Wire == wire.Wire && fabs(fDriftX[i]-x)<0.0001 && fViews[i] == view){
      return true;
    }
  }
  return false;
}

void nnet::VertexCandidate::Print() const{
  std::cout << "== Matched Vertex in " << this->GetSize() << " views" << std::endl;
  for(unsigned int i = 0; i < this->GetSize(); ++i){ 
    std::cout << "  - view " << fViews[i] << ", " << fWireIDs[i].Wire << ", " << fDriftX[i] << std::endl; 
  }
  std::cout << "  - with score = " << this->GetScore() << std::endl;
  TVector3 vtx3D = this->Get3DPos(); 
  std::cout << "  - and vertex = (" << vtx3D.X() << ", " << vtx3D.Y() << ", " << vtx3D.Z() << ")" << std::endl;
}

