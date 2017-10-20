#include <iostream>
#include <cmath>

#include "larreco/RecoAlg/ImagePatternAlgs/ShowerVertex/HitCNNOutput.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

nnet::HitCNNOutput::HitCNNOutput(geo::WireID w, float t, float rms){
  fWire = w;
  fTime = t;
  fTimeRMS = rms;

  fCNNValues = std::vector<float>(9,-999.);
}

nnet::HitCNNOutput::HitCNNOutput(const HitCNNOutput &rhs){
  fWire = rhs.fWire;
  fTime = rhs.fTime;
  fTimeRMS = rhs.fTimeRMS;

  fCNNValues = rhs.fCNNValues;
}

void nnet::HitCNNOutput::AddOutput(unsigned short bin, float CNNValue){
  if(bin < 9){
    fCNNValues[bin] = CNNValue;
  }
  else{
    mf::LogError("HitCNNOutput::AddOutput") << "bin out of range." << std::endl;
  }
}

void nnet::HitCNNOutput::AddOutput(unsigned int w, float t, float CNNValue){
  return this->AddOutput(this->GetBin(w,t),CNNValue);
}

unsigned short nnet::HitCNNOutput::GetBin(unsigned int w,float t) const{

  unsigned short bin = 999;
  float eps = 1e-6;

  if(w == fWire.Wire - 1){
    if(fabs(t-(fTime+fTimeRMS)) < eps) bin = 0;
    if(fabs(t-(fTime)) < eps) bin = 3;
    if(fabs(t-(fTime-fTimeRMS)) < eps) bin = 6; 
  }
  else if(w == fWire.Wire){
    if(fabs(t-(fTime+fTimeRMS)) < eps) bin = 1;
    if(fabs(t-(fTime)) < eps) bin = 4;
    if(fabs(t-(fTime-fTimeRMS)) < eps) bin = 7;
  }
  else if(w == fWire.Wire+1){
    if(fabs(t-(fTime+fTimeRMS)) < eps) bin = 2;
    if(fabs(t-(fTime)) < eps) bin = 5;
    if(fabs(t-(fTime-fTimeRMS)) < eps) bin = 8;
  }
 
  if(bin == 999){
    mf::LogError("HitCNNOutput::GetBin") << " Can't find bin for " << w << ", " << t << std::endl;
  }
 
  return bin;
}

void nnet::HitCNNOutput::SetCNNValue(unsigned short bin, float val){
  if(bin < 9){
    fCNNValues[bin] = val;
  }
  else{
    mf::LogError("HitCNNOutput::SetCNNValue") << "bin out of range" << std::endl;
  }
}

float nnet::HitCNNOutput::GetCNNValue(unsigned short bin) const{
  float val = -999.;
  if(bin < 9){
    val = fCNNValues[bin];
  }
  else{
    mf::LogError("HitCNNOutput::GetCNNValue") << "bin out of range, returning " << val << std::endl;
  }
  return val;
}

float nnet::HitCNNOutput::GetCNNValue(unsigned int w, float t) const{
  return this->GetCNNValue(this->GetBin(w,t));
}

float nnet::HitCNNOutput::GetMaxCNNValue() const{
  float max = -999;
  for(unsigned int v = 0; v < fCNNValues.size(); ++v){
    if(fCNNValues[v] > max) max = fCNNValues[v]; 
  }
  return max;
}

unsigned int nnet::HitCNNOutput::GetWireFromIndex(unsigned short bin) const{
  unsigned int w = 999999;
  if(bin < 9){
    
    switch(bin){
      // Cases for fWire -1
      case 0:
      case 3:
      case 6: w = fWire.Wire-1; break;
      // Cases for fWire
      case 1:
      case 4:
      case 7: w = fWire.Wire; break;
      // Cases for fWire+1
      case 2:
      case 5:
      case 8: w = fWire.Wire+1; break;
      // Default
      default: w = 999999;
    }

  }
  else{
    mf::LogError("HitCNNOutput::GetWireFromIndex") << "bin out of range, returning " << w << std::endl;
  }
  return w;
}

float nnet::HitCNNOutput::GetTimeFromIndex(unsigned short bin) const{
  float t = -9999.;
  if(bin < 9){
    switch(bin){
      // t + dt
      case 0:
      case 1:
      case 2: t = fTime + fTimeRMS; break;
      // t
      case 3:
      case 4:
      case 5: t = fTime; break;
      // t - dt
      case 6:
      case 7:
      case 8: t = fTime - fTimeRMS; break;
      // default
      default: t = -9999.;
    }
  }
  else{
    mf::LogError("HitCNNOutput::GetTimeFromIndex") << "bin out of range, returning " << t << std::endl;
  }
  return t;
}

bool nnet::HitCNNOutput::IsComplete() const{
  unsigned int count = 0;
  for(auto const v : fCNNValues){
    if(v > -998.){
      ++count;
    }
  }
  return (count == fCNNValues.size());
}

geo::WireID nnet::HitCNNOutput::GetWireID() const{
  return fWire;
}

unsigned int nnet::HitCNNOutput::GetWire() const{
  return fWire.Wire;
}

float nnet::HitCNNOutput::GetTime() const{
  return fTime;
}

float nnet::HitCNNOutput::GetTimeRMS() const{
  return fTimeRMS;
}

