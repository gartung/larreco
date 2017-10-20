#include <vector>
#include "larcore/Geometry/Geometry.h"

namespace nnet{
  class HitCNNOutput;
}

// For each hit (wire w, time t with RMS dt) we want to store the CNN output for 9 points:
// ----------------------
// |  0   |  1   |  2   |
// | w-1  |  w   | w+1  |
// | t+dt | t+dt | t+dt |
// ----------------------
// |  3   |  4   |  5   |
// | w-1  |  w   | w+1  |
// |  t   |  t   |  t   |
// ----------------------
// |  6   |  7   |  8   |
// | w-1  |  w   | w+1  |
// | t-dt | t-dt | t-dt |
// ----------------------

class nnet::HitCNNOutput{

  public:
  
  HitCNNOutput(geo::WireID w, float t, float rms);
  HitCNNOutput(const HitCNNOutput &rhs);

  void AddOutput(unsigned short bin, float CNNValue);
  void AddOutput(unsigned int w, float t, float CNNValue);

  unsigned short GetBin(unsigned int w, float t) const;

  void SetCNNValue(unsigned short bin, float val); 

  float GetCNNValue(unsigned short bin) const;
  float GetCNNValue(unsigned int w, float t) const;

  float GetMaxCNNValue() const;

  /// Allow the user to retrieve the correct wire and time giving the index 0-9
  /// This takes care of the w +/- 1 and t +/- dt
  unsigned int GetWireFromIndex(unsigned short bin) const;
  float GetTimeFromIndex(unsigned short bin) const;

  bool IsComplete() const;

  geo::WireID GetWireID() const;
  unsigned int GetWire() const;
  float GetTime() const;
  float GetTimeRMS() const;

  private:

  std::vector<float> fCNNValues;
  geo::WireID fWire;
  float fTime;
  float fTimeRMS;

};


