////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgDetermineT0
// Author:      L. Whitehead (leigh.howard.whitehead@cern.ch)
//
// Purpose:     This class is designed to collect together different methods to assign T0 to
//              tracks in the TPCs in the PMA reconstruction  
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PMAlgDetermineT0_h
#define PMAlgDetermineT0_h

#include <map>

#include "fhiclcpp/types/Atom.h"
#include "larreco/RecoAlg/PMAlg/PmaTrkCandidate.h"

#include "TVector3.h"

namespace pma{
  class PMAlgDetermineT0;
}

class pma::PMAlgDetermineT0 {

  public:

  struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> DriftWidthMargin {
          Name("DriftWidthMargin"),
          Comment("The minimum distance beyond 1 drift window required for tagging track as a cosmic background.")
      };
  };

  void DetermineT0(pma::TrkCandidateColl& tracks);

  PMAlgDetermineT0(const pma::PMAlgDetermineT0::Config &config);

  private:
  
  // For now, the only case we need to worry about is tracks crossing one drift volume
  // Those that span two volumes have T0 assigned by the stitching algorithms
  size_t GetAPAToCPAT0(pma::TrkCandidateColl& tracks);

  // Get the width in the drift direction for each TPC
  void GetTPCDriftWidths();

  double fDriftWidthMargin; // The leeway we give ourselves on saying a track traversed the entire drift length of the TPC
  std::map<unsigned int, double> fDriftMin; // The minimum value of the drift coordinate for each TPC
  std::map<unsigned int, double> fDriftMax; // The maximum value of the drift coordinate for each TPC

  short int fDriftCoord; // The drift direction coordinate (x = 0, y = 1, z = 2)
  std::map<unsigned int, short int> fDriftDir;   // The drift direction (positive or negative in the above coordinate)
};

#endif

