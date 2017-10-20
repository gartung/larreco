#include "TVector3.h"
#include "larcore/Geometry/Geometry.h"

namespace nnet{
  class VertexCandidate;
}

// Convenience class to store vertex candidates used in the shower verte finding code
class nnet::VertexCandidate {

  public:

    VertexCandidate();
    VertexCandidate(geo::WireID thisWire, float thisDrift, unsigned int thisView);

    VertexCandidate(geo::WireID w1, float d1, unsigned int v1, geo::WireID w2, float d2, unsigned int v2);

    // Add a potential vertex from one view
    void AddView(geo::WireID thisWire, float thisDrift, unsigned int thisView);

    // Copy constructor
    VertexCandidate(const VertexCandidate &rhs);

    // If there are three elements, we can see how similar the predicted
    // position is between two pairs and return the difference between
    // them as a score.
    float GetScore() const;

    TVector3 Get3DPos() const;

    bool IsCompatible(geo::WireID newWire, float drift,float limit) const;

    float GetDriftScore() const;

    unsigned int GetSize() const;

    bool HasView(unsigned int view) const;

    geo::WireID GetWire(unsigned int i) const;
    float GetDrift(unsigned int i) const;
    unsigned int GetView(unsigned int i) const;

    void SetFailed();
    bool IsFailed() const;

    bool ContainsVertex(geo::WireID wire, float x, unsigned int view);

    void Print() const;

  private:
    // Store the information to identify the individual 2D vertices
    std::vector<geo::WireID> fWireIDs;
    std::vector<float> fDriftX;
    std::vector<unsigned int> fViews;

    // Default is true, set to false when we are sure this isn't a good match.
    bool fConfirmed; 
};


