////////////////////////////////////////////////////////////////////////////
//
// \brief Definition of 3D cluster object for LArSoft
//
// \author usher@slac.stanford.edu
//
////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>

#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"
#include "lardataobj/RecoBase/Hit.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace reco{

ClusterHit2D::ClusterHit2D() : m_statusBits(0),
                               m_docaToAxis(9999.),
                               m_arcLenToPoca(0.),
                               m_xPosition(0.),
                               m_timeTicks(0.),
                               m_hit(nullptr) {}

ClusterHit2D::ClusterHit2D(unsigned           statusBits,
                           float              doca,
                           float              poca,
                           float              xPosition,
                           float              timeTicks,
                           const recob::Hit&  hit) :
                           m_statusBits(statusBits),
                           m_docaToAxis(doca),
                           m_arcLenToPoca(poca),
                           m_xPosition(xPosition),
                           m_timeTicks(timeTicks),
                           m_hit(&hit) {}

    
std::ostream& operator<< (std::ostream& o, const ClusterHit2D& c)
{
    o << c.getHit();
    
    return o;
}
    
bool operator < (const ClusterHit2D & a, const ClusterHit2D & b)
{
    return a.getHit() < b.getHit();
}


ClusterHit3D::ClusterHit3D() : m_id(std::numeric_limits<size_t>::max()),
                               m_statusBits(0),
                               m_position{0.,0.,0.},
                               m_totalCharge(0.),
                               m_avePeakTime(-1.),
                               m_deltaPeakTime(0.),
                               m_sigmaPeakTime(0.),
                               m_hitChiSquare(0.),
                               m_docaToAxis(0.),
                               m_arclenToPoca(0.)
{
    m_hitDelTSigVec.clear();
    m_wireIDVector.clear();
    m_hitVector.clear();
    
    m_hitDelTSigVec.resize(3, 0.);
    m_wireIDVector.resize(3, geo::WireID());
    m_hitVector.resize(3, NULL);
}
    
ClusterHit3D::ClusterHit3D(size_t                          id,
                           unsigned int                    statusBits,
                           const float*                    position,
                           float                           totalCharge,
                           float                           avePeakTime,
                           float                           deltaPeakTime,
                           float                           sigmaPeakTime,
                           float                           hitChiSquare,
                           float                           docaToAxis,
                           float                           arclenToPoca,
                           const ClusterHit2DVec&          hitVec,
                           const std::vector<float>&       hitDelTSigVec,
                           const std::vector<geo::WireID>& wireIDs) :
              m_id(id),
              m_statusBits(statusBits),
              m_position{position[0],position[1],position[2]},
              m_totalCharge(totalCharge),
              m_avePeakTime(avePeakTime),
              m_deltaPeakTime(deltaPeakTime),
              m_sigmaPeakTime(sigmaPeakTime),
              m_hitChiSquare(hitChiSquare),
              m_docaToAxis(docaToAxis),
              m_arclenToPoca(arclenToPoca),
              m_hitDelTSigVec(hitDelTSigVec),
              m_wireIDVector(wireIDs)
{
    m_hitVector.resize(3,NULL);
    std::copy(hitVec.begin(),hitVec.end(),m_hitVector.begin());
}
    
ClusterHit3D::ClusterHit3D(const ClusterHit3D& toCopy)
{
    m_id             = toCopy.m_id;
    m_statusBits     = toCopy.m_statusBits;
    m_position[0]    = toCopy.m_position[0];
    m_position[1]    = toCopy.m_position[1];
    m_position[2]    = toCopy.m_position[2];
    m_totalCharge    = toCopy.m_totalCharge;
    m_avePeakTime    = toCopy.m_avePeakTime;
    m_deltaPeakTime  = toCopy.m_deltaPeakTime;
    m_sigmaPeakTime  = toCopy.m_sigmaPeakTime;
    m_hitChiSquare   = toCopy.m_hitChiSquare;
    m_docaToAxis     = toCopy.m_docaToAxis;
    m_arclenToPoca   = toCopy.m_arclenToPoca;
    m_hitVector      = toCopy.m_hitVector;
    m_hitDelTSigVec  = toCopy.m_hitDelTSigVec;
    m_wireIDVector   = toCopy.m_wireIDVector;
}
    
void ClusterHit3D::initialize(size_t                          id,
                              unsigned int                    statusBits,
                              const float*                    position,
                              float                           totalCharge,
                              float                           avePeakTime,
                              float                           deltaPeakTime,
                              float                           sigmaPeakTime,
                              float                           hitChiSquare,
                              float                           docaToAxis,
                              float                           arclenToPoca,
                              const ClusterHit2DVec&          hitVec,
                              const std::vector<float>&       hitDelTSigVec,
                              const std::vector<geo::WireID>& wireIDs)
{
    m_id             = id;
    m_statusBits     = statusBits;
    m_position[0]    = position[0];
    m_position[1]    = position[1];
    m_position[2]    = position[2];
    m_totalCharge    = totalCharge;
    m_avePeakTime    = avePeakTime;
    m_deltaPeakTime  = deltaPeakTime;
    m_sigmaPeakTime  = sigmaPeakTime;
    m_hitChiSquare   = hitChiSquare;
    m_docaToAxis     = docaToAxis;
    m_arclenToPoca   = arclenToPoca;
    m_hitVector      = hitVec;
    m_hitDelTSigVec  = hitDelTSigVec;
    m_wireIDVector   = wireIDs;
    
    return;
}

void ClusterHit3D::setWireID(const geo::WireID& wid) const
{
    m_wireIDVector[wid.Plane] = wid;
}
    
std::ostream& operator<< (std::ostream& o, const ClusterHit3D& c)
{
    o << "ClusterHit3D has " << c.getHits().size() << " hits associated";
    
    return o;
}
    
  //bool operator < (const ClusterHit3D & a, const ClusterHit3D & b)
  //{
  //    if (a.m_position[2] != b.m_position[2]) return a.m_position[2] < b.m_position[2];
  //    else return a.m_position[0] < b.m_position[0];
  //}

PrincipalComponents::PrincipalComponents() :
       m_svdOK(false),
       m_numHitsUsed(0),
       m_eigenValues{0.,0.,0.},
       m_avePosition{0.,0.,0.},
       m_aveHitDoca(9999.)
{}    
    
PrincipalComponents::PrincipalComponents(bool ok, int nHits, const float* eigenValues, const EigenVectors& eigenVecs, const float* avePos, const float aveHitDoca) :
           m_svdOK(ok),
           m_numHitsUsed(nHits),
           m_eigenVectors(eigenVecs),
           m_aveHitDoca(aveHitDoca)
{
    m_eigenValues[0]     = eigenValues[0];
    m_eigenValues[1]     = eigenValues[1];
    m_eigenValues[2]     = eigenValues[2];
    m_avePosition[0]     = avePos[0];
    m_avePosition[1]     = avePos[1];
    m_avePosition[2]     = avePos[2];
}
    
void PrincipalComponents::flipAxis(size_t axisDir)
{
    std::vector<float>& axis = m_eigenVectors.at(axisDir);
    
    for(auto& val : axis) val *= -1.;
    
    return;
}
    
std::ostream&  operator << (std::ostream & o, const PrincipalComponents& a)
{
    if (a.m_svdOK)
    {
        o << std::setiosflags(std::ios::fixed) << std::setprecision(2);
        o << " PCAxis ID run with " << a.m_numHitsUsed << " space points" << std::endl;
        o << "   - center position: " << std::setw(6) << a.m_avePosition[0] << ", " << a.m_avePosition[1] << ", " << a.m_avePosition[2] << std::endl;
        o << "   - eigen values: " << std::setw(8) << std::right << a.m_eigenValues[0] << ", "
        << a.m_eigenValues[1] << ", " << a.m_eigenValues[2] << std::endl;
        o << "   - average doca: " << a.m_aveHitDoca << std::endl;
        o << "   - Principle axis: " << std::setw(7) << std::setprecision(4) << a.m_eigenVectors[0][0] << ", " << a.m_eigenVectors[0][1] << ", " << a.m_eigenVectors[0][2] << std::endl;
        o << "   - second axis:    " << std::setw(7) << std::setprecision(4) << a.m_eigenVectors[1][0] << ", " << a.m_eigenVectors[1][1] << ", " << a.m_eigenVectors[1][2] << std::endl;
        o << "   - third axis:     " << std::setw(7) << std::setprecision(4) << a.m_eigenVectors[2][0] << ", " << a.m_eigenVectors[2][1] << ", " << a.m_eigenVectors[2][2] << std::endl;
    }
    else
        o << " Principal Components Axis is not valid" << std::endl;
    
    return o;
}
    
bool operator < (const PrincipalComponents& a, const PrincipalComponents& b)
{
    if (a.m_svdOK && b.m_svdOK)
        return a.m_eigenValues[0] > b.m_eigenValues[0];
    
    return false; //They are equal
}

Cluster3D::Cluster3D() : m_statusBits(0),
                         m_pcaResults(PrincipalComponents()),
                         m_totalCharge(0.),
                         m_startPosition{0.,0.,0.},
		                 m_endPosition{0.,0.,0.},
		                 m_clusterIdx(0)
{}

Cluster3D::Cluster3D(unsigned                   statusBits,
                     const PrincipalComponents& pcaResults,
                     float                     totalCharge,
                     const float*              startPosition,
                     const float*              endPosition,
                     int                        idx) :
        m_statusBits(statusBits),
        m_pcaResults(pcaResults),
        m_totalCharge(totalCharge),
        m_startPosition{startPosition[0],startPosition[1],startPosition[2]},
        m_endPosition{endPosition[0],endPosition[1],endPosition[2]},
        m_clusterIdx(idx)
       {}
    
//----------------------------------------------------------------------
//  Addition operator.
//
Cluster3D Cluster3D::operator +(Cluster3D a)
{
/*
    // throw exception if the clusters are not from the same plane
    if( a.View() != this->View() )
      throw cet::exception("Cluster+operator") << "Attempting to sum clusters from "
                 << "different views is not allowed\n";

    // check the start and end positions - for now the
    // smallest wire number means start position, largest means end position
    std::vector<float> astart(a.StartPos());
    std::vector<float> aend  (a.EndPos()  );
    std::vector<float> start(StartPos());
    std::vector<float> end  (EndPos()  );
    std::vector<float> sigstart(SigmaStartPos());
    std::vector<float> sigend  (SigmaEndPos()  );

    if(astart[0] < fStartPos[0]){
      start = astart;
      sigstart = a.SigmaStartPos();
    }

    if(aend[0] > fEndPos[0]){
      end = aend;
      sigend = a.SigmaEndPos();
    }

    //take weighted mean in obtaining average slope and differential charge,
    //based on total charge each cluster
    float dtdw = ((this->Charge()*dTdW()) + (a.Charge()*a.dTdW()))/(this->Charge() + a.Charge());
    float dqdw = ((this->Charge()*dQdW()) + (a.Charge()*a.dQdW()))/(this->Charge() + a.Charge());

    //hits.sort();//sort the PtrVector to organize Hits of new Cluster
    float sigdtdw = TMath::Max(SigmadTdW(), a.SigmadTdW());
    float sigdqdw = TMath::Max(SigmadQdW(), a.SigmadQdW());

    Cluster sum(//hits,
    start[0], sigstart[0],
    start[1], sigstart[1],
    end[0],   sigend[0],
    end[1],   sigend[1],
    dtdw, sigdtdw,
    dqdw, sigdqdw,
    this->Charge() + a.Charge(),
    this->View(),
    ID());
*/
    //return sum;
    return a;
}

//----------------------------------------------------------------------
// ostream operator.
//
std::ostream& operator<< (std::ostream& o, const Cluster3D& c)
{
    o << std::setiosflags(std::ios::fixed) << std::setprecision(2);
    o << "Cluster ID "    << std::setw(5)  << std::right << c.getClusterIdx();
//      << " : View = "     << std::setw(3)  << std::right << c.View()
//      << " StartWire = "  << std::setw(7)  << std::right << c.StartPos()[0]
//      << " EndWire = "    << std::setw(7)  << std::right << c.EndPos()[0]
//      << " StartTime = "  << std::setw(9)  << std::right << c.StartPos()[1]
//      << " EndTime = "    << std::setw(9)  << std::right << c.EndPos()[1]
//      << " dTdW = "       << std::setw(9)  << std::right << c.dTdW()
//      << " dQdW = "       << std::setw(9)  << std::right << c.dQdW()
//      << " Charge = "     << std::setw(10) << std::right << c.Charge();

    return o;
  }


//----------------------------------------------------------------------
// < operator.
//
bool operator < (const Cluster3D & a, const Cluster3D & b)
{
/*
    if(a.View() != b.View())
      return a.View() < b.View();
    if(a.ID() != b. ID())
      return a.ID() < b.ID();
    if(a.StartPos()[0] != b.StartPos()[0])
      return a.StartPos()[0] < b.StartPos()[0];
    if(a.EndPos()[0] != b.EndPos()[0])
      return a.EndPos()[0] < b.EndPos()[0];
*/
    if (a.getStartPosition()[2] < b.getStartPosition()[2]) return true;

    return false; //They are equal
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

void RecobClusterParameters::UpdateParameters(const reco::ClusterHit2D* clusterHit)
{
    /**
     *  @brief a utility routine for building 3D clusters to keep basic info up to date
     *         (a candidate for a better way to do this)
     */
    const recob::Hit& hit = clusterHit->getHit();
    
    // Need to keep track of stuff so we can form cluster
    if (hit.WireID().Wire < m_startWire)
    {
        m_startWire      = hit.WireID().Wire;
        m_startTime      = hit.PeakTimeMinusRMS();
        m_sigmaStartTime = hit.SigmaPeakTime();
    }
    
    if (hit.WireID().Wire > m_endWire)
    {
        m_endWire      = hit.WireID().Wire;
        m_endTime      = hit.PeakTimePlusRMS();
        m_sigmaEndTime = hit.SigmaPeakTime();
    }
    
    m_totalCharge += hit.Integral();
    m_plane        = hit.WireID().Plane;
    m_view         = hit.View();
    
    m_hitVector.push_back(clusterHit);
    
    return;
}

}// namespace

