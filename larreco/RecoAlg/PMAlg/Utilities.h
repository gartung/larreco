/**
 *  @file   Utilities.h
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          Some geometrical functions and sorting helpers.
 *          See PmaTrack3D.h file for details.
 */

#ifndef Utilities_h
#define Utilities_h

#include "larcore/Geometry/Geometry.h"

#include "larreco/RecoAlg/PMAlg/GeomDefs.h"

#include <functional>

#include "TVector2.h"
#include "TVector3.h"

#include "Math/GenVector/DisplacementVector2D.h"
#include "Math/GenVector/DisplacementVector3D.h"

namespace pma
{
	typedef ROOT::Math::DisplacementVector2D< ROOT::Math::Cartesian2D<double> > Vector2D;
	typedef ROOT::Math::DisplacementVector3D< ROOT::Math::Cartesian3D<double> > Vector3D;

	typedef std::map< size_t, std::vector<double> > dedx_map;

	class Hit3D;
	class TrkCandidate;
	class bSegmentProjLess;
	class bDistCenterLess2D;
	class bDistCenterLess3D;
	struct bTrajectory3DOrderLess;
	struct bTrajectory3DDistLess;
	struct bTrack3DLonger;

	double Dist2(const TVector2& v1, const TVector2& v2);
	double Dist2(const TVector3& v1, const TVector3& v2);
	double Dist2(const Point2D_t& v1, const Point2D_t& v2);
	double Dist2(const Point3D_t& v1, const Point3D_t& v2);
	size_t GetHitsCount(const std::vector< pma::Hit3D* >& hits, unsigned int view);
	double GetSummedADC(const std::vector< pma::Hit3D* >& hits, unsigned int view = geo::kUnknown);
	double GetSummedAmpl(const std::vector< pma::Hit3D* >& hits, unsigned int view = geo::kUnknown);

	double GetHitsRadius3D(const std::vector< pma::Hit3D* >& hits, bool exact = false);
	double GetHitsRadius2D(const std::vector< pma::Hit3D* >& hits, bool exact = false);

	double GetSegmentProjVector(const TVector2& p, const TVector2& p0, const TVector2& p1);
	double GetSegmentProjVector(const TVector3& p, const TVector3& p0, const TVector3& p1);
	Point2D_t GetProjectionToSegment(const TVector2& p, const TVector2& p0, const TVector2& p1);
	Point3D_t GetProjectionToSegment(const TVector3& p, const TVector3& p0, const TVector3& p1);

	double SolveLeastSquares3D(const std::vector< std::pair<TVector3, TVector3> >& lines, TVector3& result);

	Point2D_t GetProjectionToPlane(const TVector3& p, unsigned int view, unsigned int tpc, unsigned int cryo);
	Point2D_t GetVectorProjectionToPlane(const TVector3& v, unsigned int view, unsigned int tpc, unsigned int cryo);
	Point2D_t WireDriftToCm(unsigned int wire, float drift, unsigned int view, unsigned int tpc, unsigned int cryo);
	Point2D_t CmToWireDrift(float xw, float yd, unsigned int view, unsigned int tpc, unsigned int cryo);
}


struct pma::bTrajectory3DOrderLess :
	public std::binary_function<pma::Hit3D*, pma::Hit3D*, bool>
{
	bool operator() (pma::Hit3D* h1, pma::Hit3D* h2);
};

struct pma::bTrajectory3DDistLess :
	public std::binary_function<pma::Hit3D*, pma::Hit3D*, bool>
{
	bool operator() (pma::Hit3D* h1, pma::Hit3D* h2);
};

struct pma::bTrack3DLonger :
	public std::binary_function<const pma::TrkCandidate &, const pma::TrkCandidate &, bool>
{
	bool operator() (const pma::TrkCandidate & t1, const pma::TrkCandidate & t2);
};

class pma::bSegmentProjLess :
	public std::binary_function<TVector3*, TVector3*, bool>
{
public:
	bSegmentProjLess(const TVector3& s0, const TVector3& s1);

	bool operator() (TVector3* p1, TVector3* p2)
	{
		if (p1 && p2)
		{
			double b1 = pma::GetSegmentProjVector(*p1, segStart, segStop);
			double b2 = pma::GetSegmentProjVector(*p1, segStart, segStop);
			return b1 < b2;
		}
		else return false;
	}
	
	bool operator() (Point3D_t const* p1, Point3D_t const* p2)
    {
      TVector3 v1 = makeTVector3(*p1), v2 = makeTVector3(*p2);
      return this->operator()(&v1, &v2);
    }

private:
	TVector3 segStart, segStop;
};

class pma::bDistCenterLess2D :
	public std::binary_function<TVector2, TVector2, bool>
{
public:
	bDistCenterLess2D(const TVector2& c) : center(c) {}

	bool operator() (TVector2 p1, TVector2 p2)
	{
		double b1 = pma::Dist2(p1, center);
		double b2 = pma::Dist2(p2, center);
		return b1 < b2;
	}

private:
	TVector2 center;
};

class pma::bDistCenterLess3D :
	public std::binary_function<TVector3, TVector3, bool>
{
public:
	bDistCenterLess3D(const TVector3& c) : center(c) {}

	bool operator() (TVector3 p1, TVector3 p2)
	{
		double b1 = pma::Dist2(p1, center);
		double b2 = pma::Dist2(p2, center);
		return b1 < b2;
	}

private:
	TVector3 center;
};

#endif

