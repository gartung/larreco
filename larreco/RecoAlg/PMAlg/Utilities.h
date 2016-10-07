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

#include "larreco/RecoAlg/PMAlg/LegacyGeomDefs.h"

#include <functional>

namespace pma
{
	typedef std::map< size_t, std::vector<double> > dedx_map;

	class Hit3D;
	class TrkCandidate;
	class bSegmentProjLess;
	class bDistCenterLess2D;
	class bDistCenterLess3D;
	struct bTrajectory3DOrderLess;
	struct bTrajectory3DDistLess;
	struct bTrack3DLonger;

	double Dist2(const Point2D_t& v1, const Point2D_t& v2);
	double Dist2(const Point3D_t& v1, const Point3D_t& v2);
	size_t GetHitsCount(const std::vector< pma::Hit3D* >& hits, unsigned int view);
	double GetSummedADC(const std::vector< pma::Hit3D* >& hits, unsigned int view = geo::kUnknown);
	double GetSummedAmpl(const std::vector< pma::Hit3D* >& hits, unsigned int view = geo::kUnknown);

	double GetHitsRadius3D(const std::vector< pma::Hit3D* >& hits, bool exact = false);
	double GetHitsRadius2D(const std::vector< pma::Hit3D* >& hits, bool exact = false);

	double GetSegmentProjVector(const Point2D_t& p, const Point2D_t& p0, const Point2D_t& p1);
	double GetSegmentProjVector(const Point3D_t& p, const Point3D_t& p0, const Point3D_t& p1);
	Point2D_t GetProjectionToSegment(const Point2D_t& p, const Point2D_t& p0, const Point2D_t& p1);
	Point3D_t GetProjectionToSegment(const Point3D_t& p, const Point3D_t& p0, const Point3D_t& p1);

	double SolveLeastSquares3D(const std::vector< std::pair<Point3D_t, Point3D_t> >& lines, Point3D_t& result);

	Point2D_t GetProjectionToPlane(const Point3D_t& p, unsigned int view, unsigned int tpc, unsigned int cryo);
	Vector2D_t GetVectorProjectionToPlane(const Point3D_t& v, unsigned int view, unsigned int tpc, unsigned int cryo);
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
	public std::binary_function<Point3D_t const*, Point3D_t const*, bool>
{
public:
	bSegmentProjLess(const Point3D_t& s0, const Point3D_t& s1);

	bool operator() (Point3D_t const* p1, Point3D_t const* p2) const
	{
		if (p1 && p2)
		{
			double b1 = pma::GetSegmentProjVector(*p1, segStart, segStop);
			double b2 = pma::GetSegmentProjVector(*p1, segStart, segStop);
			return b1 < b2;
		}
		else return false;
	}

private:
	Point3D_t segStart, segStop;
};

class pma::bDistCenterLess2D :
	public std::binary_function<Point2D_t const&, Point2D_t const&, bool>
{
public:
	bDistCenterLess2D(const Point2D_t& c) : center(c) {}

	bool operator() (Point2D_t const& p1, Point2D_t const& p2)
	{
		double b1 = pma::Dist2(p1, center);
		double b2 = pma::Dist2(p2, center);
		return b1 < b2;
	}

private:
	Point2D_t center;
};

class pma::bDistCenterLess3D :
	public std::binary_function<Point3D_t const&, Point3D_t const&, bool>
{
public:
	bDistCenterLess3D(const Point3D_t& c) : center(c) {}

	bool operator() (Point3D_t const& p1, Point3D_t const& p2)
	{
		double b1 = pma::Dist2(p1, center);
		double b2 = pma::Dist2(p2, center);
		return b1 < b2;
	}

private:
	Point3D_t center;
};


inline double pma::Dist2(const Point2D_t& v1, const Point2D_t& v2)
{
	double dx = v1.X() - v2.X(), dy = v1.Y() - v2.Y();
	return dx * dx + dy * dy;
}

inline double pma::Dist2(const Point3D_t& v1, const Point3D_t& v2)
{
	double dx = v1.X() - v2.X(), dy = v1.Y() - v2.Y(), dz = v1.Z() - v2.Z();
	return dx * dx + dy * dy + dz * dz;
}

#endif

