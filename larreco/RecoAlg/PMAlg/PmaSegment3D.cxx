/**
 *  @file   PmaSegment3D.cxx
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          3D track segment. See PmaTrack3D.h file for details.
 */

#include "larreco/RecoAlg/PMAlg/PmaSegment3D.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

pma::Segment3D::Segment3D(pma::Track3D* trk, pma::Node3D* vstart, pma::Node3D* vstop) :
	SortedObjectBase(vstart, vstop),
	fParent(trk)
{
	if (vstart->TPC() == vstop->TPC()) fTPC = vstart->TPC();
	if (vstart->Cryo() == vstop->Cryo()) fCryo = vstart->Cryo();
}

double pma::Segment3D::GetDistance2To(const Point3D_t& p3d) const
{
	pma::Node3D* v0 = static_cast< pma::Node3D* >(prev);
	pma::Node3D* v1 = static_cast< pma::Node3D* >(next);
	return GetDist2(p3d, v0->Point3D(), v1->Point3D());
}

double pma::Segment3D::GetDistance2To(const Point2D_t& p2d, unsigned int view) const
{
	pma::Node3D* v0 = static_cast< pma::Node3D* >(prev);
	pma::Node3D* v1 = static_cast< pma::Node3D* >(next);
	return GetDist2(p2d, v0->Projection2D(view), v1->Projection2D(view));
}

pma::Vector3D_t pma::Segment3D::GetDirection3D(void) const
{
	pma::Node3D const* v0 = static_cast< pma::Node3D const* >(prev);
	pma::Node3D const* v1 = static_cast< pma::Node3D const* >(next);
	return (v1->Point3D() - v0->Point3D()).Unit();
}

pma::Point3D_t pma::Segment3D::GetProjection(const Point2D_t& p, unsigned int view) const
{
	pma::Node3D const* vStart = static_cast< pma::Node3D const* >(prev);
	pma::Node3D const* vStop = static_cast< pma::Node3D const* >(next);

	Vector2D_t v0 = (p - vStart->Projection2D(view));

	Vector2D_t v1 = (vStop->Projection2D(view) - vStart->Projection2D(view));
	
	Vector3D_t v3d = (vStop->Point3D() - vStart->Point3D());
	
	Point3D_t v3dStart = vStart->Point3D();
	Point3D_t v3dStop = vStop->Point3D();

	double v0Norm = v0.R();
	double v1Norm = v1.R();

	Point3D_t result;
	double eps = 1.0E-6; // 0.01mm
	if (v1Norm > eps)
	{
		double mag = v0Norm * v1Norm;
		double cosine = 0.0;
		if (mag != 0.0) cosine = v0.Dot(v1) / mag;
		double b = v0Norm * cosine / v1Norm;

		if (b < 1.0)
		{
			result = v3dStart;
			if (b > 0.0) result += (v3d * b);
		}
		else result = v3dStop;
	}
	else // segment 2D projection is almost a point
	{
		mf::LogWarning("pma::Segment3D") << "Short segment projection.";

		result = middlePoint(v3dStart, v3dStop);
	}
	return result;
}

pma::Point3D_t pma::Segment3D::GetUnconstrainedProj3D(const Point2D_t& p2d, unsigned int view) const
{
	pma::Node3D const* vStart = static_cast< pma::Node3D const* >(prev);
	pma::Node3D const* vStop = static_cast< pma::Node3D const* >(next);

	Vector2D_t const v0 = (p2d - vStart->Projection2D(view));

	Vector2D_t const v1 = (vStop->Projection2D(view) - vStart->Projection2D(view));

	Vector3D_t const v3d = (vStop->Point3D() - vStart->Point3D());

	double v0Norm = v0.R();
	double v1Norm = v1.R();

	double eps = 1.0E-6; // 0.01mm
	if (v1Norm > eps)
	{
		double mag = v0Norm * v1Norm;
		double cosine = 0.0;
		if (mag != 0.0) cosine = v0.Dot(v1) / mag;
		double b = v0Norm * cosine / v1Norm;

		return vStart->Point3D() + (v3d * b);
	}
	else // segment 2D projection is almost a point
	{
		mf::LogWarning("pma::Segment3D") << "Short segment projection." << std::endl;

		return middlePoint(vStart->Point3D(), vStop->Point3D());
	}
}

void pma::Segment3D::SetProjection(pma::Hit3D& h) const
{
	pma::Node3D* vStart = static_cast< pma::Node3D* >(prev);
	pma::Node3D* vStop = static_cast< pma::Node3D* >(next);

	auto const & pointStart = vStart->Point3D();
	auto const & pointStop = vStop->Point3D();

	auto const & projStart = vStart->Projection2D(h.View2D());
	auto const & projStop = vStop->Projection2D(h.View2D());

	pma::Vector2D_t v0 = (h.Point2D() - projStart);

	pma::Vector2D_t v1 = (projStop - projStart);

	pma::Vector3D_t v3d = (pointStop - pointStart);

	double v0Norm = v0.R();
	double v1Norm = v1.R();

	double eps = 1.0E-6; // 0.01mm
	if (v1Norm > eps)
	{
		double mag = v0Norm * v1Norm;
		double cosine = 0.0;
		if (mag != 0.0) cosine = v0.Dot(v1) / mag;
		double b = v0Norm * cosine / v1Norm;

		pma::Vector2D_t p(projStart.X(), projStart.Y());
		p += (v1 * b);
		v3d *= b;

		h.SetProjection(p.X(), p.Y(), (float)b);
		h.SetPoint3D(
			vStart->Point3D().X() + v3d.X(),
			vStart->Point3D().Y() + v3d.Y(),
			vStart->Point3D().Z() + v3d.Z());
	}
	else // segment 2D projection is almost a point
	{
		mf::LogWarning("pma::Segment3D") << "Short segment projection.";

		h.SetProjection(middlePoint(projStart, projStop), 0.0F);

		h.SetPoint3D(middlePoint(pointStart, pointStop));
	}
}

double pma::Segment3D::Length2(void) const
{
	if (prev && next)
		return pma::Dist2( ((pma::Node3D*)prev)->Point3D(), ((pma::Node3D*)next)->Point3D() );
	else
	{
		mf::LogError("pma::Segment3D") << "Segment endpoints not set.";
		return 0.0;
	}
}


double pma::Segment3D::GetDist2(const Point3D_t& psrc, const Point3D_t& p0, const Point3D_t& p1)
{
	pma::Vector3D_t v0(psrc.X() - p0.X(), psrc.Y() - p0.Y(), psrc.Z() - p0.Z());
	pma::Vector3D_t v1(p1.X() - p0.X(), p1.Y() - p0.Y(), p1.Z() - p0.Z());

	pma::Vector3D_t v2(psrc.X() - p1.X(), psrc.Y() - p1.Y(), psrc.Z() - p1.Z());
	pma::Vector3D_t v3(v1); v3 *= -1.0;

	double v0Norm2 = v0.Mag2();
	double v1Norm2 = v1.Mag2();

	double eps = 1.0E-6; // 0.01mm
	if (v1Norm2 < eps)
	{
		mf::LogWarning("pma::Segment3D") << "Short segment or its projection.";

		double dx = 0.5 * (p0.X() + p1.X()) - psrc.X();
		double dy = 0.5 * (p0.Y() + p1.Y()) - psrc.Y();
		double dz = 0.5 * (p0.Z() + p1.Z()) - psrc.Z();
		return dx * dx + dy * dy + dz * dz;
	}

	double v0v1 = v0.Dot(v1);
	double v2v3 = v2.Dot(v3);
	double v2Norm2 = v2.Mag2();

	double result = 0.0;
	if ((v0v1 > 0.0) && (v2v3 > 0.0))
	{
		double cosine01_square = 0.0;
		double mag01_square = v0Norm2 * v1Norm2;
		if (mag01_square != 0.0) cosine01_square = v0v1 * v0v1 / mag01_square;

		result = (1.0 - cosine01_square) * v0Norm2;
	}
	else // increase distance to prefer hit assigned to the vertex, not segment
	{
		if (v0v1 <= 0.0) result = 1.0001 * v0Norm2;
		else result = 1.0001 * v2Norm2;
	}

	if (result >= 0.0) return result;
	else return 0.0;
}

double pma::Segment3D::GetDist2(const Point2D_t& psrc, const Point2D_t& p0, const Point2D_t& p1)
{
	pma::Vector2D_t v0(psrc.X() - p0.X(), psrc.Y() - p0.Y());
	pma::Vector2D_t v1(p1.X() - p0.X(), p1.Y() - p0.Y());

	pma::Vector2D_t v2(psrc.X() - p1.X(), psrc.Y() - p1.Y());
	pma::Vector2D_t v3(v1); v3 *= -1.0;

	double v0Norm2 = v0.Mag2();
	double v1Norm2 = v1.Mag2();

	double eps = 1.0E-6; // 0.01mm
	if (v1Norm2 < eps)
	{
		mf::LogVerbatim("pma::Segment3D") << "Short segment or its projection.";

		double dx = 0.5 * (p0.X() + p1.X()) - psrc.X();
		double dy = 0.5 * (p0.Y() + p1.Y()) - psrc.Y();
		return dx * dx + dy * dy;
	}

	double v0v1 = v0.Dot(v1);
	double v2v3 = v2.Dot(v3);
	double v2Norm2 = v2.Mag2();

	double result = 0.0;
	if ((v0v1 > 0.0) && (v2v3 > 0.0))
	{
		double cosine01_square = 0.0;
		double mag01_square = v0Norm2 * v1Norm2;
		if (mag01_square != 0.0) cosine01_square = v0v1 * v0v1 / mag01_square;

		result = (1.0 - cosine01_square) * v0Norm2;
	}
	else // increase distance to prefer hit assigned to the vertex, not segment
	{
		if (v0v1 <= 0.0) result = 1.0001 * v0Norm2;
		else result = 1.0001 * v2Norm2;
	}

	if (result >= 0.0) return result;
	else return 0.0;
}


double pma::Segment3D::GetDist2(const Vector2D_t& psrc, const Vector2D_t& p0, const Vector2D_t& p1)
{
	pma::Vector2D_t v0(psrc.X() - p0.X(), psrc.Y() - p0.Y());
	pma::Vector2D_t v1(p1.X() - p0.X(), p1.Y() - p0.Y());

	pma::Vector2D_t v2(psrc.X() - p1.X(), psrc.Y() - p1.Y());
	pma::Vector2D_t v3(v1); v3 *= -1.0;

	double v0Norm2 = v0.Mag2();
	double v1Norm2 = v1.Mag2();

	double eps = 1.0E-6; // 0.01mm
	if (v1Norm2 < eps)
	{
		mf::LogVerbatim("pma::Segment3D") << "Short segment or its projection.";

		double dx = 0.5 * (p0.X() + p1.X()) - psrc.X();
		double dy = 0.5 * (p0.Y() + p1.Y()) - psrc.Y();
		return dx * dx + dy * dy;
	}

	double v0v1 = v0.Dot(v1);
	double v2v3 = v2.Dot(v3);
	double v2Norm2 = v2.Mag2();

	double result = 0.0;
	if ((v0v1 > 0.0) && (v2v3 > 0.0))
	{
		double cosine01_square = 0.0;
		double mag01_square = v0Norm2 * v1Norm2;
		if (mag01_square != 0.0) cosine01_square = v0v1 * v0v1 / mag01_square;

		result = (1.0 - cosine01_square) * v0Norm2;
	}
	else // increase distance to prefer hit assigned to the vertex, not segment
	{
		if (v0v1 <= 0.0) result = 1.0001 * v0Norm2;
		else result = 1.0001 * v2Norm2;
	}

	if (result >= 0.0) return result;
	else return 0.0;
}
