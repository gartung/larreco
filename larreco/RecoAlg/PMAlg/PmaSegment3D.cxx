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

double pma::Segment3D::GetDistance2To(const TVector3& p3d) const
{
	pma::Node3D* v0 = static_cast< pma::Node3D* >(prev);
	pma::Node3D* v1 = static_cast< pma::Node3D* >(next);
	return GetDist2(makePoint3D(p3d), v0->Point3D(), v1->Point3D());
}

double pma::Segment3D::GetDistance2To(const TVector2& p2d, unsigned int view) const
{
	pma::Node3D* v0 = static_cast< pma::Node3D* >(prev);
	pma::Node3D* v1 = static_cast< pma::Node3D* >(next);
	return GetDist2(p2d, makeTVector2(v0->Projection2D(view)), makeTVector2(v1->Projection2D(view)));
}

pma::Vector3D_t pma::Segment3D::GetDirection3D(void) const
{
	pma::Node3D* v0 = static_cast< pma::Node3D* >(prev);
	pma::Node3D* v1 = static_cast< pma::Node3D* >(next);
	TVector3 dir = makeTVector3(v1->Point3D() - v0->Point3D());
	dir *= 1.0 / dir.Mag();
	return makeVector3D(dir);
}

pma::Point3D_t pma::Segment3D::GetProjection(const TVector2& p, unsigned int view) const
{
	pma::Node3D* vStart = static_cast< pma::Node3D* >(prev);
	pma::Node3D* vStop = static_cast< pma::Node3D* >(next);

	TVector2 v0(p);
	v0 -= makeTVector2(vStart->Projection2D(view));

	TVector2 v1 = makeTVector2(vStop->Projection2D(view));
	v1 -= makeTVector2(vStart->Projection2D(view));
	
	TVector3 v3d = makeTVector3(vStop->Point3D());
	v3d -= makeTVector3(vStart->Point3D());
	
	TVector3 v3dStart = makeTVector3(vStart->Point3D());
	TVector3 v3dStop = makeTVector3(vStop->Point3D());

	double v0Norm = v0.Mod();
	double v1Norm = v1.Mod();

	TVector3 result(0, 0, 0);
	double eps = 1.0E-6; // 0.01mm
	if (v1Norm > eps)
	{
		double mag = v0Norm * v1Norm;
		double cosine = 0.0;
		if (mag != 0.0) cosine = v0 * v1 / mag;
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

		result = v3dStart;
		result += v3dStop;
		result *= 0.5;
	}
	return makePoint3D(result);
}

pma::Vector3D_t pma::Segment3D::GetUnconstrainedProj3D(const TVector2& p2d, unsigned int view) const
{
	pma::Node3D* vStart = static_cast< pma::Node3D* >(prev);
	pma::Node3D* vStop = static_cast< pma::Node3D* >(next);

	TVector2 v0(p2d);
	v0 -= makeTVector2(vStart->Projection2D(view));

	TVector2 v1 = makeTVector2(vStop->Projection2D(view));
	v1 -= makeTVector2(vStart->Projection2D(view));

	TVector3 v3d = makeTVector3(vStop->Point3D());
	v3d -= makeTVector3(vStart->Point3D());

	double v0Norm = v0.Mod();
	double v1Norm = v1.Mod();

	double eps = 1.0E-6; // 0.01mm
	if (v1Norm > eps)
	{
		double mag = v0Norm * v1Norm;
		double cosine = 0.0;
		if (mag != 0.0) cosine = v0 * v1 / mag;
		double b = v0Norm * cosine / v1Norm;

		return Vector3D_t(vStart->Point3D() + makeVector3D(v3d * b));
	}
	else // segment 2D projection is almost a point
	{
		mf::LogWarning("pma::Segment3D") << "Short segment projection." << std::endl;

		v3d = makeTVector3(vStart->Point3D());
		v3d += makeTVector3(vStop->Point3D());
		v3d *= 0.5;

		return makeVector3D(v3d);
	}
}

void pma::Segment3D::SetProjection(pma::Hit3D& h) const
{
	pma::Node3D* vStart = static_cast< pma::Node3D* >(prev);
	pma::Node3D* vStop = static_cast< pma::Node3D* >(next);

	auto /* const & */ pointStart = makeTVector3(vStart->Point3D());
	auto /* const & */ pointStop = makeTVector3(vStop->Point3D());

	auto const & projStart = vStart->Projection2D(h.View2D());
	auto const & projStop = vStop->Projection2D(h.View2D());

	pma::Vector2D v0(
		h.Point2D().X() - projStart.X(),
		h.Point2D().Y() - projStart.Y());

	pma::Vector2D v1(
		projStop.X() - projStart.X(),
		projStop.Y() - projStart.Y());

	pma::Vector3D v3d(
		pointStop.X() - pointStart.X(),
		pointStop.Y() - pointStart.Y(),
		pointStop.Z() - pointStart.Z());

	double v0Norm = sqrt(v0.Mag2());
	double v1Norm = sqrt(v1.Mag2());

	double eps = 1.0E-6; // 0.01mm
	if (v1Norm > eps)
	{
		double mag = v0Norm * v1Norm;
		double cosine = 0.0;
		if (mag != 0.0) cosine = v0.Dot(v1) / mag;
		double b = v0Norm * cosine / v1Norm;

		pma::Vector2D p(projStart.X(), projStart.Y());
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

		h.SetProjection(
			0.5 * (projStart.X() + projStop.X()),
			0.5 * (projStart.Y() + projStop.Y()), 0.0F);

		h.SetPoint3D(
			0.5 * (pointStart.X() + pointStop.X()),
			0.5 * (pointStart.Y() + pointStop.Y()),
			0.5 * (pointStart.Z() + pointStop.Z()));
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


double pma::Segment3D::GetDist2(const TVector3& psrc, const TVector3& p0, const TVector3& p1)
{
	pma::Vector3D v0(psrc.X() - p0.X(), psrc.Y() - p0.Y(), psrc.Z() - p0.Z());
	pma::Vector3D v1(p1.X() - p0.X(), p1.Y() - p0.Y(), p1.Z() - p0.Z());

	pma::Vector3D v2(psrc.X() - p1.X(), psrc.Y() - p1.Y(), psrc.Z() - p1.Z());
	pma::Vector3D v3(v1); v3 *= -1.0;

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

double pma::Segment3D::GetDist2(const TVector2& psrc, const TVector2& p0, const TVector2& p1)
{
	pma::Vector2D v0(psrc.X() - p0.X(), psrc.Y() - p0.Y());
	pma::Vector2D v1(p1.X() - p0.X(), p1.Y() - p0.Y());

	pma::Vector2D v2(psrc.X() - p1.X(), psrc.Y() - p1.Y());
	pma::Vector2D v3(v1); v3 *= -1.0;

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


double pma::Segment3D::GetDist2(const Point3D_t& psrc, const Point3D_t& p0, const Point3D_t& p1) {
  return GetDist2(makeTVector3(psrc), makeTVector3(p0), makeTVector3(p1));
}

double pma::Segment3D::GetDist2(const Point2D_t& psrc, const Point2D_t& p0, const Point2D_t& p1) {
  return GetDist2(makeTVector2(psrc), makeTVector2(p0), makeTVector2(p1));
}
