/**
 *  @file   PmaNode3D.h
 *
 *  @author D.Stefan and R.Sulej
 * 
 *  @brief  Implementation of the Projection Matching Algorithm
 *
 *          3D track node. See PmaTrack3D.h file for details.
 */

#ifndef PmaNode3D_h
#define PmaNode3D_h

#include "larreco/RecoAlg/PMAlg/PmaElement3D.h"
#include "larreco/RecoAlg/PMAlg/SortedObjects.h"

#include "larreco/RecoAlg/PMAlg/GeomDefs.h"

// #include "larcore/Geometry/Geometry.h"

namespace pma
{
	class Node3D;
}

class pma::Node3D : public pma::Element3D, public pma::SortedBranchBase
{
public:
	Node3D(void);
	Node3D(const Point3D_t& p3d, unsigned int tpc, unsigned int cryo, bool vtx = false);
	virtual ~Node3D(void) {}

	Point3D_t const & Point3D(void) const { return fPoint3D; }

	/// Returns true if the new position was accepted; returns false if the new position
	/// was trimmed to fit insite TPC volume + fMargin.
	bool SetPoint3D(const Point3D_t& p3d);

	Point2D_t const& Projection2D(unsigned int view) const { return fProj2D[view]; }

	double GetDistToWall(void) const;

	/// Check if p3d is in the same TPC as the node.
	bool SameTPC(const Point3D_t& p3d, float margin = 0.0F) const;

	/// Belongs to more than one track?
	bool IsBranching(void) const;

	/// Is the first/last in this TPC?
	bool IsTPCEdge(void) const;

	/// Check fIsVertex flag.
	bool IsVertex(void) const { return fIsVertex; }
	void SetVertex(bool state) { fIsVertex = state; }
	void SetVertexToBranching(bool setAllNodes)
	{
		if (setAllNodes || !fIsVertex) fIsVertex = IsBranching();
	}

	std::vector< pma::Track3D* > GetBranches(void) const;

	/// Distance [cm] from the 3D point to the point 3D.
	virtual double GetDistance2To(const Point3D_t& p3d) const override;

	/// Distance [cm] from the 2D point to the object's 2D projection in one of wire views.
	virtual double GetDistance2To(const Point2D_t& p2d, unsigned int view) const override;

	/// In case of a node it is simply 3D position of the node.
	virtual Vector3D_t GetUnconstrainedProj3D(const Point2D_t& , unsigned int) const override { return Vector3D_t(fPoint3D); }

	/// Set hit 3D position and its 2D projection to the vertex.
	virtual void SetProjection(pma::Hit3D& h) const override;

	/// Squared sum of half-lengths of connected 3D segments
	/// (used in the vertex position optimization).
	virtual double Length2(void) const override;

	/// Cosine of 3D angle between connected segments.
	double SegmentCos(void) const;
	/// Cosine of 2D angle (in plane parallel to wire planes) between connected segments.
	/// Should be changed / generalized for horizontal wire planes (e.g. 2-phase LAr).
	double SegmentCosWirePlane(void) const;
	/// Cosine of 2D angle (in horizontal plane, parallel to drift) between connected segments.
	/// Should be changed / generalized for horizontal wire planes (e.g. 2-phase LAr).
	double SegmentCosTransverse(void) const;

	/// Objective function minimized during oprimization.
	double GetObjFunction(float penaltyValue, float endSegWeight) const;

	/// Optimize vertex 3D position with given penalty on connected
	/// segments angle and weight assigned to the outermost segments.
	/// Only MSE is used in case of branching nodes.
	void Optimize(float penaltyValue, float endSegWeight);

	virtual void ClearAssigned(pma::Track3D* trk = 0) override;

	/// Set allowed node position margin around TPC.
	static void SetMargin(double m) { if (m >= 0.0) fMargin = m; }

private:
	/// Returns true if node position was trimmed to its TPC volume + fMargin
	bool LimitPoint3D(void);
	void UpdateProj2D(void);

	double EndPtCos2Transverse(void) const;
	double PiInWirePlane(void) const;
	double PenaltyInWirePlane(void) const;

	double Pi(float endSegWeight, bool doAsymm) const;
	double Penalty(float endSegWeight) const;
	double Mse(void) const;

	double MakeGradient(float penaltyValue, float endSegWeight);
	double StepWithGradient(float alfa, float tol, float penalty, float weight);

	art::ServiceHandle<geo::Geometry> fGeom;

	double fMinX, fMaxX, fMinY, fMaxY, fMinZ, fMaxZ; // TPC boundaries to limit the node position (+margin)
	double fWirePitch[3];                            // TPC params to scale do [cm] domain

	Point3D_t fPoint3D;       // node position in 3D space in [cm]
	Point2D_t fProj2D[3];     // node projections to 2D views, scaled to [cm], updated on each change of 3D position

	Vector3D_t fGradient;
	bool fIsVertex;          // no penalty on segments angle if branching or kink detected

	static bool fGradFixed[3];
	static double fMargin;
};

#endif

