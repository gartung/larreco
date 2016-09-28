/**
 * @file   GeomDefs.h
 * @brief  Definition of some geometry-related data structures
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   September 8, 2016
 * 
 * This is a pure header: no implementation file is required.
 */

#ifndef LARRECO_RECOALG_PMALG_GEOMDEFS_H
#define LARRECO_RECOALG_PMALG_GEOMDEFS_H

// ROOT
#include "Math/GenVector/Cartesian2D.h"
#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/PositionVector2D.h"
#include "Math/GenVector/PositionVector3D.h"
#include "Math/GenVector/DisplacementVector2D.h"
#include "Math/GenVector/DisplacementVector3D.h"


namespace pma {
  
  /// Type of real number used in the computations (double precision)
  using Real = double;
  
  /// coordinate system used for our 3D objects
  using CoordSystem3D = ROOT::Math::Cartesian3D<Real>;
  
  /// coordinate system used for our 2D objects
  using CoordSystem2D = ROOT::Math::Cartesian2D<Real>;
  
  /// A point in 2D space
  using Point2D_t = ROOT::Math::PositionVector2D<CoordSystem2D>;
  
  /// A point in 3D space
  using Point3D_t = ROOT::Math::PositionVector3D<CoordSystem3D>;
  
  
  /// A vector in 2D space
  using Vector2D_t = ROOT::Math::DisplacementVector2D<CoordSystem2D>;
  
  /// A vector in 3D space
  using Vector3D_t = ROOT::Math::DisplacementVector3D<CoordSystem3D>;
  
} // namespace pma


// FIXME DELME
#include "TVector.h"
#include "TVector3.h"

inline TVector2 makeTVector2(pma::Point2D_t const& point)
  { return { point.X(), point.Y() }; }

inline TVector2 makeTVector2(pma::Vector2D_t const& vector)
  { return { vector.X(), vector.Y() }; }

inline TVector3 makeTVector3(pma::Point3D_t const& point)
  { return { point.X(), point.Y(), point.Z() }; }

inline TVector3 makeTVector3(pma::Vector3D_t const& vector)
  { return { vector.X(), vector.Y(), vector.Z() }; }

inline pma::Vector2D_t makeVector2D(TVector2 const& vector)
  { return { vector.X(), vector.Y() }; }

inline pma::Point2D_t makePoint2D(TVector2 const& point)
  { return { point.X(), point.Y() }; }

inline pma::Vector3D_t makeVector3D(TVector3 const& vector)
  { return { vector.X(), vector.Y(), vector.Z() }; }

inline pma::Point3D_t makePoint3D(TVector3 const& point)
  { return { point.X(), point.Y(), point.Z() }; }


#endif // LARRECO_RECOALG_PMALG_GEOMDEFS_H
