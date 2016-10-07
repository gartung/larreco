/**
 * @file   LegacyGeomDefs.h
 * @brief  Definition of some geometry-related data structures
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   September 8, 2016
 * 
 * This is a pure header: no implementation file is required.
 */

#ifndef LARRECO_RECOALG_PMALG_LEGACYGEOMDEFS_H
#define LARRECO_RECOALG_PMALG_LEGACYGEOMDEFS_H

#include "larreco/RecoAlg/PMAlg/GeomDefs.h"

// ROOT
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


#endif // LARRECO_RECOALG_PMALG_LEGACYGEOMDEFS_H
