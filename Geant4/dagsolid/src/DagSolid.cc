//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration and of QinetiQ Ltd,   *
// * subject to DEFCON 705 IPR conditions.                            *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: DagSolid.hh,v 1.10 2010/12/10 16:30:13 gunter Exp $
// GEANT4 tag $Name: geant4-09-05 $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              DagSolid.hh
//
// Date:                27/03/2014
// Author:              A. Davis M. C. Han, C. H. Kim, J. H. Jeong, Y. S. Yeom, S.
//                      Kim, Paul. P. H. Wilson, J. Apostolakis
// Organisation:        The University of Wisconsin-Madison, USA & Hanyang Univ., KR
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 27 March 2014, A Davis, UW - Updated description text
// 31 October 2010, J. H. Jeong, Hanyang Univ., KR
//  - Created.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//#include "G4TessellatedSolid.hh"
#include "DagSolid.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4TriangularFacet.hh"
#include "G4VFacet.hh"
#include "G4TessellatedSolid.hh"

#include <iostream>

#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/CartVect.hpp"
#include "DagMC.hpp"
using namespace moab;

#include "DagSolid.hh"

//#define G4SPECSDEBUG 1
///////////////////////////////////////////////////////////////////////////////
//
// Standard contructor has blank name and defines no facets.
//
DagSolid::DagSolid ()
  : G4TessellatedSolid("dummy"), cubicVolume(0.), surfaceArea(0.)
{

  geometryType = "DagSolid";
  
  xMinExtent =  kInfinity;
  xMaxExtent = -kInfinity;
  yMinExtent =  kInfinity;
  yMaxExtent = -kInfinity;
  zMinExtent =  kInfinity;
  zMaxExtent = -kInfinity;
//G4TessellatedSolid
//  SetRandomVectors();

  // SetRndVectors();

}


///////////////////////////////////////////////////////////////////////////////
//
// Alternative constructor. Simple define name and geometry type - no facets
// to detine.
//
DagSolid::DagSolid (const G4String &name, DagMC* dagmc, int volID)
  : G4TessellatedSolid(name), cubicVolume(0.), surfaceArea(0.)
{
  geometryType = "DagSolid";
  Myname=name;
  
  fdagmc=dagmc;
  fvolID=volID;
  fvolEntity = fdagmc->entity_by_index(3, volID);

  xMinExtent =  kInfinity;
  xMaxExtent = -kInfinity;
  yMinExtent =  kInfinity;
  yMaxExtent = -kInfinity;
  zMinExtent =  kInfinity;
  zMaxExtent = -kInfinity;

  int num_entities;
  std::vector<EntityHandle> surfs;  
  std::vector<EntityHandle> tris;  
  const EntityHandle *tri_conn;

  std::vector<CartVect> coords(3);
  G4ThreeVector vertex[3];
  int n_verts=0;
  
  Interface* moab = dagmc->moab_instance();
  moab->get_child_meshsets(fvolEntity, surfs, 1 );

  for(unsigned i=0 ; i<surfs.size() ; i++)
    {
      My_sulf_hit=surfs[i]; 
      moab->get_number_entities_by_type( surfs[i], MBTRI, num_entities);
      G4cout<<"Number of triangles = "<<num_entities<<" in surface index: "<<fdagmc->index_by_handle(surfs[i])<<G4endl;
      G4cout<<"please wait for visualization... "<<G4endl;
    
      moab->get_entities_by_type( surfs[i], MBTRI, tris);

      for (unsigned j=0 ; j <tris.size() ; j++) 
	{ 
	  moab->get_connectivity( tris[j], tri_conn, n_verts );    
	  moab->get_coords( tri_conn, n_verts, coords[0].array() );

	  //	  G4cout<<"add facet for vis = "<<coords[0]<<" "<<coords[1]<<" "<<coords[2]<<G4endl;

	  vertex[0] = G4ThreeVector( coords[0][0], coords[0][1], coords[0][2] );
	  vertex[1] = G4ThreeVector( coords[1][0], coords[1][1], coords[1][2] );
	  vertex[2] = G4ThreeVector( coords[2][0], coords[2][1], coords[2][2] );
      
	  G4TriangularFacet *facet = new G4TriangularFacet (vertex[0], vertex[1], vertex[2], ABSOLUTE);
	  AddFacet((G4VFacet*)facet);

	  for (G4int k=0 ; k < 3 ; k++) 
	    { 	
	      if ( vertex[k].x() < xMinExtent ) xMinExtent = vertex[k].x();
	      if ( vertex[k].x() > xMaxExtent ) xMaxExtent = vertex[k].x();
	      if ( vertex[k].y() < yMinExtent ) yMinExtent = vertex[k].y();
	      if ( vertex[k].y() > yMaxExtent ) yMaxExtent = vertex[k].y();
	      if ( vertex[k].z() < zMinExtent ) zMinExtent = vertex[k].z();
	      if ( vertex[k].z() > zMaxExtent ) zMaxExtent = vertex[k].z();      
	    }
	}
      tris.clear();
    }
  

 //SetRandomVectorSet();

  SetSolidClosed(true);
  
//  G4cout<<"Number Of Facets = "<<GetNumberOfFacets() <<G4endl;
//  G4cout <<"maximum point = "<< xMinExtent <<" "<< yMinExtent <<" "<< zMinExtent << G4endl
//         <<"minimum point = "<< xMaxExtent <<" "<< yMaxExtent <<" "<< zMaxExtent << G4endl;

}



///////////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
DagSolid::DagSolid( __void__& a )
  : G4TessellatedSolid(a), 
    geometryType("DagSolid"), cubicVolume(0.), surfaceArea(0.),
    xMinExtent(0.), xMaxExtent(0.),
    yMinExtent(0.), yMaxExtent(0.), 
    zMinExtent(0.), zMaxExtent(0.)
{
  //SetRandomVectorSet();
}



///////////////////////////////////////////////////////////////////////////////
//
// Destructor.
//
DagSolid::~DagSolid ()
{}


///////////////////////////////////////////////////////////////////////////////
//
// EInside DagSolid::Inside (const G4ThreeVector &p) const
//
EInside DagSolid::Inside (const G4ThreeVector &p) const
{
  G4double point[3]={p.x(), p.y(), p.z()};
  double u = rand();
  double v = rand();
  double w = rand();
  const double magnitude = sqrt( u*u + v*v + w*w );
  u /= magnitude;
  v /= magnitude;
  w /= magnitude;

  G4double direction[3]={u,v,w};

  int result;
  ErrorCode ec;
  ec = fdagmc->point_in_volume(fvolEntity, point, result,
			       direction);  // if uvw is not given, this function generate uvw ran

  if (ec != MB_SUCCESS)
    {
      G4cout << "failed to get point in volume" << std::endl;
      exit(1);
    }

  if ( result == 1 ) // point is contained within fvolEntity
    {
      return kInside;
    }
  else
    // point may still be within +kCarTolerance*0.5 of boundary
    {
      double minDist;
      ec = fdagmc->closest_to_location(fvolEntity, point, minDist);
      if( ec != MB_SUCCESS)
	{
	  G4cout << "failed to determine closed to location" << G4endl;
	  exit(1);
	}
      if (minDist <= 0.5*kCarTolerance) 
	return kSurface;
      else
	return kOutside;
    }


  // fast check to rule out points that are clearly beyond the volume
  /*
  if ( p.x() < xMinExtent - kCarTolerance || p.x() > xMaxExtent + kCarTolerance ||
       p.y() < yMinExtent - kCarTolerance || p.y() > yMaxExtent + kCarTolerance ||
       p.z() < zMinExtent - kCarTolerance || p.z() > zMaxExtent + kCarTolerance )
  {
    return kOutside;
  }  
 

  //  G4double minDist = kInfinity;
  G4double point[3]={p.x(), p.y(), p.z()};

  ErrorCode ec = fdagmc->closest_to_location(fvolEntity, point, minDist);
  if( ec != MB_SUCCESS)
    {
      exit(1);
    }

  if (minDist <= 0.5*kCarTolerance) 
    return kSurface;

  double u = rand();
  double v = rand();
  double w = rand();
  const double magnitude = sqrt( u*u + v*v + w*w );
  u /= magnitude;
  v /= magnitude;
  w /= magnitude;

  G4double direction[3]={u,v,w};

  int result;
  ec = fdagmc->point_in_volume(fvolEntity, point, result,
			       direction);  // if uvw is not given, this function generate uvw randomly
  G4cout << "inside: result =  " << result << std::endl;
  if (ec != MB_SUCCESS)
    {
      G4cout << "failed to get point in volume" << std::endl;
      exit(1);
    }

  if ( result == 0 )
    {
      G4cout << "point_in_cell not in vol " << fvolEntity << " " << p << G4endl;
      return kOutside;
    }
  else if (result == 1)
    {
      G4cout << "point_in_cell in vol " << fvolEntity << " " << p << G4endl;
      return kInside;
    }
  else
    {
      G4cout << "pint is nowhere" << G4endl;
      exit(1);
    }
  */
}

///////////////////////////////////////////////////////////////////////////////
//
// G4ThreeVector DagSolid::SurfaceNormal (const G4ThreeVector &p) const
//
// Return the outwards pointing unit normal of the shape for the
// surface closest to the point at offset p.

G4ThreeVector DagSolid::SurfaceNormal (const G4ThreeVector &p) const
{

  G4double ang[3]={0,0,1};
  G4double position[3]={p.x(),p.y(),p.z()};

  fdagmc->get_angle(My_sulf_hit, position, ang);

  G4ThreeVector normal = G4ThreeVector(ang[0],ang[1],ang[2]);

  //G4cout<<"SurfaceNormal = "<<normal<<G4endl;

  return normal;
}






///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v)
//
G4double DagSolid::DistanceToIn (const G4ThreeVector &p,
  const G4ThreeVector &v) const
{

  G4double minDist = kInfinity;
  G4double position[3] = {p.x(),p.y(),p.z()};
  G4double dir[3] = {v.x(),v.y(),v.z()};
  EntityHandle next_surf;
  G4double distance;

  DagMC::RayHistory history;
  
  // perform the ray fire with modified dag call
  fdagmc->ray_fire(fvolEntity,position,dir,next_surf,distance,&history,0,-1);
  history.reset();
  
  if ( next_surf == 0 ) // no intersection
    return kInfinity;
  else if ( -kCarTolerance*0.5 >= distance && distance <= kCarTolerance*0.5 )
    return 0.0;
  else
    return distance;

}



///////////////////////////////////////////////////////////////////////////////
//

// G4double DistanceToIn(const G4ThreeVector& p)
//
// Calculate distance to nearest surface of shape from an outside point p.

G4double DagSolid::DistanceToIn (const G4ThreeVector &p) const
{

  G4double minDist = kInfinity;
  G4double point[3]={p.x(), p.y(), p.z()};

  
  fdagmc->closest_to_location(fvolEntity, point, minDist);
  
  if ( minDist <= kCarTolerance*0.5 )
    return kInfinity;
  else
    {
      return minDist;
    }
}


///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
//                        const G4bool calcNorm=false,
//                        G4bool *validNorm=0, G4ThreeVector *n=0);
//
// Return distance along the normalised vector v to the shape, from a
// point at an offset p inside or on the surface of the shape. 
// Intersections with surfaces, when the point is not greater
// than kCarTolerance/2 from a surface, must be ignored.
//     If calcNorm is true, then it must also set validNorm to either
//     * true, if the solid lies entirely behind or on the exiting
//        surface. Then it must set n to the outwards normal vector
//        (the Magnitude of the vector is not defined).
//     * false, if the solid does not lie entirely behind or on the
//       exiting surface.
// If calcNorm is false, then validNorm and n are unused.

G4double DagSolid::DistanceToOut (const G4ThreeVector &p,
                    const G4ThreeVector &v, const G4bool calcNorm,
                          G4bool *validNorm, G4ThreeVector *n) const
{
  G4double minDist = kInfinity;
  double position[3]={p.x(),p.y(),p.z()};
  double dir[3]={v.x(),v.y(),v.z()};

  EntityHandle next_surf;
  double next_dist;
  DagMC::RayHistory history;

  fdagmc->ray_fire(fvolEntity,position,dir,next_surf,next_dist,&history,0,1);
  history.reset();

  if(next_surf == 0 )
    return kInfinity;
  
  if (calcNorm)
    {
      *n         = SurfaceNormal(p+minDist*v);
    }

  if (next_dist < minDist )
    minDist = next_dist;
 
          
  if (minDist > 0.0 && minDist <= 0.5*kCarTolerance ) 
    {
      return kInfinity;
    }
  else if ( minDist < kInfinity)
    {
      return minDist;
    }
  else
    {
      return kInfinity;
    }
}

///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToOut(const G4ThreeVector& p)
//
// Calculate distance to nearest surface of shape from an inside point. 

G4double DagSolid::DistanceToOut (const G4ThreeVector &p) const
{
  G4double minDist = kInfinity;
  G4double point[3]={p.x(), p.y(), p.z()};

  fdagmc->closest_to_location(fvolEntity, point, minDist);
  if ( minDist < kCarTolerance/2.0 )
    return kInfinity;
  else
    {
      return minDist;
    }
}

///////////////////////////////////////////////////////////////////////////////
//
// G4GeometryType GetEntityType() const;
//
// Provide identification of the class of an object (required for persistency and STEP interface).
//
G4GeometryType DagSolid::GetEntityType () const
{
  return geometryType;
}

///////////////////////////////////////////////////////////////////////////////
//
void DagSolid::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddSolid (*this);
}


std::ostream &DagSolid::StreamInfo(std::ostream &os) const
{
  os << G4endl;
  os << "Geometry Type    = " << geometryType  << G4endl;
  return os;
}


///////////////////////////////////////////////////////////////////////////////
//
// CalculateExtent
//
// Based on correction provided by Stan Seibert, University of Texas.
//
G4bool
DagSolid::CalculateExtent(const EAxis pAxis,
                          const G4VoxelLimits& pVoxelLimit,
                          const G4AffineTransform& pTransform,
                                G4double& pMin, G4double& pMax) const
{


// Calculate the minimum and maximum extent of the solid, 
// when under the specified transform, and within the specified limits. 
// If the solid is not intersected by the region, return false, else return true.


    G4ThreeVector minExtent(xMinExtent, yMinExtent, zMinExtent);
    G4ThreeVector maxExtent(xMaxExtent, yMaxExtent, zMaxExtent);


    // Check for containment and clamp to voxel boundaries
    for (G4int axis=G4ThreeVector::X; axis < G4ThreeVector::SIZE; axis++)
    {
      EAxis geomAxis = kXAxis; // G4 geom classes use different index type
      switch(axis)
      {
        case G4ThreeVector::X: geomAxis = kXAxis; break;
        case G4ThreeVector::Y: geomAxis = kYAxis; break;
        case G4ThreeVector::Z: geomAxis = kZAxis; break;
      }
      G4bool isLimited = pVoxelLimit.IsLimited(geomAxis);
      G4double voxelMinExtent = pVoxelLimit.GetMinExtent(geomAxis);
      G4double voxelMaxExtent = pVoxelLimit.GetMaxExtent(geomAxis);

      if (isLimited)
      {
        if ( minExtent[axis] > voxelMaxExtent+kCarTolerance ||
             maxExtent[axis] < voxelMinExtent-kCarTolerance    )
        {
          return false ;
        }
        else
        {
          if (minExtent[axis] < voxelMinExtent)
          {
            minExtent[axis] = voxelMinExtent ;
          }
          if (maxExtent[axis] > voxelMaxExtent)
          {
            maxExtent[axis] = voxelMaxExtent;
          }
        }
      }
    }

    // Convert pAxis into G4ThreeVector index
    G4int vecAxis=0;
    switch(pAxis)
    {
      case kXAxis: vecAxis = G4ThreeVector::X; break;
      case kYAxis: vecAxis = G4ThreeVector::Y; break;
      case kZAxis: vecAxis = G4ThreeVector::Z; break;
      default: break;
    }

    pMin = minExtent[vecAxis] - kCarTolerance;
    pMax = maxExtent[vecAxis] + kCarTolerance;

    return true;
}

G4double DagSolid::GetMinXExtent () const
  {return xMinExtent;}

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetMaxXExtent () const
  {return xMaxExtent;}

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetMinYExtent () const
  {return yMinExtent;}

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetMaxYExtent () const
  {return yMaxExtent;}

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetMinZExtent () const
  {return zMinExtent;}

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetMaxZExtent () const
  {return zMaxExtent;}

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetCubicVolume()
{
  G4double result;
  fdagmc->measure_volume(fvolEntity, result);
  return result;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double DagSolid::GetSurfaceArea()
{
  G4double result;
  fdagmc->measure_area(fvolEntity, result);
  return result;
}
