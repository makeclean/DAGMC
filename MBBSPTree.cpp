/*
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2008 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

/**\file MBBSPTree.cpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2008-05-13
 */

#include "MBBSPTree.hpp"
#include "MBGeomUtil.hpp"
#include "MBRange.hpp"
#include "MBInternals.hpp"
#include "BSPTreePoly.hpp"

#include <assert.h>
#include <algorithm>
#include <limits>

#if defined(_MSC_VER) || defined(__MINGW32__)
#  include <float.h>
#  define finite(A) _finite(A)
#elif defined(HAVE_IEEEFP_H)
#  include <ieeefp.h>
#endif

#define MB_BSP_TREE_DEFAULT_TAG_NAME "BSPTree"

static void corners_from_box( const double box_min[3],
                              const double box_max[3],
                              double corners[8][3] )
{
  const double* ranges[] = { box_min, box_max };
  for (int z = 0; z < 2; ++z) {
    corners[4*z  ][0] = box_min[0];
    corners[4*z  ][1] = box_min[1];
    corners[4*z  ][2] = ranges[z][2];

    corners[4*z+1][0] = box_max[0];
    corners[4*z+1][1] = box_min[1];
    corners[4*z+1][2] = ranges[z][2];

    corners[4*z+2][0] = box_max[0];
    corners[4*z+2][1] = box_max[1];
    corners[4*z+2][2] = ranges[z][2];

    corners[4*z+3][0] = box_min[0];
    corners[4*z+3][1] = box_max[1];
    corners[4*z+3][2] = ranges[z][2];
  }
}

// assume box has planar sides
// test if point is contained in box
static bool point_in_box( const double corners[8][3],
                          const double point[3] )
{
  const unsigned side_verts[6][3] = { { 0, 3, 1 },
                                      { 4, 5, 7 },
                                      { 0, 1, 4 },
                                      { 1, 2, 5 },
                                      { 2, 3, 6 },
                                      { 3, 0, 7 } };
    // If we assume planar sides, then the box is the intersection
    // of 6 half-spaces defined by the planes of the sides.
  const MBCartVect pt(point);
  for (unsigned s = 0; s < 6; ++s) {
    MBCartVect v0( corners[side_verts[s][0]] );
    MBCartVect v1( corners[side_verts[s][1]] );
    MBCartVect v2( corners[side_verts[s][2]] );
    MBCartVect N = (v1 - v0) * (v2 - v0);
    if ((v0 - pt) % N < 0)
      return false;
  }
  return true;
}

void MBBSPTree::Plane::set( const double pt1[3], const double pt2[3], const double pt3[3] )
{
  const double v1[] = { pt2[0] - pt1[0], pt2[1] - pt1[1], pt2[2] - pt1[2] };
  const double v2[] = { pt3[0] - pt1[0], pt3[1] - pt1[1], pt3[2] - pt1[2] };
  const double norm[] = { v1[1]*v2[2] - v1[2]*v2[1],
                          v1[2]*v2[0] - v1[0]*v2[2],
                          v1[0]*v2[1] - v1[1]*v2[0] };
  set( norm, pt1 );
}

MBErrorCode MBBSPTree::init_tags( const char* tagname )
{
  if (!tagname) 
    tagname = MB_BSP_TREE_DEFAULT_TAG_NAME;
  
  std::string rootname(tagname);
  rootname += "_box";
  
  MBErrorCode rval = moab()->tag_create( tagname, 4*sizeof(double), MB_TAG_DENSE, MB_TYPE_DOUBLE, planeTag, 0, true );
  if (MB_SUCCESS != rval)
    planeTag = 0;
  else
    rval = moab()->tag_create( rootname.c_str(), 24*sizeof(double), MB_TAG_SPARSE, MB_TYPE_DOUBLE, rootTag, 0, true );
  if (MB_SUCCESS != rval)
    rootTag = 0;
  return rval;
}

MBBSPTree::MBBSPTree( MBInterface* mb, 
                      const char* tagname, 
                      unsigned set_flags )
  : mbInstance(mb), meshSetFlags(set_flags), cleanUpTrees(false)
{ init_tags( tagname ); }

MBBSPTree::MBBSPTree( MBInterface* mb, 
                      bool destroy_created_trees,
                      const char* tagname, 
                      unsigned set_flags )
  : mbInstance(mb), meshSetFlags(set_flags), cleanUpTrees(destroy_created_trees)
{ init_tags( tagname ); }

MBBSPTree::~MBBSPTree()
{
  if (!cleanUpTrees)
    return;
    
  while (!createdTrees.empty()) {
    MBEntityHandle tree = createdTrees.back();
      // make sure this is a tree (rather than some other, stale handle)
    const void* data_ptr = 0;
    MBErrorCode rval = moab()->tag_get_data( rootTag, &tree, 1, &data_ptr );
    if (MB_SUCCESS == rval)
      rval = delete_tree( tree );
    if (MB_SUCCESS != rval)
      createdTrees.pop_back();
  }
}

MBErrorCode MBBSPTree::set_split_plane( MBEntityHandle node, const Plane& p )
{ 
    // check for unit-length normal
  const double lensqr = p.norm[0]*p.norm[0] 
                      + p.norm[1]*p.norm[1] 
                      + p.norm[2]*p.norm[2];
  if (fabs(lensqr - 1.0) < std::numeric_limits<double>::epsilon())
    return moab()->tag_set_data( planeTag, &node, 1, &p ); 
    
  const double inv_len = 1.0/sqrt(lensqr);
  Plane p2(p);
  p2.norm[0] *= inv_len;
  p2.norm[1] *= inv_len;
  p2.norm[2] *= inv_len;
  p2.coeff   *= inv_len;
  
    // check for zero-length normal
  if (!finite(p2.norm[0]+p2.norm[1]+p2.norm[2]+p2.coeff))
    return MB_FAILURE;

    // store plane
  return moab()->tag_set_data( planeTag, &node, 1, &p2 ); 
}

MBErrorCode MBBSPTree::set_tree_box( MBEntityHandle root_handle,
                                     const double box_min[3],
                                     const double box_max[3] )
{
  double corners[8][3];
  corners_from_box( box_min, box_max, corners );
  return set_tree_box( root_handle, corners );
}

MBErrorCode MBBSPTree::set_tree_box( MBEntityHandle root_handle,
                                     const double corners[8][3] )
{
  return moab()->tag_set_data( rootTag, &root_handle, 1, corners );
}

MBErrorCode MBBSPTree::get_tree_box( MBEntityHandle root_handle,
                                     double corners[8][3] )
{
  return moab()->tag_get_data( rootTag, &root_handle, 1, corners );
}

MBErrorCode MBBSPTree::get_tree_box( MBEntityHandle root_handle,
                                     double corners[24] )
{
  return moab()->tag_get_data( rootTag, &root_handle, 1, corners );
}

MBErrorCode MBBSPTree::create_tree( MBEntityHandle& root_handle )
{
  const double min[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
  const double max[3] = {  HUGE_VAL,  HUGE_VAL,  HUGE_VAL };
  return create_tree( min, max, root_handle );
}

MBErrorCode MBBSPTree::create_tree( const double corners[8][3],
                                    MBEntityHandle& root_handle )
{
  MBErrorCode rval = moab()->create_meshset( meshSetFlags, root_handle );
  if (MB_SUCCESS != rval)
    return rval;
  
  rval = set_tree_box( root_handle, corners );
  if (MB_SUCCESS != rval) {
    moab()->delete_entities( &root_handle, 1 );
    root_handle = 0;
    return rval;
  }
  
  createdTrees.push_back( root_handle );
  return MB_SUCCESS;
}
                                    

MBErrorCode MBBSPTree::create_tree( const double box_min[3],
                                    const double box_max[3],
                                    MBEntityHandle& root_handle )
{
  double corners[8][3];
  corners_from_box( box_min, box_max, corners );
  return create_tree( corners, root_handle );
}

MBErrorCode MBBSPTree::delete_tree( MBEntityHandle root_handle )
{
  MBErrorCode rval;
  
  std::vector<MBEntityHandle> children, dead_sets, current_sets;
  current_sets.push_back( root_handle );
  while (!current_sets.empty()) {
    MBEntityHandle set = current_sets.back();
    current_sets.pop_back();
    dead_sets.push_back( set );
    rval = moab()->get_child_meshsets( set, children );
    if (MB_SUCCESS != rval)
      return rval;
    std::copy( children.begin(), children.end(), std::back_inserter(current_sets) );
    children.clear();
  }
  
  rval = moab()->tag_delete_data( rootTag, &root_handle, 1 );
  if (MB_SUCCESS != rval)
    return rval;
  
  createdTrees.erase(
    std::remove( createdTrees.begin(), createdTrees.end(), root_handle ),
    createdTrees.end() );
  return moab()->delete_entities( &dead_sets[0], dead_sets.size() );
}

MBErrorCode MBBSPTree::find_all_trees( MBRange& results )
{
  return moab()->get_entities_by_type_and_tag( 0, MBENTITYSET, 
                                               &rootTag, 0, 1,
                                               results );
}

MBErrorCode MBBSPTree::get_tree_iterator( MBEntityHandle root,
                                          MBBSPTreeIter& iter )
{
  MBErrorCode rval = iter.initialize( this, root );
  if (MB_SUCCESS != rval)
    return rval;
  return iter.step_to_first_leaf( MBBSPTreeIter::LEFT );
}

MBErrorCode MBBSPTree::get_tree_end_iterator( MBEntityHandle root,
                                          MBBSPTreeIter& iter )
{
  MBErrorCode rval = iter.initialize( this, root );
  if (MB_SUCCESS != rval)
    return rval;
  return iter.step_to_first_leaf( MBBSPTreeIter::RIGHT );
}

MBErrorCode MBBSPTree::split_leaf( MBBSPTreeIter& leaf,
                                   Plane plane,
                                   MBEntityHandle& left,
                                   MBEntityHandle& right )
{
  MBErrorCode rval;
  
  rval = moab()->create_meshset( meshSetFlags, left );
  if (MB_SUCCESS != rval)
    return rval;
  
  rval = moab()->create_meshset( meshSetFlags, right );
  if (MB_SUCCESS != rval) {
    moab()->delete_entities( &left, 1 );
    return rval;
  }
  
  if (MB_SUCCESS != set_split_plane( leaf.handle(), plane ) ||
      MB_SUCCESS != moab()->add_child_meshset( leaf.handle(), left ) ||
      MB_SUCCESS != moab()->add_child_meshset( leaf.handle(), right) ||
      MB_SUCCESS != leaf.step_to_first_leaf(MBBSPTreeIter::LEFT)) {
    MBEntityHandle children[] = { left, right };
    moab()->delete_entities( children, 2 );
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBBSPTree::split_leaf( MBBSPTreeIter& leaf, Plane plane )
{
  MBEntityHandle left, right;
  return split_leaf( leaf, plane, left, right );
}

MBErrorCode MBBSPTree::split_leaf( MBBSPTreeIter& leaf, 
                                   Plane plane,
                                   const MBRange& left_entities,
                                   const MBRange& right_entities )
{
  MBEntityHandle left, right, parent = leaf.handle();
  MBErrorCode rval = split_leaf( leaf, plane, left, right );
  if (MB_SUCCESS != rval)
    return rval;
  
  if (MB_SUCCESS == moab()->add_entities( left, left_entities ) &&
      MB_SUCCESS == moab()->add_entities(right,right_entities ) &&
      MB_SUCCESS == moab()->clear_meshset( &parent, 1 ))
    return MB_SUCCESS;
  
  moab()->remove_child_meshset( parent, left );
  moab()->remove_child_meshset( parent, right );
  MBEntityHandle children[] = { left, right };
  moab()->delete_entities( children, 2 );
  return MB_FAILURE;
}

MBErrorCode MBBSPTree::split_leaf( MBBSPTreeIter& leaf, Plane plane,
                                   const std::vector<MBEntityHandle>& left_entities,
                                   const std::vector<MBEntityHandle>& right_entities )
{
  MBEntityHandle left, right, parent = leaf.handle();
  MBErrorCode rval = split_leaf( leaf, plane, left, right );
  if (MB_SUCCESS != rval)
    return rval;
  
  if (MB_SUCCESS == moab()->add_entities( left, &left_entities[0], left_entities.size() ) &&
      MB_SUCCESS == moab()->add_entities(right,&right_entities[0],right_entities.size() ) &&
      MB_SUCCESS == moab()->clear_meshset( &parent, 1 ))
    return MB_SUCCESS;
  
  moab()->remove_child_meshset( parent, left );
  moab()->remove_child_meshset( parent, right );
  MBEntityHandle children[] = { left, right };
  moab()->delete_entities( children, 2 );
  return MB_FAILURE;
}

MBErrorCode MBBSPTree::merge_leaf( MBBSPTreeIter& iter )
{
  MBErrorCode rval;
  if (iter.depth() == 1) // at root
    return MB_FAILURE;
  
    // Move iter to parent
  iter.up();

    // Get sets to merge
  MBEntityHandle parent = iter.handle();
  iter.childVect.clear();
  rval = moab()->get_child_meshsets( parent, iter.childVect );
  if (MB_SUCCESS != rval)
    return rval;
    
    // Remove child links
  moab()->remove_child_meshset( parent, iter.childVect[0] );
  moab()->remove_child_meshset( parent, iter.childVect[1] );
  std::vector<MBEntityHandle> stack( iter.childVect );
  
    // Get all entities from children and put them in parent
  MBRange range;
  while (!stack.empty()) {
    MBEntityHandle h = stack.back();
    stack.pop_back();
    range.clear();
    rval = moab()->get_entities_by_handle( h, range );
    if (MB_SUCCESS != rval)
      return rval;
    rval = moab()->add_entities( parent, range );
    if (MB_SUCCESS != rval)
      return rval;
    
    iter.childVect.clear();
    moab()->get_child_meshsets( h, iter.childVect );
    if (!iter.childVect.empty()) {
     moab()->remove_child_meshset( h, iter.childVect[0] );
     moab()->remove_child_meshset( h, iter.childVect[1] );
     stack.push_back( iter.childVect[0] );
     stack.push_back( iter.childVect[1] );
    }
  
    rval = moab()->delete_entities( &h, 1 );
    if (MB_SUCCESS != rval)
      return rval;
  }
  
  return MB_SUCCESS;
}

  

MBErrorCode MBBSPTreeIter::initialize( MBBSPTree* tool,
                                       MBEntityHandle root,
                                       const double* point )
{
  treeTool = tool;
  mStack.clear();
  mStack.push_back( root );
  return MB_SUCCESS;
}


MBErrorCode MBBSPTreeIter::step_to_first_leaf( Direction direction )
{
  MBErrorCode rval;
  for (;;) {
    childVect.clear();
    rval = tool()->moab()->get_child_meshsets( mStack.back(), childVect );
    if (MB_SUCCESS != rval)
      return rval;
    if (childVect.empty()) // leaf
      break;
  
    mStack.push_back( childVect[direction] );
  }
  return MB_SUCCESS;
}

MBErrorCode MBBSPTreeIter::step( Direction direction )
{
  MBEntityHandle node, parent;
  MBErrorCode rval;
  const Direction opposite = static_cast<Direction>(1-direction);
  
    // If stack is empty, then either this iterator is uninitialized
    // or we reached the end of the iteration (and return 
    // MB_ENTITY_NOT_FOUND) already.
  if (mStack.empty())
    return MB_FAILURE;
    
    // Pop the current node from the stack.
    // The stack should then contain the parent of the current node.
    // If the stack is empty after this pop, then we've reached the end.
  node = mStack.back();
  mStack.pop_back();
  
  while(!mStack.empty()) {
      // Get data for parent entity
    parent = mStack.back();
    childVect.clear();
    rval = tool()->moab()->get_child_meshsets( parent, childVect );
    if (MB_SUCCESS != rval)
      return rval;
    
      // If we're at the left child
    if (childVect[opposite] == node) {
        // push right child on stack
      mStack.push_back( childVect[direction] );
        // descend to left-most leaf of the right child
      return step_to_first_leaf(opposite);
    }
    
      // The current node is the right child of the parent,
      // continue up the tree.
    assert( childVect[direction] == node );
    node = parent;
    mStack.pop_back();
  }
  
  return MB_ENTITY_NOT_FOUND;
}

MBErrorCode MBBSPTreeIter::up() 
{
  if (mStack.size() < 2)
    return MB_ENTITY_NOT_FOUND;
  mStack.pop_back();
  return MB_SUCCESS;
}

MBErrorCode MBBSPTreeIter::down( const MBBSPTree::Plane& plane, Direction dir ) 
{
  childVect.clear();
  MBErrorCode rval = tool()->moab()->get_child_meshsets( mStack.back(), childVect );
  if (MB_SUCCESS != rval)
    return rval;
  if (childVect.empty())
    return MB_ENTITY_NOT_FOUND;
  
  mStack.push_back( childVect[dir] );
  return MB_SUCCESS;
}

MBErrorCode MBBSPTreeIter::get_parent_split_plane( MBBSPTree::Plane& plane ) const
{
  if (mStack.size() < 2) // at tree root
    return MB_ENTITY_NOT_FOUND;
  
  MBEntityHandle parent = mStack[mStack.size()-2];
  return tool()->get_split_plane( parent, plane );
}

double MBBSPTreeIter::volume() const
{
  BSPTreePoly polyhedron;
  MBErrorCode rval = calculate_polyhedron( polyhedron );
  return MB_SUCCESS == rval ? polyhedron.volume() : -1.0;
}

bool MBBSPTreeIter::is_sibling( const MBBSPTreeIter& other_leaf ) const
{
  const size_t s = mStack.size();
  return (s > 1) && (s == other_leaf.mStack.size()) &&
         (other_leaf.mStack[s-2] == mStack[s-2]) &&
         other_leaf.handle() != handle();
}

bool MBBSPTreeIter::is_sibling( MBEntityHandle other_leaf ) const
{
  if (mStack.size() < 2 || other_leaf == handle())
    return false;
  MBEntityHandle parent = mStack[mStack.size()-2];
  childVect.clear();
  MBErrorCode rval = tool()->moab()->get_child_meshsets( parent, childVect );
  if (MB_SUCCESS != rval || childVect.size() != 2) {
    assert(false);
    return false;
  }
  return childVect[0] == other_leaf || childVect[1] == other_leaf;
}

bool MBBSPTreeIter::sibling_is_forward() const
{
  if (mStack.size() < 2) // if root
    return false;
  MBEntityHandle parent = mStack[mStack.size()-2];
  childVect.clear();
  MBErrorCode rval = tool()->moab()->get_child_meshsets( parent, childVect );
  if (MB_SUCCESS != rval || childVect.size() != 2) {
    assert(false);
    return false;
  }
  return childVect[0] == handle();
}  

MBErrorCode MBBSPTreeIter::calculate_polyhedron( BSPTreePoly& poly_out ) const
{
  MBErrorCode rval;
  
  assert( sizeof(MBCartVect) == 3*sizeof(double) );
  MBCartVect corners[8];
  rval = treeTool->get_tree_box( mStack.front(), corners[0].array() );
  if (MB_SUCCESS != rval)
    return rval;
  
  rval = poly_out.set( corners );
  if (MB_SUCCESS != rval)
    return rval;
  
  MBBSPTree::Plane plane;
  std::vector<MBEntityHandle>::const_iterator i = mStack.begin();
  std::vector<MBEntityHandle>::const_iterator here = mStack.end() - 1;
  while (i != here) {
    rval = treeTool->get_split_plane( *i, plane );
    if (MB_SUCCESS != rval)
      return rval;
    
    childVect.clear();
    rval = treeTool->moab()->get_child_meshsets( *i, childVect );
    if (MB_SUCCESS != rval)
      return rval;
    if (childVect.size() != 2)
      return MB_FAILURE;
      
    ++i;
    if (childVect[1] == *i)
      plane.flip();
    
    MBCartVect norm( plane.norm );
    poly_out.cut_polyhedron( norm, plane.coeff );
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBBSPTreeBoxIter::initialize( MBBSPTree* tool_ptr,
                                          MBEntityHandle root,
                                          const double* point )
{
  MBErrorCode rval = MBBSPTreeIter::initialize( tool_ptr, root );
  if (MB_SUCCESS != rval)
    return rval;
  
  tool()->get_tree_box( root, leafCoords );
  if (MB_SUCCESS != rval)
    return rval;

  if (point && !point_in_box( leafCoords, point ))
    return MB_ENTITY_NOT_FOUND;

  stackData.resize(1);
  return MB_SUCCESS;
}

MBBSPTreeBoxIter::SideBits
MBBSPTreeBoxIter::side_above_plane( const double hex_coords[8][3],
                                    const MBBSPTree::Plane& plane )
{
  unsigned result  = 0;
  for (unsigned i = 0; i < 8u; ++i) 
    result |= plane.above(hex_coords[i]) << i;
  return (MBBSPTreeBoxIter::SideBits)result;
}

MBBSPTreeBoxIter::SideBits
MBBSPTreeBoxIter::side_on_plane( const double hex_coords[8][3],
                                 const MBBSPTree::Plane& plane )
{
  unsigned result  = 0;
  for (unsigned i = 0; i < 8u; ++i) {
    bool on = plane.distance(hex_coords[i]) <= MBBSPTree::epsilon();
    result |= on << i;
  }
  return (MBBSPTreeBoxIter::SideBits)result;
}

static inline void copy_coords( const double src[3], double dest[3] )
{
  dest[0] = src[0];
  dest[1] = src[1];
  dest[2] = src[2];
}

MBErrorCode MBBSPTreeBoxIter::face_corners( const SideBits face,
                                            const double hex_corners[8][3],
                                            double face_corners[4][3] )
{
  switch (face) {
    case MBBSPTreeBoxIter::B0154:
      copy_coords( hex_corners[0], face_corners[0] );
      copy_coords( hex_corners[1], face_corners[1] );
      copy_coords( hex_corners[5], face_corners[2] );
      copy_coords( hex_corners[4], face_corners[3] );
      break;
    case MBBSPTreeBoxIter::B1265:
      copy_coords( hex_corners[1], face_corners[0] );
      copy_coords( hex_corners[2], face_corners[1] );
      copy_coords( hex_corners[6], face_corners[2] );
      copy_coords( hex_corners[5], face_corners[3] );
      break;
    case MBBSPTreeBoxIter::B2376:
      copy_coords( hex_corners[2], face_corners[0] );
      copy_coords( hex_corners[3], face_corners[1] );
      copy_coords( hex_corners[7], face_corners[2] );
      copy_coords( hex_corners[6], face_corners[3] );
      break;
    case MBBSPTreeBoxIter::B3047:
      copy_coords( hex_corners[3], face_corners[0] );
      copy_coords( hex_corners[0], face_corners[1] );
      copy_coords( hex_corners[4], face_corners[2] );
      copy_coords( hex_corners[7], face_corners[3] );
      break;
    case MBBSPTreeBoxIter::B3210:
      copy_coords( hex_corners[3], face_corners[0] );
      copy_coords( hex_corners[2], face_corners[1] );
      copy_coords( hex_corners[1], face_corners[2] );
      copy_coords( hex_corners[0], face_corners[3] );
      break;
    case MBBSPTreeBoxIter::B4567:
      copy_coords( hex_corners[4], face_corners[0] );
      copy_coords( hex_corners[5], face_corners[1] );
      copy_coords( hex_corners[6], face_corners[2] );
      copy_coords( hex_corners[7], face_corners[3] );
      break;
    default:
      return MB_FAILURE; // child is not a box
  }
  
  return MB_SUCCESS;

}

/** \brief Clip an edge using a plane
 *
 * Given an edge from keep_end_coords to cut_end_coords,
 * cut the edge using the passed plane, such that cut_end_coords
 * is updated with a new location on the plane, and old_coords_out
 * contains the original value of cut_end_coords.
 */
static inline
void plane_cut_edge( double old_coords_out[3],
                     const double keep_end_coords[3],
                     double cut_end_coords[3],
                     const MBBSPTree::Plane& plane )
{
  const MBCartVect start( keep_end_coords ), end( cut_end_coords );
  const MBCartVect norm( plane.norm );
  MBCartVect xsect_point;
  
  const MBCartVect m = end - start;
  const double t = -(norm % start + plane.coeff) / (norm % m);
  assert( t > 0.0 && t < 1.0 );
  xsect_point = start + t * m;
  
  end.get( old_coords_out );
  xsect_point.get( cut_end_coords );
}

/** Given the corners of a hexahedron in corners_input and a 
 *  plane, cut the hex with the plane, updating corners_input
 *  and storing the original,cut-off side of the hex in cut_face_out.
 *
 *  The portion of the hex below the plane is retained.  cut_face_out
 *  will contain the side of the hex that is entirely above the plane.
 *\return MB_FAILURE if plane/hex intersection is not a quadrilateral.
 */
static MBErrorCode plane_cut_box( double cut_face_out[4][3],
                                  double corners_inout[8][3],
                                  const MBBSPTree::Plane& plane )
{
  switch (MBBSPTreeBoxIter::side_above_plane( corners_inout, plane )) {
    case MBBSPTreeBoxIter::B0154:
      plane_cut_edge( cut_face_out[0], corners_inout[3], corners_inout[0], plane );
      plane_cut_edge( cut_face_out[1], corners_inout[2], corners_inout[1], plane );
      plane_cut_edge( cut_face_out[2], corners_inout[6], corners_inout[5], plane );
      plane_cut_edge( cut_face_out[3], corners_inout[7], corners_inout[4], plane );
      break;
    case MBBSPTreeBoxIter::B1265:
      plane_cut_edge( cut_face_out[0], corners_inout[0], corners_inout[1], plane );
      plane_cut_edge( cut_face_out[1], corners_inout[3], corners_inout[2], plane );
      plane_cut_edge( cut_face_out[2], corners_inout[7], corners_inout[6], plane );
      plane_cut_edge( cut_face_out[3], corners_inout[4], corners_inout[5], plane );
      break;
    case MBBSPTreeBoxIter::B2376:
      plane_cut_edge( cut_face_out[0], corners_inout[1], corners_inout[2], plane );
      plane_cut_edge( cut_face_out[1], corners_inout[0], corners_inout[3], plane );
      plane_cut_edge( cut_face_out[2], corners_inout[4], corners_inout[7], plane );
      plane_cut_edge( cut_face_out[3], corners_inout[5], corners_inout[6], plane );
      break;
    case MBBSPTreeBoxIter::B3047:
      plane_cut_edge( cut_face_out[0], corners_inout[2], corners_inout[3], plane );
      plane_cut_edge( cut_face_out[1], corners_inout[1], corners_inout[0], plane );
      plane_cut_edge( cut_face_out[2], corners_inout[5], corners_inout[4], plane );
      plane_cut_edge( cut_face_out[3], corners_inout[6], corners_inout[7], plane );
      break;
    case MBBSPTreeBoxIter::B3210:
      plane_cut_edge( cut_face_out[0], corners_inout[7], corners_inout[3], plane );
      plane_cut_edge( cut_face_out[1], corners_inout[6], corners_inout[2], plane );
      plane_cut_edge( cut_face_out[2], corners_inout[5], corners_inout[1], plane );
      plane_cut_edge( cut_face_out[3], corners_inout[4], corners_inout[0], plane );
      break;
    case MBBSPTreeBoxIter::B4567:
      plane_cut_edge( cut_face_out[0], corners_inout[0], corners_inout[4], plane );
      plane_cut_edge( cut_face_out[1], corners_inout[1], corners_inout[5], plane );
      plane_cut_edge( cut_face_out[2], corners_inout[2], corners_inout[6], plane );
      plane_cut_edge( cut_face_out[3], corners_inout[3], corners_inout[7], plane );
      break;
    default:
      return MB_FAILURE; // child is not a box
  }
  
  return MB_SUCCESS;
}

static inline
void copy_coords( double dest[3], const double source[3] )
{
  dest[0] = source[0];
  dest[1] = source[1];
  dest[2] = source[2];
}

/** reverse of plane_cut_box */
static inline
MBErrorCode plane_uncut_box( const double cut_face_in[4][3],
                             double corners_inout[8][3],
                             const MBBSPTree::Plane& plane )
{
  switch (MBBSPTreeBoxIter::side_on_plane( corners_inout, plane )) {
    case MBBSPTreeBoxIter::B0154:
      copy_coords( corners_inout[0], cut_face_in[0] );
      copy_coords( corners_inout[1], cut_face_in[1] );
      copy_coords( corners_inout[5], cut_face_in[2] );
      copy_coords( corners_inout[4], cut_face_in[3] );
      break;
    case MBBSPTreeBoxIter::B1265:
      copy_coords( corners_inout[1], cut_face_in[0] );
      copy_coords( corners_inout[2], cut_face_in[1] );
      copy_coords( corners_inout[6], cut_face_in[2] );
      copy_coords( corners_inout[5], cut_face_in[3] );
      break;
    case MBBSPTreeBoxIter::B2376:
      copy_coords( corners_inout[2], cut_face_in[0] );
      copy_coords( corners_inout[3], cut_face_in[1] );
      copy_coords( corners_inout[7], cut_face_in[2] );
      copy_coords( corners_inout[6], cut_face_in[3] );
      break;
    case MBBSPTreeBoxIter::B3047:
      copy_coords( corners_inout[3], cut_face_in[0] );
      copy_coords( corners_inout[0], cut_face_in[1] );
      copy_coords( corners_inout[4], cut_face_in[2] );
      copy_coords( corners_inout[7], cut_face_in[3] );
      break;
    case MBBSPTreeBoxIter::B3210:
      copy_coords( corners_inout[3], cut_face_in[0] );
      copy_coords( corners_inout[2], cut_face_in[1] );
      copy_coords( corners_inout[1], cut_face_in[2] );
      copy_coords( corners_inout[0], cut_face_in[3] );
      break;
    case MBBSPTreeBoxIter::B4567:
      copy_coords( corners_inout[4], cut_face_in[0] );
      copy_coords( corners_inout[5], cut_face_in[1] );
      copy_coords( corners_inout[6], cut_face_in[2] );
      copy_coords( corners_inout[7], cut_face_in[3] );
      break;
    default:
      return MB_FAILURE; // child is not a box
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBBSPTreeBoxIter::step_to_first_leaf( Direction direction )
{
  MBErrorCode rval;
  MBBSPTree::Plane plane;
  Corners clipped_corners;
  
  for (;;) {
    childVect.clear();
    rval = tool()->moab()->get_child_meshsets( mStack.back(), childVect );
    if (MB_SUCCESS != rval)
      return rval;
    if (childVect.empty()) // leaf
      break;
  
    rval = tool()->get_split_plane( mStack.back(), plane );
    if (MB_SUCCESS != rval)
      return rval;
    
    if (direction == RIGHT)
      plane.flip();
    rval = plane_cut_box( clipped_corners.coords, leafCoords, plane );
    if (MB_SUCCESS != rval)
      return rval; 
    mStack.push_back( childVect[direction] );
    stackData.push_back( clipped_corners );
  }
  return MB_SUCCESS;
}

MBErrorCode MBBSPTreeBoxIter::up()
{
  MBErrorCode rval;
  if (mStack.size() == 1)
    return MB_ENTITY_NOT_FOUND;
  
  MBEntityHandle node = mStack.back();
  Corners clipped_face = stackData.back();
  mStack.pop_back();
  stackData.pop_back();
  
  MBBSPTree::Plane plane;
  rval = tool()->get_split_plane( mStack.back(), plane );
  if (MB_SUCCESS != rval) {
    mStack.push_back( node );
    stackData.push_back( clipped_face );
    return rval;
  }
  
  rval = plane_uncut_box( clipped_face.coords, leafCoords, plane );
  if (MB_SUCCESS != rval) {
    mStack.push_back( node );
    stackData.push_back( clipped_face );
    return rval;
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBBSPTreeBoxIter::down( const MBBSPTree::Plane& plane_ref, Direction direction )
{
  childVect.clear();
  MBErrorCode rval = tool()->moab()->get_child_meshsets( mStack.back(), childVect );
  if (MB_SUCCESS != rval)
    return rval;
  if (childVect.empty())
    return MB_ENTITY_NOT_FOUND;
  
  MBBSPTree::Plane plane(plane_ref);
  if (direction == RIGHT)
    plane.flip();
  
  Corners clipped_face;
  rval = plane_cut_box( clipped_face.coords, leafCoords, plane );
  if (MB_SUCCESS != rval)
    return rval;
  
  mStack.push_back( childVect[direction] );
  stackData.push_back( clipped_face );
  return MB_SUCCESS;
}

MBErrorCode MBBSPTreeBoxIter::step( Direction direction )
{
  MBEntityHandle node, parent;
  Corners clipped_face;
  MBErrorCode rval;
  MBBSPTree::Plane plane;
  const Direction opposite = static_cast<Direction>(1-direction);
  
    // If stack is empty, then either this iterator is uninitialized
    // or we reached the end of the iteration (and return 
    // MB_ENTITY_NOT_FOUND) already.
  if (mStack.empty())
    return MB_FAILURE;
    
    // Pop the current node from the stack.
    // The stack should then contain the parent of the current node.
    // If the stack is empty after this pop, then we've reached the end.
  node = mStack.back();
  mStack.pop_back();
  clipped_face = stackData.back();
  stackData.pop_back();
  
  while(!mStack.empty()) {
      // Get data for parent entity
    parent = mStack.back();
    childVect.clear();
    rval = tool()->moab()->get_child_meshsets( parent, childVect );
    if (MB_SUCCESS != rval)
      return rval;
    rval = tool()->get_split_plane( parent, plane );
    if (MB_SUCCESS != rval)
      return rval;
    if (direction == LEFT)
      plane.flip();
    
      // If we're at the left child
    if (childVect[opposite] == node) {
        // change from box of left child to box of parent
      plane_uncut_box( clipped_face.coords, leafCoords, plane );
        // change from box of parent to box of right child
      plane.flip();
      plane_cut_box( clipped_face.coords, leafCoords, plane );
        // push right child on stack
      mStack.push_back( childVect[direction] );
      stackData.push_back( clipped_face );
        // descend to left-most leaf of the right child
      return step_to_first_leaf(opposite);
    }
    
      // The current node is the right child of the parent,
      // continue up the tree.
    assert( childVect[direction] == node );
    plane.flip();
    plane_uncut_box( clipped_face.coords, leafCoords, plane );
    node = parent;
    clipped_face = stackData.back();
    mStack.pop_back();
    stackData.pop_back();
  }
  
  return MB_ENTITY_NOT_FOUND;
}

MBErrorCode MBBSPTreeBoxIter::get_box_corners( double coords[8][3] ) const
{
  memcpy( coords, leafCoords, 24*sizeof(double) );
  return MB_SUCCESS;
}

// result = a - b
static void subtr( double result[3], const double a[3], const double b[3] )
{
  result[0] = a[0] - b[0];
  result[1] = a[1] - b[1];
  result[2] = a[2] - b[2];
}

// result = a + b + c + d
static void sum( double result[3], 
                 const double a[3], 
                 const double b[3],
                 const double c[3],
                 const double d[3] )
{
  result[0] = a[0] + b[0] + c[0] + d[0];
  result[1] = a[1] + b[1] + c[1] + d[1];
  result[2] = a[2] + b[2] + c[2] + d[2];
}

// result = a cross b
static void cross( double result[3], const double a[3], const double b[3] )
{
  result[0] = a[1]*b[2] - a[2]*b[1];
  result[1] = a[2]*b[0] - a[0]*b[2];
  result[2] = a[0]*b[1] - a[1]*b[0];
}

static double dot( const double a[3], const double b[3] )
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double MBBSPTreeBoxIter::volume() const
{
    // have planar sides, so use mid-face tripple product
  double f1[3], f2[3], f3[3], f4[3], f5[3], f6[3];
  sum( f1, leafCoords[0], leafCoords[1], leafCoords[4], leafCoords[5] );
  sum( f2, leafCoords[1], leafCoords[2], leafCoords[5], leafCoords[6] );
  sum( f3, leafCoords[2], leafCoords[3], leafCoords[6], leafCoords[7] );
  sum( f4, leafCoords[0], leafCoords[3], leafCoords[4], leafCoords[7] );
  sum( f5, leafCoords[0], leafCoords[1], leafCoords[2], leafCoords[3] );
  sum( f6, leafCoords[4], leafCoords[5], leafCoords[6], leafCoords[7] );
  double v13[3], v24[3], v65[3];
  subtr( v13, f1, f3 );
  subtr( v24, f2, f4 );
  subtr( v65, f6, f5 );
  double cr[3];
  cross( cr, v13, v24 );
  return (1./64) * dot( cr, v65 );
}

MBBSPTreeBoxIter::XSect 
MBBSPTreeBoxIter::splits( const MBBSPTree::Plane& plane ) const
{
  // test each corner relative to the plane
  unsigned result  = 0;
  for (unsigned i = 0; i < 8u; ++i) {
    double d = plane.signed_distance( leafCoords[i] );
      // if corner is on plane, than intersection 
      // will result in a degenerate hex
    if (fabs(d) < MBBSPTree::epsilon())
      return NONHEX;
      // if mark vertices above plane
    if (d > 0.0)
      result |= 1<<i;
  }
  
  switch (result) {
      // if all vertices or no vertices above plane,
      // then plane doesn't intersect
    case 0:
    case 0xFF:
      return MISS;
  
      // if there are four vertices above the plane
      // and they compose a single face of the hex,
      // then the cut will result in two hexes
    case B0154:
    case B1265:
    case B2376:
    case B3047:
    case B3210:
    case B4567:
      return SPLIT;
      
      // otherwise intersects, but split would not result
      // in two hexahedrons
    default:
      return NONHEX;
  }
}

bool MBBSPTreeBoxIter::intersects( const MBBSPTree::Plane& plane ) const
{
  // test each corner relative to the plane
  unsigned count  = 0;
  for (unsigned i = 0; i < 8u; ++i) 
    count += plane.above( leafCoords[i] );
  return count > 0 && count < 8u;
}

MBErrorCode MBBSPTreeBoxIter::sibling_side( SideBits& side_out ) const
{
  if (mStack.size() < 2) // at tree root
    return MB_ENTITY_NOT_FOUND;
  
  MBEntityHandle parent = mStack[mStack.size()-2];
  MBBSPTree::Plane plane;
  MBErrorCode rval = tool()->get_split_plane( parent, plane );
  if (MB_SUCCESS != rval)
    return MB_FAILURE;
  
  side_out = side_on_plane( leafCoords, plane );
  return MB_SUCCESS;
}

MBErrorCode MBBSPTreeBoxIter::get_neighbors( 
                      SideBits side,
                      std::vector<MBBSPTreeBoxIter>& results,
                      double epsilon ) const
{
  MBEntityHandle handle;
  MBBSPTree::Plane plane;
  MBErrorCode rval;
  int n;
   
  Corners face;
  rval = face_corners( side, leafCoords, face.coords );
  if (MB_SUCCESS != rval)
    return rval;
  
    // Move up tree until we find the split that created the specified side.
    // Push the sibling at that level onto the iterator stack as
    // all neighbors will be rooted at that node.
  MBBSPTreeBoxIter iter( *this ); // temporary iterator (don't modifiy *this)
  for (;;) {
    handle = iter.handle();
  
    rval = iter.up();
    if (MB_SUCCESS != rval) // reached root - no neighbors on that side
      return (rval == MB_ENTITY_NOT_FOUND) ? MB_SUCCESS : rval;
    
    iter.childVect.clear();
    rval = tool()->moab()->get_child_meshsets( iter.handle(), iter.childVect );
    if (MB_SUCCESS!= rval)
      return rval;
    
    rval = tool()->get_split_plane( iter.handle(), plane );
    if (MB_SUCCESS != rval)
      return rval;
    SideBits s = side_above_plane( iter.leafCoords, plane );

    if (handle == iter.childVect[0] && s == side) {
      rval = iter.down( plane, RIGHT );
      if (MB_SUCCESS != rval)
        return rval;
      break;
    }
    else if (handle == iter.childVect[1] && opposite_face(s) == side) {
      rval = iter.down( plane, LEFT );
      if (MB_SUCCESS != rval)
        return rval;
      break;
    }
  }

    // now move down tree, searching for adjacent boxes
  std::vector<MBBSPTreeBoxIter> list;
    // loop over all potential paths to neighbors (until list is empty)
  for (;;) {
      // follow a single path to a leaf, append any other potential
      // paths to neighbors to 'list'
    for (;;) { 
      rval = tool()->moab()->num_child_meshsets( iter.handle(), &n );
      if (MB_SUCCESS != rval)
        return rval;
        
        // if leaf
      if (!n) {
        results.push_back( iter );
        break; 
      }
      
      rval = tool()->get_split_plane( iter.handle(), plane );
      if (MB_SUCCESS != rval)
        return rval;
     
      bool some_above = false, some_below = false;
      for (int i = 0; i < 4; ++i) {
        double signed_d = plane.signed_distance( face.coords[i] );
        if (signed_d > -epsilon)
          some_above = true;
        if (signed_d < epsilon)
          some_below = true;
      }
     
      if (some_above && some_below) {
        list.push_back( iter );
        list.back().down( plane, RIGHT );
        iter.down( plane, LEFT );
      }
      else if (some_above) {
        iter.down( plane, RIGHT );
      }
      else if (some_below) {
        iter.down( plane, LEFT );
      }
      else {
        // tolerance issue -- epsilon to small? 2D box?
        return MB_FAILURE;
      }
    }
    
    if (list.empty())
      break;
    
    iter = list.back();
    list.pop_back();
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBBSPTreeBoxIter::calculate_polyhedron( BSPTreePoly& poly_out ) const
{
  const MBCartVect* ptr = reinterpret_cast<const MBCartVect*>(leafCoords);
  return poly_out.set( ptr );
}

MBErrorCode MBBSPTree::leaf_containing_point( MBEntityHandle tree_root,
                                              const double point[3],
                                              MBEntityHandle& leaf_out )
{
  std::vector<MBEntityHandle> children;
  Plane plane;
  MBEntityHandle node = tree_root;
  MBErrorCode rval = moab()->get_child_meshsets( node, children );
  if (MB_SUCCESS != rval)
    return rval;
  while (!children.empty()) {
    rval = get_split_plane( node, plane );
    if (MB_SUCCESS != rval)
      return rval;
      
    node = children[plane.above(point)];
    children.clear();
    rval = moab()->get_child_meshsets( node, children );
    if (MB_SUCCESS != rval)
      return rval;
  }
  leaf_out = node;
  return MB_SUCCESS;
}

MBErrorCode MBBSPTree::leaf_containing_point( MBEntityHandle root,
                                              const double point[3],
                                              MBBSPTreeIter& iter )
{
  MBErrorCode rval;
  
  rval = iter.initialize( this, root, point );
  if (MB_SUCCESS != rval)
    return rval;
  
  for (;;) {
    iter.childVect.clear();
    rval = moab()->get_child_meshsets( iter.handle(), iter.childVect );
    if (MB_SUCCESS != rval || iter.childVect.empty())
      return rval;

    Plane plane;
    rval = get_split_plane( iter.handle(), plane );
    if (MB_SUCCESS != rval)
      return rval;

    rval = iter.down( plane, (MBBSPTreeIter::Direction)(plane.above( point )) );
    if (MB_SUCCESS != rval)
      return rval;
  }
}



template <typename PlaneIter> static inline
bool ray_intersect_halfspaces( const MBCartVect& ray_pt,
                               const MBCartVect& ray_dir,
                               const PlaneIter& begin, 
                               const PlaneIter& end,
                               double& t_enter, 
                               double& t_exit )
{
  const double epsilon = 1e-12;

    // begin with inifinite ray
  t_enter = 0.0;
  t_exit  = std::numeric_limits<double>::infinity();

    // cut ray such that we keep only the portion contained
    // in each halfspace
  for (PlaneIter i = begin; i != end; ++i) {
    MBCartVect norm( i->norm );
    double coeff = i->coeff;
    double den = norm % ray_dir;
    if (fabs(den) < epsilon) { // ray is parallel to plane
      if (i->above( ray_pt.array() ))
        return false; // ray entirely outside half-space
    }
    else {
      double t_xsect = (-coeff - (norm % ray_pt)) / den;
        // keep portion of ray/segment below plane
      if (den > 0) {
        if (t_xsect < t_exit)
          t_exit = t_xsect;
      }
      else {
        if (t_xsect > t_enter)
          t_enter = t_xsect;
      }
    }
  }
  
  return t_exit >= t_enter;
}
                               
class BoxPlaneIter {
      int faceNum;
      MBBSPTree::Plane facePlanes[6];
  
  public:  
    BoxPlaneIter( const double coords[8][3] );
    BoxPlaneIter( ) : faceNum(6) {} // initialize to 'end'
    const MBBSPTree::Plane* operator->() const
      { return facePlanes + faceNum; }
    bool operator==( const BoxPlaneIter& other ) const
      { return faceNum == other.faceNum; }
    bool operator!=( const BoxPlaneIter& other ) const
      { return faceNum != other.faceNum; }
    BoxPlaneIter& operator++()
      { ++faceNum; return *this; }
};

static const int box_face_corners[6][4] = { { 0, 1, 5, 4 },
                                            { 1, 2, 6, 5 },
                                            { 2, 3, 7, 6 },
                                            { 3, 0, 4, 7 },
                                            { 3, 2, 1, 0 },
                                            { 4, 5, 6, 7 } };
 
BoxPlaneIter::BoxPlaneIter( const double coords[8][3] )
  : faceNum(0)
{
    // NOTE:  In the case of a BSP tree, all sides of the
    //        leaf will planar.
  assert( sizeof(MBCartVect) == sizeof(coords[0]) );
  const MBCartVect* corners = reinterpret_cast<const MBCartVect*>(coords);
  for (int i = 0; i < 6; ++i) {
    const int* indices = box_face_corners[i];
    MBCartVect v1 = corners[indices[1]] - corners[indices[0]];
    MBCartVect v2 = corners[indices[3]] - corners[indices[0]];
    MBCartVect n = v1 * v2;
    facePlanes[i] = MBBSPTree::Plane( n.array(), -(n % corners[indices[2]]) );
  }
}


bool MBBSPTreeBoxIter::intersect_ray( const double ray_point[3],
                                      const double ray_vect[3],
                                      double& t_enter, double& t_exit ) const
{
  BoxPlaneIter iter( this->leafCoords ), end;
  return ray_intersect_halfspaces( MBCartVect(ray_point),
                                   MBCartVect(ray_vect),
                                   iter, end,
                                   t_enter, t_exit );
}

class BSPTreePlaneIter {
    MBBSPTree* toolPtr;
    const MBEntityHandle* const pathToRoot;
    int pathPos;
    MBBSPTree::Plane tmpPlane;
    std::vector<MBEntityHandle> tmpChildren;
  public:
    BSPTreePlaneIter( MBBSPTree* tool, const MBEntityHandle* path, int path_len )
      : toolPtr(tool), pathToRoot(path), pathPos(path_len-1)
      { operator++(); }
    BSPTreePlaneIter() // initialize to 'end'
      : toolPtr(0), pathToRoot(0), pathPos(-1) {}
  
    const MBBSPTree::Plane* operator->() const
      { return &tmpPlane; }
    bool operator==( const BSPTreePlaneIter& other ) const
      { return pathPos == other.pathPos; }
    bool operator!=( const BSPTreePlaneIter& other ) const
      { return pathPos != other.pathPos; }
    BSPTreePlaneIter& operator++();
};

BSPTreePlaneIter& BSPTreePlaneIter::operator++()
{
  if (--pathPos < 0)
    return *this;

  MBEntityHandle prev = pathToRoot[pathPos+1];
  MBEntityHandle curr = pathToRoot[pathPos];
  
  MBErrorCode rval = toolPtr->get_split_plane( curr, tmpPlane );
  if (MB_SUCCESS != rval) {
    assert(false);
    pathPos = 0;
    return *this;
  }
  
  tmpChildren.clear();
  rval = toolPtr->moab()->get_child_meshsets( curr, tmpChildren );
  if (MB_SUCCESS != rval || tmpChildren.size() != 2) {
    assert(false);
    pathPos = 0;
    return *this;
  }
  
  if (tmpChildren[1] == prev) 
    tmpPlane.flip();
  return *this;
}


bool MBBSPTreeIter::intersect_ray( const double ray_point[3],
                                   const double ray_vect[3],
                                   double& t_enter, double& t_exit ) const
{
    // intersect with half-spaces defining tree
  BSPTreePlaneIter iter1( tool(), &mStack[0], mStack.size() ), end1;
  if (!ray_intersect_halfspaces( MBCartVect(ray_point),
                                 MBCartVect(ray_vect),
                                 iter1, end1,
                                 t_enter, t_exit ))
    return false;
  

    // itersect with box bounding entire tree
  double corners[8][3];
  MBErrorCode rval = tool()->get_tree_box( mStack.front(), corners );
  if (MB_SUCCESS != rval) {
    assert(false); 
    return false;
  }
  
  BoxPlaneIter iter2( corners ), end2;
  double t2_enter, t2_exit;
  if (!ray_intersect_halfspaces( MBCartVect(ray_point),
                                 MBCartVect(ray_vect),
                                 iter2, end2,
                                 t2_enter, t2_exit ))
    return false;
  
    // if intersect both box and halfspaces, check that
    // two intersections overlap
  if (t_enter < t2_enter)
    t_enter = t2_enter;
  if (t_exit > t2_exit)
    t_exit = t2_exit;
  return t_enter <= t_exit;
}
