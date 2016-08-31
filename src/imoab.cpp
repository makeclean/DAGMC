/** \file iMOAB.cpp
*/

#include "moab/Core.hpp"
using namespace moab;
#include "mpi.h"

#include "moab/imoab.h"
/*
this mhdf.h is not part of moab installation, but it is part of moab library
copy it in this folder (imoab/src/mhdf) temporarily; after imoab is part of moab, it is not neded 
*/
#include "moab/mhdf.h"
#include <stdio.h>
/*
 this is needed so far because of direct access to hdf5/mhdf
  */

#include <H5Tpublic.h>

#include <iostream>
#include "moab/ParallelComm.hpp"
#include "MBTagConventions.hpp"
#include "moab/MeshTopoUtil.hpp"
#include <sstream>

// global variables ; should they be organized in a structure, for easier references?
// or how do we keep them global?

Interface * MBI = 0;
// we should also have the default tags stored, initialized
Tag gtags[5]; // material, neumann, dirichlet,  globalID, partition tag
// should this be part of init moab?
// gtags[4]: partition tag is not yet used/initialized

struct appData {
  EntityHandle file_set;
  Range all_verts;
  Range local_verts; // it could include shared, but not owned at the interface
                     // these vertices would be all_verts if no ghosting was required
  Range ghost_vertices; // locally ghosted from other processors
  Range primary_elems;
  Range owned_elems;
  Range ghost_elems;
  int dimension; // 2 or 3, dimension of primary elements (redundant?)
  Range mat_sets;
  std::map<int, int> matIndex; // map from global block id to index in mat_sets
  Range neu_sets;
  Range diri_sets;
  std::map< std::string, Tag> tagMap;
  std::vector<Tag>  tagList;
 };

// are there reasons to have multiple moab inits? Is ref count needed?
int refCountMB( 0) ;
int iArgc;
iMOAB_String * iArgv;

/*
 list of moab entity sets corresponding to each application and pid
 */
int unused_pid =0;
// std::vector<EntityHandle>  app_FileSets; // in order of creation
std::map<std::string, int> appIdMap;     // from app string (uppercase) to app id
std::vector<ParallelComm*> pcomms; // created in order of applications, one moab::ParallelComm for each
std::vector<appData> appDatas; // the same order as pcomms

/** 
  \fn ErrorCode iMOABInitialize( int argc, iMOAB_String* argv )
  \brief Initialize the iMOAB interface implementation and create the MOAB instance, if not created already (reference counted).

  <B>Operations:</B> Collective

  \param[in] argc (int)           Number of command line arguments 
  \param[in] argv (iMOAB_String*) Command line arguments
*/
ErrCode iMOABInitialize( int argc, iMOAB_String* argv )
{
   iArgc = argc;
   iArgv = argv; // shalow copy
   if (0==refCountMB)
   {
     MBI = new Core();
     // retrieve the default tags
     const char* const shared_set_tag_names[] = {MATERIAL_SET_TAG_NAME,
                                                 NEUMANN_SET_TAG_NAME,
                                                 DIRICHLET_SET_TAG_NAME,
                                                 GLOBAL_ID_TAG_NAME};
     // blocks, visible surfaceBC(neumann), vertexBC (Dirichlet), global id, parallel partition
     for (int i = 0; i < 4; i++) {

       ErrorCode rval = MBI->tag_get_handle(shared_set_tag_names[i], 1, MB_TYPE_INTEGER,
                                           gtags[i], MB_TAG_ANY);
       if (MB_SUCCESS!=rval)
         return 1;
     }
   }
   refCountMB++;
   return MB_SUCCESS;
}

#if 0
/** 
  \fn ErrorCode iMOABInitializeFortran( )
  \brief Initialize the iMOAB interface implementation from Fortran driver and create the MOAB instance, if not created already (reference counted).

  <B>Operations:</B> Collective
*/

ErrorCode iMOABInitializeFortran( );
#endif 

/**
  \fn ErrorCode iMOABFinalize()
  \brief Finalize the iMOAB interface implementation and delete the internally reference counted MOAB instance.

  <B>Operations:</B> Collective
*/
ErrCode iMOABFinalize()
{
   refCountMB--;
   if (0==refCountMB)
      delete MBI; 
   return MB_SUCCESS;
}

/**
  \fn ErrorCode RegisterApplication( iMOAB_String app_name, MPI_Comm* comm, iMOAB_AppID pid )
  \brief Register application - Create a unique application ID and bootstrap interfaces for further queries.
  
  \note
  Internally, a mesh set will be associated with the application ID and all subsequent queries on the MOAB
  instance will be directed to this mesh/file set.

  <B>Operations:</B> Collective

  \param[in]  app_name (iMOAB_String) Application name (PROTEUS, NEK5000, etc)
  \param[in]  comm (MPI_Comm*)        MPI communicator to be used for all mesh-releated queries originating from this application
  \param[out] pid (iMOAB_AppID)       The unique pointer to the application ID
*/


ErrCode RegisterApplication( iMOAB_String app_name, MPI_Comm* comm, iMOAB_AppID pid )
{
  // will create a parallel comm for this application too, so there will be a
  // mapping from *pid to file set and to parallel comm instances
  std::string name(app_name);
  if (appIdMap.find(name)!=appIdMap.end())
  {
    std::cout << " application already registered \n";
    return 1;
  }
  *pid =  unused_pid++;
  appIdMap[name] = *pid;
  // now create ParallelComm and a file set for this application
  ParallelComm * pco = new ParallelComm(MBI, *comm);

#if 1
  int index = pco->get_id(); // t could be useful to get app id from pcomm instance ...
  assert(index==*pid);
#endif
  pcomms.push_back(pco);

  // create now the file set that will be used for loading the model in
  EntityHandle file_set;
  ErrorCode rval = MBI->create_meshset(MESHSET_SET, file_set);
  if (MB_SUCCESS != rval )
    return 1;
  appData app_data;
  app_data.file_set=file_set;
  appDatas.push_back(app_data); // it will correspond to app_FileSets[*pid] will be the file set of interest
  return 0;
}
#if 0
/**
  \fn ErrorCode RegisterFortranApplication( iMOAB_String app_name, int* comm, iMOAB_AppID pid, int app_name_length )
  \brief Register a Fortran-basedapplication - Create a unique application ID and bootstrap interfaces for further queries.
  
  \note
  Internally, the Comm object will be converted and stored as MPI_Comm. Additionally, a mesh set will be associated with the 
  application ID and all subsequent queries on the MOAB instance will be directed to this mesh/file set.

  <B>Operations:</B> Collective

  \param[in]  app_name (iMOAB_String) Application name (PROTEUS, NEK5000, etc)
  \param[in]  comm (int*)             MPI communicator to be used for all mesh-releated queries originating from this application
  \param[out] pid (iMOAB_AppID)       The unique pointer to the application ID
  \param[in]  app_name_length (int)   Length of application name string
*/


ErrorCode RegisterFortranApplication( iMOAB_String app_name, int* comm, iMOAB_AppID pid, int app_name_length );
#endif
/**
  \fn ErrorCode DeregisterApplication( iMOAB_AppID pid )
  \brief De-Register application: delete mesh (set) associated with the application ID

  <B>Operations:</B> Collective

  \param[in] pid (iMOAB_AppID) The unique pointer to the application ID
*/
ErrCode DeregisterApplication( iMOAB_AppID pid )
{
	// the file set , parallel comm are all in vectors indexed by *pid
  // assume we did not delete anything yet
  // *pid will not be reused if we register another application
  ParallelComm * pco = pcomms[*pid];
  // we could get the pco also with
  // ParallelComm * pcomm = ParallelComm::get_pcomm(MBI, *pid);
  EntityHandle fileSet = appDatas[*pid].file_set;
  // get all entities part of the file set
  Range fileents;
  ErrorCode rval = MBI->get_entities_by_handle(fileSet, fileents, /*recursive */true);
  if (MB_SUCCESS != rval )
    return 1;

  fileents.insert(fileSet);

  rval = MBI->get_entities_by_type(fileSet, MBENTITYSET, fileents); // append all mesh sets
  if (MB_SUCCESS != rval )
    return 1;
  delete pco;
  rval = MBI->delete_entities(fileents);

  if (MB_SUCCESS != rval )
    return 1;

  return 0;
}

/**
  \fn ErrorCode ReadHeaderInfo ( iMOAB_String filename, int* num_global_vertices, int* num_global_elements, int* num_dimension, int* num_parts, int filename_length )
  \brief Get global information from the file

  <B>Operations:</B> Not collective

  \param[in]  filename (iMOAB_String)    The MOAB mesh file (H5M) to probe for header information
  \param[out] num_global_vertices (int*) The total number of vertices in the mesh file
  \param[out] num_global_elements (int*) The total number of elements (of highest dimension only) 
  \param[out] num_dimension (int*)       The highest dimension of elements in the mesh (Edge=1, Tri/Quad=2, Tet/Hex/Prism/Pyramid=3)
  \param[out] num_parts (int*)           The total number of partitions available in the mesh file, typically partitioned with mbpart during pre-processing
  \param[in]  filename_length (int)      Length of the file name string
*/

ErrCode ReadHeaderInfo ( iMOAB_String filename, int* num_global_vertices, int* num_global_elements, int* num_dimension, int* num_parts, int filename_length )
{
  std::string filen(filename);
  if (filename_length< (int)filen.length())
  {
    filen = filen.substr(0,filename_length);
  }
  *num_global_vertices = 0;
  int edges = 0;
  int faces = 0;
  int regions = 0;
  *num_global_elements =0;
  *num_dimension = 0;
  *num_parts = 0;

  mhdf_FileHandle file;
  mhdf_Status status;
  unsigned long max_id;
  struct mhdf_FileDesc* data;
  /* find PARALLEL_PARTITION tag index */
  const char * pname = "PARALLEL_PARTITION";

  long int nval, junk;
  hid_t table[3];


  file = mhdf_openFile( filen.c_str() , 0, &max_id, -1, &status );
  if (mhdf_isError( &status )) {
    fprintf( stderr,"%s: %s\n", filename, mhdf_message( &status ) );
    return 1;
  }

  data = mhdf_getFileSummary( file, H5T_NATIVE_ULONG, &status );
  if (mhdf_isError( &status )) {
    fprintf( stderr,"%s: %s\n", filename, mhdf_message( &status ) );
    return 1;
  }
  *num_dimension = data->nodes.vals_per_ent;
  *num_global_vertices = (int)data->nodes.count;

  for (int i=0; i<data->num_elem_desc; i++)
  {
    struct mhdf_ElemDesc * el_desc = &(data->elems[i]);
    struct mhdf_EntDesc * ent_d = &(el_desc->desc);
    if (0==strcmp(el_desc->type, mhdf_EDGE_TYPE_NAME)) edges += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_TRI_TYPE_NAME))  faces += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_QUAD_TYPE_NAME)) faces += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_POLYGON_TYPE_NAME)) faces += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_TET_TYPE_NAME)) regions += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_PYRAMID_TYPE_NAME)) regions += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_PRISM_TYPE_NAME)) regions += ent_d->count;
    if (0==strcmp(el_desc->type, mdhf_KNIFE_TYPE_NAME)) regions += ent_d->count;
    if (0==strcmp(el_desc->type, mdhf_HEX_TYPE_NAME)) regions += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_POLYHEDRON_TYPE_NAME)) regions += ent_d->count;
    if (0==strcmp(el_desc->type, mhdf_SEPTAHEDRON_TYPE_NAME)) regions += ent_d->count;
  }
  for (int i=0; i<data->num_tag_desc; i++)
  {
    struct mhdf_TagDesc * tag_desc = &(data->tags[i]);
    if (strcmp(pname,tag_desc->name)==0)
    {
      /*printf(" tag index %d is parallel partition tag\n", i);*/
      if (tag_desc->have_sparse) {
        mhdf_openSparseTagData(file, pname, &nval, &junk, table, &status);
        if (mhdf_isError( &status )) {
          fprintf( stderr,"%s: %s\n", filename, mhdf_message( &status ) );
          return 1;
        }
      }
      else
      {
        /* could be dense tags on sets */
        mhdf_openDenseTagData(file, pname, mhdf_set_type_handle(), &nval, &status);
        if (mhdf_isError( &status )) {
          fprintf( stderr,"%s: %s\n", filename, mhdf_message( &status ) );
          return 1;
        }
      }

      *num_parts = (int)nval;
    }
  }

  // is this required?
  if (edges >0 ){
    *num_dimension = 1; // I don't think it will ever return 1
    *num_global_elements = edges;
  }
  if (faces >0 ){
    *num_dimension = 2;
    *num_global_elements = faces;
  }
  if (regions>0){
    *num_dimension = 3;
    *num_global_elements = regions;
  }
  mhdf_closeFile( file, &status );

  free( data );
  return 0;
}



/**
  \fn ErrorCode LoadMesh( iMOAB_AppID pid, iMOAB_String filename, iMOAB_String read_options, int num_ghost_layers, int filename_length, int read_options_length )
  \brief Load a MOAB mesh file in parallel and exchange ghost layers as requested

  \note
  This will exchange ghosts and the dense/sparse tags that are specified in the mesh.
  Do we need an interface to exchange tags explicitly that user specifies separately ?
  In which case, do we assume that implicit tags like GLOBAL_ID, MATERIAL_SET, NEUMANN_SET, DIRICHLET_SET are exchanged by default ?

  <B>Operations:</B> Collective

  \param[in] pid (iMOAB_AppID)            The unique pointer to the application ID
  \param[in] filename (iMOAB_String)      The MOAB mesh file (H5M) to load onto the internal application mesh set
  \param[in] read_options (iMOAB_String)  Additional options for reading the MOAB mesh file in parallel 
  \param[in] num_ghost_layers (int)       The total number of ghost layers to exchange during mesh loading
  \param[in] filename_length (int)        Length of the filename string
  \param[in] read_options_length (int)    Length of the read options string  
*/


ErrCode LoadMesh( iMOAB_AppID pid, iMOAB_String filename, iMOAB_String read_options, int * num_ghost_layers, int filename_length, int read_options_length )
{


  // make sure we use the file set and pcomm associated with the *pid
  std::ostringstream newopts;
  newopts  << read_options;
  newopts << ";PARALLEL_COMM="<<*pid;
  if (*num_ghost_layers>=1)
  {
    // if we want ghosts, we will want additional entities, the last .1
    // because the addl ents can be edges, faces that are part of the neumann sets
    newopts << ";PARALLEL_GHOSTS=3.0."<<*num_ghost_layers<<".3";
  }
  ErrorCode rval = MBI->load_file(filename, &appDatas[*pid].file_set, newopts.str().c_str());
  if (MB_SUCCESS!=rval)
    return 1;
  int rank = pcomms[*pid]->rank();
  int nprocs=pcomms[*pid]->size();

#if 1
  // some debugging stuff
  std::ostringstream outfile;
  outfile <<"TaskMesh_n" <<nprocs<<"."<< rank<<".h5m";
  // the mesh contains ghosts too, but they are not part of mat/neumann set
  // write in serial the file, to see what tags are missing
  rval = MBI->write_file(outfile.str().c_str()); // everything on root
  if (MB_SUCCESS!=rval)
    return 1;
#endif
  return 0;
}


/**
  \fn ErrorCode WriteMesh( iMOAB_AppID pid, iMOAB_String filename, iMOAB_String write_options, int filename_length, int write_options_length )
  \brief Write a MOAB mesh along with the solution tags to a file

  \note
  The interface will write one single file (H5M) and for serial files (VTK/Exodus), it will write one file per task

  <B>Operations:</B> Collective

  \param[in] pid (iMOAB_AppID)            The unique pointer to the application ID
  \param[in] filename (iMOAB_String)      The MOAB mesh file (H5M) to write all the entities contained in the internal application mesh set
  \param[in] write_options (iMOAB_String) Additional options for writing the MOAB mesh in parallel
  \param[in] filename_length (int*)       Length of the filename string
  \param[in] write_options_length (int*)  Length of the write options string
*/
ErrCode WriteMesh( iMOAB_AppID pid, iMOAB_String filename, iMOAB_String write_options, int filename_length, int write_options_length )
{
  // maybe do some processing of strings and lengths
  // maybe do some options processing?
  ErrorCode rval = MBI->write_file(filename,0, write_options,  &appDatas[*pid].file_set, 1);
  if (MB_SUCCESS!=rval)
    return 1;
  return 0;
}



/**
  \fn ErrorCode GetMeshInfo( iMOAB_AppID pid, int* num_visible_vertices, int* num_visible_elements, int *num_visible_blocks, int* num_visible_surfaceBC, int* num_visible_vertexBC )
  \brief Obtain local mesh size information based on the loaded file

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)            The unique pointer to the application ID
  \param[out] num_visible_vertices (int*)  The number of vertices in the current partition/process arranged as: owned/shared only, ghosted, total_visible (array allocated by client, <TT>size := 3</TT>)
  \param[out] num_visible_elements (int*)  The number of elements in current partition/process arranged as: owned only, ghosted, total_visible (array allocated by client, <TT>size := 3</TT>)
  \param[out] num_visible_blocks (int*)    The number of material sets in local mesh in current partition/process arranged as: owned only, ghosted/shared, total_visible (array allocated by client, <TT>size := 3</TT>)
  \param[out] num_visible_surfaceBC (int*) The number of surfaces that have a NEUMANN_SET B.C defined in local mesh in current partition/process arranged as: owned only, ghosted/shared, total_visible (array allocated by client, <TT>size := 3</TT>)
  \param[out] num_visible_vertexBC (int*)  The number of vertices that have a DIRICHLET_SET B.C defined in local mesh in current partition/process arranged as: owned only, ghosted/shared, total_visible (array allocated by client, <TT>size := 3</TT>)
*/
ErrCode GetMeshInfo( iMOAB_AppID pid, int* num_visible_vertices, int* num_visible_elements, int *num_visible_blocks, int* num_visible_surfaceBC, int* num_visible_vertexBC )
{

  // this will include ghost elements
  // we should keep a data structure with mesh, sets, etc, for each pid
  //
  appData & data = appDatas[*pid];
  EntityHandle fileSet=data.file_set;
  ErrorCode rval = MBI->get_entities_by_type(fileSet, MBVERTEX, data.all_verts, true); // recursive
  if (MB_SUCCESS!=rval)
    return 1;
  num_visible_vertices[2] = (int) data.all_verts.size();
  // we need to differentiate pure ghosted vertices from owned/shared
  // is dimension 3?
  rval = MBI->get_entities_by_dimension(fileSet, 3, data.primary_elems, true); // recursive
  if (MB_SUCCESS!=rval)
    return 1;
  data.dimension = 3;
  if (data.primary_elems.empty())
  {
    appDatas[*pid].dimension = 2;
    rval = MBI->get_entities_by_dimension(fileSet, 2, data.primary_elems, true); // recursive
    if (MB_SUCCESS!=rval)
      return 1;
    if (data.primary_elems.empty())
    {
      appDatas[*pid].dimension = 1;
      rval = MBI->get_entities_by_dimension(fileSet, 1, data.primary_elems, true); // recursive
      if (MB_SUCCESS!=rval)
        return 1;
      if (data.primary_elems.empty())
        return 1; // no elements of dimension 1 or 2 or 3
    }
      return 1;

  }
  num_visible_elements[2] = (int) data.primary_elems.size();
  // separate ghost and local/owned primary elements
  ParallelComm * pco = pcomms[*pid];

  // filter ghost vertices, from local
  rval = pco -> filter_pstatus(data.all_verts, PSTATUS_GHOST, PSTATUS_NOT, -1, &data.local_verts);
  if (MB_SUCCESS!=rval)
    return 1;
  data.ghost_vertices = subtract(data.all_verts, data.local_verts);
  num_visible_vertices[0] = (int) data.local_verts.size();
  num_visible_vertices[1] = (int) data.ghost_vertices.size();
  // get all blocks, BCs, etc

  // filter ghost elements, from local
  rval = pco -> filter_pstatus(data.primary_elems, PSTATUS_GHOST, PSTATUS_NOT, -1, &data.owned_elems);
  if (MB_SUCCESS!=rval)
    return 1;
  data.ghost_elems = subtract(data.primary_elems, data.owned_elems);
    // get all blocks, BCs, etc
  num_visible_elements[0] = (int)data.owned_elems.size();
  num_visible_elements[1] = (int)data.ghost_elems.size();

  rval = MBI->get_entities_by_type_and_tag(fileSet, MBENTITYSET, &(gtags[0]), 0, 1, data.mat_sets , Interface::UNION);
  if (MB_SUCCESS!=rval)
    return 1;
  num_visible_blocks[2] = data.mat_sets.size();
  rval = MBI->get_entities_by_type_and_tag(fileSet, MBENTITYSET, &(gtags[1]), 0, 1, data.neu_sets , Interface::UNION);
  if (MB_SUCCESS!=rval)
    return 1;
  num_visible_surfaceBC[2] = 0;
  // count how many faces are in each neu set, and how many regions are
  // adjacent to them;
  int numNeuSets = (int)data.neu_sets.size();
  for (int i=0; i<numNeuSets; i++)
  {
    Range subents;
    EntityHandle nset = data.neu_sets[i];
    rval = MBI->get_entities_by_dimension(nset, data.dimension-1, subents);
    if (MB_SUCCESS!=rval)
      return 1;
    for (Range::iterator it=subents.begin(); it!=subents.end(); ++it)
    {
      EntityHandle subent = *it;
      Range adjPrimaryEnts;
      rval = MBI->get_adjacencies(&subent, 1, data.dimension, false, adjPrimaryEnts);
      if (MB_SUCCESS!=rval)
        return 1;
      num_visible_surfaceBC[2] += (int)adjPrimaryEnts.size();
    }
  }
  rval = MBI->get_entities_by_type_and_tag(fileSet, MBENTITYSET, &(gtags[2]), 0, 1, data.diri_sets , Interface::UNION);
  if (MB_SUCCESS!=rval)
    return 1;
  num_visible_vertexBC[2]= 0;
  int numDiriSets = (int)data.diri_sets.size();
  for (int i=0; i<numDiriSets; i++)
  {
    Range verts;
    EntityHandle diset = data.diri_sets[i];
    rval = MBI->get_entities_by_dimension(diset, 0, verts);
    if (MB_SUCCESS!=rval)
      return 1;
    num_visible_vertexBC[2] += (int)verts.size();
  }


  return 0;
}


/**
  \fn ErrorCode GetVertexID( iMOAB_AppID pid, int vertices_length, iMOAB_GlobalID* global_vertex_ID, iMOAB_LocalID* local_vertex_ID )
  \brief Get the global vertex ID for all locally visible (owned and shared/ghosted) vertices

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                   The unique pointer to the application ID
  \param[in]  vertices_length (int)               The allocated size of array (typical <TT>size := num_visible_vertices</TT>)
  \param[out] global_vertex_ID (iMOAB_GlobalID*)  The global IDs for all locally visible vertices (array allocated by client)
*/
ErrCode GetVertexID( iMOAB_AppID pid, int * vertices_length, iMOAB_GlobalID* global_vertex_ID)
{
//
  Range & verts = appDatas[*pid].all_verts;
  if ((int)verts.size()!=*vertices_length)
      return 1; // problem with array length
  // global id tag is gtags[3]
  ErrorCode rval = MBI->tag_get_data(gtags[3], verts, global_vertex_ID);
  if (MB_SUCCESS!=rval)
    return 1;

  return 0;
}
/**
  \fn ErrorCode GetVertexOwnership( iMOAB_AppID pid, int vertices_length, int* visible_global_rank_ID )
  \brief Get vertex ownership information i.e., for each vertex based on the local ID, return the process that owns the vertex (local, shared or ghost)

  \note
  Do we need to implement this for owned ? That doesn't make sense.
  If we query only for shared, how do we relate the ordering ?

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)             The unique pointer to the application ID
  \param[in]  vertices_length (int)         The allocated size of array (typically <TT>size := num_visible_vertices</TT>)
  \param[out] visible_global_rank_ID (int*) The processor rank owning each of the local vertices 
*/
ErrCode GetVertexOwnership( iMOAB_AppID pid, int *vertices_length, int* visible_global_rank_ID )
{
  Range & verts = appDatas[*pid].all_verts;
  ParallelComm * pco = pcomms[*pid];
  int i=0;
  for (Range::iterator vit=verts.begin(); vit!=verts.end(); vit++, i++)
  {
    ErrorCode rval = pco->  get_owner(*vit, visible_global_rank_ID[i]);
    if (MB_SUCCESS!=rval)
      return 1;
  }
  if (i!=*vertices_length)
    return 1; // warning array allocation problem

  return 0;
}

/**
  \fn ErrorCode GetVisibleVerticesCoordinates( iMOAB_AppID pid, int coords_length, double* coordinates )
  \brief Get vertex coordinates for all local (owned and ghosted) vertices

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)     The unique pointer to the application ID
  \param[in]  coords_length (int)   The size of the allocated coordinate array (array allocated by client, <TT>size := 3*num_visible_vertices</TT>)
  \param[out] coordinates (double*) The pointer to client allocated memory that will be filled with interleaved coordinates (do need an option for blocked coordinates ?)
*/

ErrCode GetVisibleVerticesCoordinates( iMOAB_AppID pid, int * coords_length, double* coordinates )
{
  Range & verts = appDatas[*pid].all_verts;
  // interleaved coordinates, so that means deep copy anyway
  if (*coords_length!=3*(int)verts.size())
    return 1;
  ErrorCode rval = MBI->get_coords(verts, coordinates);
  if (MB_SUCCESS!=rval)
    return 1;
  return 0;
}

/**
  \fn ErrorCode GetBlockID( iMOAB_AppID pid, int block_length, iMOAB_GlobalID* global_block_IDs, iMOAB_LocalID* local_block_IDs )
  \brief Get the global vertex ID for all locally visible (owned and shared/ghosted) vertices

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                  The unique pointer to the application ID
  \param[in]  block_length (int)                 The allocated size of array (typical <TT>size := num_visible_blocks</TT>)
  \param[out] global_block_IDs (iMOAB_GlobalID*) The global IDs for all locally visible blocks (array allocated by client)
*/
ErrCode GetBlockID( iMOAB_AppID pid, int * block_length, iMOAB_GlobalID* global_block_IDs)
{
  // local id blocks? they are counted from 0 to number of visible blocks ...
  // will actually return material set tag value for global
  Range & matSets = appDatas[*pid].mat_sets;
  if (*block_length!=(int)matSets.size())
    return 1;
  // return material set tag gtags[0 is material set tag
  ErrorCode rval = MBI->tag_get_data(gtags[0], matSets, global_block_IDs);
  if (MB_SUCCESS!=rval)
    return 1;
  // populate map with index
  std::map <int, int> & matIdx = appDatas[*pid].matIndex;
  //
  for (int i=0; i<(int)matSets.size(); i++)
  {
    matIdx[global_block_IDs[i]] = i;
  }
  return 0;
}


/**
  \fn ErrorCode  GetBlockInfo(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int* vertices_per_element, int* num_elements_in_block)
  \brief Get the global block information and elements of certain type or belonging to MATERIAL_SET

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                 The unique pointer to the application ID
  \param[in]  global_block_ID (iMOAB_GlobalID)  The global block ID of the set to be queried
  \param[out] vertices_per_element (int*)       The number of vertices per element
  \param[out] num_elements_in_block (int*)      The number of elements in block
*/
ErrCode GetBlockInfo(iMOAB_AppID pid, iMOAB_GlobalID * global_block_ID,
    int* vertices_per_element, int* num_elements_in_block)
{
  std::map<int, int> & matMap = appDatas[*pid].matIndex;
  std::map<int,int>::iterator it = matMap.find(*global_block_ID);
  if (it==matMap.end())
    return 1; // error in finding block with id
  int blockIndex = matMap[*global_block_ID];
  EntityHandle matMeshSet = appDatas[*pid].mat_sets[blockIndex];
  Range blo_elems;
  ErrorCode rval = MBI-> get_entities_by_handle(matMeshSet, blo_elems);
  if (MB_SUCCESS!=rval ||  blo_elems.empty() )
    return 1;

  EntityType type = MBI->type_from_handle(blo_elems[0]);
  if (!blo_elems.all_of_type(type))
    return 1; //not all of same  type

  const EntityHandle * conn=NULL;
  int num_verts=0;
  rval = MBI->get_connectivity(blo_elems[0], conn, num_verts);
  if (MB_SUCCESS!=rval)
    return 1;
  *vertices_per_element=num_verts;
  *num_elements_in_block = (int)blo_elems.size();

  return 0;
}

/**
  \fn ErrCode  GetVisibleElementsInfo(iMOAB_AppID pid, int* num_visible_elements, iMOAB_GlobalID * element_global_IDs, int * ranks, iMOAB_GlobalID * block_IDs)
  \brief Get the elements information, global ids, ranks they belong to, block ids they belong to

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                     The unique pointer to the application ID
  \param[in]  num_visible_elements (int*)           The global block ID of the set to be queried
  \param[out] element_global_IDs (iMOAB_GlobalID*)  The number of vertices per element
  \param[out] ranks (int*)                          The owning ranks of elements
  \param[out] block_IDs (iMOAB_GlobalID*)           The block ids the elements belong
*/
ErrCode GetVisibleElementsInfo(iMOAB_AppID pid, int* num_visible_elements,
    iMOAB_GlobalID * element_global_IDs, int * ranks, iMOAB_GlobalID * block_IDs)
{
  appData & data =  appDatas[*pid];
  ParallelComm * pco = pcomms[*pid];
  ErrorCode rval = MBI-> tag_get_data(gtags[3], data.primary_elems, element_global_IDs);
  if (MB_SUCCESS!=rval)
    return 1;

  int i=0;
  for (Range::iterator eit=data.primary_elems.begin(); eit!=data.primary_elems.end(); ++eit, ++i)
  {
    rval = pco->get_owner(*eit, ranks[i]);
    if (MB_SUCCESS!=rval)
      return 1;
  }
  for (Range::iterator mit=data.mat_sets.begin(); mit!=data.mat_sets.end(); ++mit)
  {
    EntityHandle matMeshSet = *mit;
    Range elems;
    rval = MBI-> get_entities_by_handle(matMeshSet, elems);
    if (MB_SUCCESS!=rval )
      return 1;
    int valMatTag;
    rval = MBI->tag_get_data(gtags[0], &matMeshSet, 1, &valMatTag);
    if (MB_SUCCESS!=rval )
      return 1;

    for (Range::iterator eit=elems.begin(); eit!=elems.end(); ++eit)
    {
      EntityHandle eh=*eit;
      int index=data.primary_elems.index(eh);
      if (-1==index)
        return 1;
      if (-1>= *num_visible_elements)
        return 1;
      block_IDs[index]=valMatTag;
    }
  }


  return 0;
}

/** 
  \fn ErrorCode GetBlockElementConnectivities(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int connectivity_length, int* element_connectivity)
  \brief Get the connectivity for elements within a certain block, ordered based on global element IDs

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                 The unique pointer to the application ID
  \param[in]  global_block_ID (iMOAB_GlobalID)  The global block ID of the set being queried
  \param[in]  connectivity_length (int)         The allocated size of array (typical <TT>size := vertices_per_element*num_visible_elements</TT>)
  \param[out] element_connectivity (int*)       The connectivity array to store element ordering in MOAB canonical numbering scheme (array allocated by client); array contains vertex identifiers with global ID numbering
*/
ErrCode GetBlockElementConnectivities(iMOAB_AppID pid, iMOAB_GlobalID * global_block_ID, int * connectivity_length, int* element_connectivity)
{
  appData & data =  appDatas[*pid];
  std::map<int, int> & matMap = data.matIndex;
  std::map<int,int>::iterator it = matMap.find(*global_block_ID);
  if (it==matMap.end())
    return 1; // error in finding block with id
  int blockIndex = matMap[*global_block_ID];
  EntityHandle matMeshSet = data.mat_sets[blockIndex];
  std::vector<EntityHandle> elems;

  ErrorCode rval = MBI-> get_entities_by_handle(matMeshSet, elems);
  if (MB_SUCCESS!=rval ||  elems.empty() )
    return 1;


  std::vector<EntityHandle> vconnect;
  rval = MBI->get_connectivity(&elems[0], elems.size(), vconnect);
  if (MB_SUCCESS!=rval)
    return 1;
  if (*connectivity_length!=(int)vconnect.size())
    return 1; // mismatched sizes

  //gtags[3] is global id tag;
  /*rval = MBI->tag_get_data(gtags[3], &vconnect[0], connectivity_length, element_connectivity);
  if (MB_SUCCESS!=rval)
    return 1;*/
  // will return the index in data.all_verts;

  for (int i=0; i<*connectivity_length; i++)
  {
    int inx = data.all_verts.index(vconnect[i]);
    if (-1==inx)
      return 1; // error, vertex not in local range
    element_connectivity[i] = inx;
  }
  return 0;
}

/**
  \fn ErrCode GetElementConnectivity(iMOAB_AppID pid, iMOAB_LocalID * elem_index, int * connectivity_length, int* element_connectivity)
  \brief Get the connectivity for one element

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                 The unique pointer to the application ID
  \param[in]  elem_index (iMOAB_LocalID *)      Local element index
  \param[in]  connectivity_length (int)         The allocated size of array (max 27)
  \param[out] element_connectivity (int*)       The connectivity array to store connectivity in MOAB canonical numbering scheme
    array contains vertex indices in the local numbering order for vertices
*/
ErrCode GetElementConnectivity(iMOAB_AppID pid, iMOAB_LocalID * elem_index, int * connectivity_length, int* element_connectivity)
{
  appData & data =  appDatas[*pid];
  assert((*elem_index >=0)  && (*elem_index< (int)data.primary_elems.size()) );
  EntityHandle eh = data.primary_elems[*elem_index];
  int num_nodes;
  const EntityHandle * conn;
  ErrorCode rval = MBI->get_connectivity(eh, conn, num_nodes);
  if (MB_SUCCESS!=rval)
    return 1;
  if (* connectivity_length < num_nodes)
    return 1; // wrong number of vertices

  for (int i=0; i<num_nodes; i++)
  {
    int index = data.all_verts.index(conn[i]);
    if (-1==index)
      return 1;
    element_connectivity[i] = index;
  }
  * connectivity_length = num_nodes;
  return 0;
}

/**
  \fn ErrorCode GetElementOwnership(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int num_elements_in_block, int* element_ownership)
  \brief Get the element ownership within a certain block i.e., processor ID of the element owner

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                 The unique pointer to the application ID
  \param[in]  global_block_ID (iMOAB_GlobalID)  The global block ID of the set being queried
  \param[in]  num_elements_in_block (int)       The allocated size of ownership array, same as <TT>num_elements_in_block</TT> returned from GetBlockInfo()
  \param[out] element_ownership (int*)          The ownership array to store processor ID for all elements (array allocated by client) 
*/
ErrCode GetElementOwnership(iMOAB_AppID pid, iMOAB_GlobalID * global_block_ID, int * num_elements_in_block, int* element_ownership)
{
  std::map<int, int> & matMap = appDatas[*pid].matIndex;
  ParallelComm * pco = pcomms[*pid];

  std::map<int,int>::iterator it = matMap.find(*global_block_ID);
  if (it==matMap.end())
    return 1; // error in finding block with id
  int blockIndex = matMap[*global_block_ID];
  EntityHandle matMeshSet = appDatas[*pid].mat_sets[blockIndex];
  Range elems;

  ErrorCode rval = MBI-> get_entities_by_handle(matMeshSet, elems);
  if (MB_SUCCESS!=rval ||  elems.empty() )
    return 1;

  if (*num_elements_in_block!=(int)elems.size())
    return 1; // bad memory allocation
  int i=0;
  for (Range::iterator vit=elems.begin(); vit!=elems.end(); vit++, i++)
  {
    rval = pco->  get_owner(*vit, element_ownership[i]);
    if (MB_SUCCESS!=rval)
      return 1;
  }
  return 0;
}

/**
  \fn ErrorCode GetElementID(iMOAB_AppID pid, iMOAB_GlobalID global_block_ID, int num_elements_in_block, iMOAB_GlobalID* global_element_ID, iMOAB_LocalID* local_element_ID)
  \brief Get the global IDs for all locally visible elements belonging to a particular block 

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                   The unique pointer to the application ID
  \param[in]  global_block_ID (iMOAB_GlobalID)    The global block ID of the set being queried
  \param[in]  num_elements_in_block (int)         The allocated size of global element ID array, same as <TT>num_elements_in_block</TT> returned from GetBlockInfo()
  \param[out] global_element_ID (iMOAB_GlobalID*) The global IDs for all locally visible elements (array allocated by client)
  \param[out] local_element_ID (iMOAB_LocalID*)   (<I><TT>Optional</TT></I>) The local IDs for all locally visible elements (array allocated by client)
*/
ErrCode GetElementID(iMOAB_AppID pid, iMOAB_GlobalID * global_block_ID, int * num_elements_in_block, iMOAB_GlobalID* global_element_ID, iMOAB_LocalID* local_element_ID)
{
  appData & data = appDatas[*pid];
  std::map<int, int> & matMap = data.matIndex;

  std::map<int,int>::iterator it = matMap.find(*global_block_ID);
  if (it==matMap.end())
    return 1; // error in finding block with id
  int blockIndex = matMap[*global_block_ID];
  EntityHandle matMeshSet = data.mat_sets[blockIndex];
  Range elems;
  ErrorCode rval = MBI-> get_entities_by_handle(matMeshSet, elems);
  if (MB_SUCCESS!=rval ||  elems.empty() )
    return 1;



  if (*num_elements_in_block!=(int)elems.size())
    return 1; // bad memory allocation

  rval = MBI->tag_get_data(gtags[3], elems, global_element_ID);
  if (MB_SUCCESS!=rval )
    return 1;

  // check that elems are among primary_elems in data
  for (int i=0; i<*num_elements_in_block; i++)
  {
    local_element_ID[i]=data.primary_elems.index(elems[i]);
    if (-1==local_element_ID[i])
      return 1;// error, not in local primary elements
  }

  return 0;
}

/**
  \fn ErrorCode GetPointerToSurfaceBC(iMOAB_AppID pid, int surface_BC_length, iMOAB_GlobalID* global_element_ID, int* reference_surface_ID, int* boundary_condition_value)
  \brief Get the surface boundary condition information

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                   The unique pointer to the application ID
  \param[in]  surface_BC_length (int)             The allocated size of surface boundary condition array, same as <TT>num_visible_surfaceBC</TT> returned by GetMeshInfo()
  \param[out] local_element_ID (iMOAB_LocalID*)   The local element IDs that contains the side with the surface BC
  \param[out] reference_surface_ID (int*)         The surface number with the BC in the canonical reference element (e.g., 1 to 6 for HEX, 1-4 for TET)
  \param[out] boundary_condition_value (int*)     The boundary condition type as obtained from the mesh description (value of the NeumannSet defined on the element)
*/
ErrCode GetPointerToSurfaceBC(iMOAB_AppID pid, int * surface_BC_length, iMOAB_LocalID* local_element_ID, int* reference_surface_ID, int* boundary_condition_value)
{
  // we have to fill bc data for neumann sets;/

  // it was counted above, in GetMeshInfo
  appData & data = appDatas[*pid];
  int numNeuSets = (int)data.neu_sets.size();

  int index = 0; // index [0, surface_BC_length) for the arrays returned
  for (int i=0; i<numNeuSets; i++)
  {
    Range subents;
    EntityHandle nset = data.neu_sets[i];
    ErrorCode rval = MBI->get_entities_by_dimension(nset, data.dimension-1, subents);
    if (MB_SUCCESS!=rval)
      return 1;
    int neuVal ;
    rval = MBI->tag_get_data(gtags[1], &nset, 1, &neuVal);
    if (MB_SUCCESS!=rval)
      return 1;
    for (Range::iterator it=subents.begin(); it!=subents.end(); ++it)
    {
      EntityHandle subent = *it;
      Range adjPrimaryEnts;
      rval = MBI->get_adjacencies(&subent, 1, data.dimension, false, adjPrimaryEnts);
      if (MB_SUCCESS!=rval)
        return 1;
      // get global id of the primary ents, and side number of the quad/subentity
      // this is moab ordering
      for (Range::iterator pit=adjPrimaryEnts.begin(); pit!=adjPrimaryEnts.end(); pit++)
      {
        EntityHandle primaryEnt = *pit;
        // get global id
        /*int globalID;
        rval = MBI->tag_get_data(gtags[3], &primaryEnt, 1, &globalID);
        if (MB_SUCCESS!=rval)
          return 1;
        global_element_ID[index] = globalID;*/
        // get local element id
        local_element_ID[index] = data.primary_elems.index(primaryEnt);
        if (-1 == local_element_ID[index] )
          return 1; // did not find the element locally

        int side_number, sense, offset;
        rval = MBI->side_number(primaryEnt, subent,  side_number, sense, offset);
        if (MB_SUCCESS!=rval)
           return 1;
        reference_surface_ID[index] = side_number+1; // moab is from 0 to 5, it needs 1 to 6
        boundary_condition_value[index] = neuVal;
        index++;
      }
    }
  }
  if (index != *surface_BC_length)
    return 1; // error in array allocations

  return 0;
}


/**
  \fn ErrorCode GetPointerToVertexBC(iMOAB_AppID pid, int vertex_BC_length, iMOAB_GlobalID* global_vertext_ID, int* num_vertex_BC, int* boundary_condition_value)
  \brief Get the vertex boundary condition information

  <B>Operations:</B> Collective

  \param[in]  pid (iMOAB_AppID)                   The unique pointer to the application ID
  \param[in]  vertex_BC_length (int)              The allocated size of vertex boundary condition array, same as <TT>num_visible_vertexBC</TT> returned by GetMeshInfo()
  \param[out] local_vertex_ID (iMOAB_LocalID*)    The local vertex ID that has Dirichlet BC defined
  \param[out] boundary_condition_value (int*)     The boundary condition type as obtained from the mesh description (value of the DirichletSet defined on the vertex)
*/
ErrCode GetPointerToVertexBC(iMOAB_AppID pid, int * vertex_BC_length,
    iMOAB_LocalID* local_vertex_ID, int* boundary_condition_value)
{
  // it was counted above, in GetMeshInfo
  appData & data = appDatas[*pid];
  int numDiriSets = (int)data.diri_sets.size();
  int index = 0; // index [0, *vertex_BC_length) for the arrays returned
  for (int i=0; i<numDiriSets; i++)
  {
    Range verts;
    EntityHandle diset = data.diri_sets[i];
    ErrorCode rval = MBI->get_entities_by_dimension(diset, 0, verts);
    if (MB_SUCCESS!=rval)
      return 1;
    int diriVal;
    rval = MBI->tag_get_data(gtags[2], &diset, 1, &diriVal);
    if (MB_SUCCESS!=rval)
      return 1;

    for (Range::iterator vit=verts.begin(); vit!=verts.end(); ++vit)
    {
      EntityHandle vt =*vit;
      /*int vgid;
      rval = MBI->tag_get_data(gtags[3], &vt, 1, &vgid);
      if (MB_SUCCESS!=rval)
        return 1;
      global_vertext_ID[index] = vgid;*/
      local_vertex_ID[index] = data.all_verts.index(vt);
      if (-1==local_vertex_ID[index])
        return 1; // vertex was not found
      boundary_condition_value[index] = diriVal;
      index++;
    }
  }
  if (*vertex_BC_length!=index)
    return 1; // array allocation issue

  return 0;
}


/**
  \fn ErrorCode DefineTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int* tag_type, int* components_per_entity, int tag_storage_name_length)
  \brief Define a MOAB Tag corresponding to the application depending on requested types.

  \note
  In MOAB, for most solution vectors, we only need to create a "Dense", "Double" Tag

  \todo 1) Do we care about sparse/bit/integer/handle tags, and variable-length tags ?

  <B>Operations:</B> Collective

   \param[in] pid (iMOAB_AppID)               The unique pointer to the application ID
   \param[in] tag_storage_name (iMOAB_String) The tag name to store/retreive the data in MOAB
   \param[in] tag_type (int*)                 The type of MOAB tag (Dense/Sparse on Vertices/Elements, Double/Int/EntityHandle)
   \param[in] components_per_entity (int*)    The total size of vector dimension per entity for the tag (e.g., number of doubles per entity)
   \param [out] tag_index (int*)              A tag unique identifier, can be used later for tag sync
   \param[in] tag_storage_name_length (int)   The length of the tag_storage_name string
*/
ErrCode DefineTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int* tag_type, int* components_per_entity, int *tag_index,  int tag_storage_name_length)
{
  // see if the tag is already existing, and if yes, check the type, length
  if (*tag_type <0 || *tag_type>5)
    return 1; // we have 6 types of tags supported so far

  DataType tagDataType;
  TagType tagType;
  void * defaultVal = NULL;
  int * defInt = new int [*components_per_entity];
  double * defDouble = new double [*components_per_entity];
  EntityHandle * defHandle = new EntityHandle[*components_per_entity];
  for (int i=0; i<*components_per_entity; i++)
  {
    defInt[i] = 0;
    defDouble[i] = 0.;
    defHandle[i] = (EntityHandle)0;
  }
  switch (*tag_type) {
    case 0: tagDataType = MB_TYPE_INTEGER; tagType = MB_TAG_DENSE; defaultVal=defInt; break;
    case 1: tagDataType = MB_TYPE_DOUBLE;  tagType = MB_TAG_DENSE; defaultVal=defDouble; break;
    case 2: tagDataType = MB_TYPE_HANDLE;  tagType = MB_TAG_DENSE; defaultVal=defHandle; break;
    case 3: tagDataType = MB_TYPE_INTEGER; tagType = MB_TAG_SPARSE; defaultVal=defInt; break;
    case 4: tagDataType = MB_TYPE_DOUBLE;  tagType = MB_TAG_SPARSE; defaultVal=defDouble; break;
    case 5: tagDataType = MB_TYPE_HANDLE;  tagType = MB_TAG_SPARSE; defaultVal=defHandle; break;
    default : return 1; // error
  }
  std::string tag_name(tag_storage_name);
  if (tag_storage_name_length< (int)tag_name.length())
  {
    tag_name = tag_name.substr(0, tag_storage_name_length);
  }

  Tag tagHandle;
  ErrorCode rval = MBI->tag_get_handle(tag_name.c_str(), *components_per_entity,
      tagDataType,
      tagHandle, tagType, defaultVal);

  appData & data = appDatas[*pid];
  if (MB_ALREADY_ALLOCATED==rval)
  {
    std::map<std::string, Tag> & mTags = data.tagMap;
    std::map<std::string, Tag>::iterator mit = mTags.find(tag_name);
    if (mit==mTags.end())
    {
      // add it to the map
      mTags[tag_name] = tagHandle;
      // push it to the list of tags, too
      *tag_index = (int)data.tagList.size();
      data.tagList.push_back(tagHandle) ;
    }
    return 0; // OK, we found it, and we have it stored in the map tag
  }
  else if (MB_SUCCESS == rval)
  {
    data.tagMap[tag_name] = tagHandle;
    *tag_index = (int)data.tagList.size();
    data.tagList.push_back(tagHandle) ;
    return 0;
  }
  return 1; // some error, maybe the tag was not created
}

ErrCode SetIntTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name,
    int * num_tag_storage_length, int * ent_type, int* tag_storage_data,
    int tag_storage_name_length)
{
  std::string tag_name(tag_storage_name);
  if (tag_storage_name_length< (int)tag_name.length())
  {
    tag_name = tag_name.substr(0, tag_storage_name_length);
  }
  appData & data = appDatas[*pid];
  if (data.tagMap.find(tag_name)== data.tagMap.end())
    return 1; // tag not defined
  Tag tag =  data.tagMap[tag_name];

  int tagLength =0;
  ErrorCode rval = MBI->tag_get_length(tag, tagLength);
  if (MB_SUCCESS!=rval)
    return 1;
  DataType  dtype;
  rval = MBI->tag_get_data_type(tag, dtype);
  if (MB_SUCCESS!=rval || dtype!=MB_TYPE_INTEGER)
    return 1;
  // set it on a subset of entities, based on type and length
  Range * ents_to_set;
  if (* ent_type == 0)// vertices
    ents_to_set = &data.all_verts;
  else if (* ent_type == 1)
    ents_to_set = &data.primary_elems;

  int nents_to_be_set = *num_tag_storage_length /tagLength;

  if (nents_to_be_set > (int)ents_to_set->size() || nents_to_be_set<1)
    return 1; // to many entities to be set or too few
  // restrict the range; everything is contiguous; or not?

  Range contig_range( *(ents_to_set->begin()), *(ents_to_set->begin()+nents_to_be_set-1));
  rval = MBI->tag_set_data(tag, contig_range, tag_storage_data);
  if (MB_SUCCESS!=rval)
    return 1;

  return 0; // no error
}

/**
   \fn ErrCode GetIntTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, int* tag_storage_data, int tag_storage_name_length)
   \brief Store the specified values in a MOAB Tag corresponding to the application

   <B>Operations:</B> Collective

   \param[in]  pid (iMOAB_AppID)                       The unique pointer to the application ID
   \param[in]  tag_storage_name (iMOAB_String)         The tag name to store/retreive the data in MOAB
   \param[in]  num_tag_storage_length (int*)            The size of tag storage data (e.g., num_visible_vertices*components_per_entity or num_visible_elements*components_per_entity)
   \param[in]  entity_type (int*)                      type 0 for vertices, 1 for primary elements
   \param[out] tag_storage_data (int*)                 The array data of type <I>int</I> to be copied from the internal tag memory; The data is assumed to be contiguous over the local set of visible entities (either vertices or elements)
   \param[in]  tag_storage_name_length (iMOAB_String)  The length of the tag_storage_name string
*/
ErrCode GetIntTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int *num_tag_storage_length, int * ent_type, int* tag_storage_data, int tag_storage_name_length)
{
  std::string tag_name(tag_storage_name);
  if (tag_storage_name_length< (int)tag_name.length())
  {
    tag_name = tag_name.substr(0, tag_storage_name_length);
  }
  appData & data = appDatas[*pid];
  if (data.tagMap.find(tag_name)== data.tagMap.end())
    return 1; // tag not defined
  Tag tag =  data.tagMap[tag_name];

  int tagLength =0;
  ErrorCode rval = MBI->tag_get_length(tag, tagLength);
  if (MB_SUCCESS!=rval)
    return 1;
  DataType  dtype;
  rval = MBI->tag_get_data_type(tag, dtype);
  if (MB_SUCCESS!=rval || dtype!=MB_TYPE_INTEGER)
    return 1;

  // set it on a subset of entities, based on type and length
  Range * ents_to_get;
  if (* ent_type == 0)// vertices
    ents_to_get = &data.all_verts;
  else if (* ent_type == 1)
    ents_to_get = &data.primary_elems;

  int nents_to_get = *num_tag_storage_length /tagLength;

  if (nents_to_get > (int)ents_to_get->size() || nents_to_get<1)
    return 1; // to many entities to get, or too little
  // restrict the range; everything is contiguous; or not?

  Range contig_range( *(ents_to_get->begin()), *(ents_to_get->begin()+nents_to_get-1));

  rval = MBI->tag_get_data(tag, contig_range, tag_storage_data);
  if (MB_SUCCESS!=rval)
    return 1;

  return 0; // no error
}

/**
   \fn ErrCode SetDoubleTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, double* tag_storage_data, int tag_storage_name_length)
   \brief Store the specified values in a MOAB Tag corresponding to the application

   <B>Operations:</B> Collective

   \param[in]  pid (iMOAB_AppID)                       The unique pointer to the application ID
   \param[in]  tag_storage_name (iMOAB_String)         The tag name to store/retreive the data in MOAB
   \param[in]  num_tag_storage_length (int*)            The size of tag storage data (e.g., num_visible_vertices*components_per_entity or num_visible_elements*components_per_entity)
   \param[in]  entity_type (int*)                      type 0 for vertices, 1 for primary elements
   \param[out] tag_storage_data (double*)              The array data of type <I>double</I> to replace the internal tag memory; The data is assumed to be contiguous over the local set of visible entities (either vertices or elements)
   \param[in]  tag_storage_name_length (iMOAB_String)  The length of the tag_storage_name string
*/
ErrCode SetDoubleTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int * num_tag_storage_length, int * ent_type, double* tag_storage_data, int tag_storage_name_length)
{
  // exactly the same code as for int tag :) maybe should check the type of tag too
  std::string tag_name(tag_storage_name);
  if (tag_storage_name_length< (int)tag_name.length())
  {
    tag_name = tag_name.substr(0, tag_storage_name_length);
  }
  appData & data = appDatas[*pid];
  if (data.tagMap.find(tag_name)== data.tagMap.end())
    return 1; // tag not defined
  Tag tag =  data.tagMap[tag_name];

  int tagLength =0;
  ErrorCode rval = MBI->tag_get_length(tag, tagLength);
  if (MB_SUCCESS!=rval)
    return 1;

  DataType  dtype;
  rval = MBI->tag_get_data_type(tag, dtype);
  if (MB_SUCCESS!=rval || dtype!=MB_TYPE_DOUBLE)
    return 1;

  // set it on a subset of entities, based on type and length
  Range * ents_to_set;
  if (* ent_type == 0)// vertices
    ents_to_set = &data.all_verts;
  else if (* ent_type == 1)
    ents_to_set = &data.primary_elems;

  int nents_to_be_set = *num_tag_storage_length /tagLength;

  if (nents_to_be_set > (int)ents_to_set->size() || nents_to_be_set<1)
    return 1; // to many entities to be set
  // restrict the range; everything is contiguous; or not?

  Range contig_range( *(ents_to_set->begin()), *(ents_to_set->begin()+nents_to_be_set-1));

  rval = MBI->tag_set_data(tag, contig_range, tag_storage_data);
  if (MB_SUCCESS!=rval)
    return 1;

  return 0; // no error
}

/**
   \fn ErrCode GetDoubleTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int num_tag_storage_length, double* tag_storage_data, int tag_storage_name_length)
   \brief Store the specified values in a MOAB Tag corresponding to the application

   <B>Operations:</B> Collective

   \param[in]  pid (iMOAB_AppID)                The unique pointer to the application ID
   \param[in]  tag_storage_name (iMOAB_String)  The tag name to store/retreive the data in MOAB
   \param[in]  num_tag_storage_length (int)     The size of tag storage data (e.g., num_visible_vertices*components_per_entity or num_visible_elements*components_per_entity)
   \param[in]  entity_type (int*)                      type 0 for vertices, 1 for primary elements
   \param[out] tag_storage_data (double*)       The array data of type <I>double</I> to be copied from the internal tag memory; The data is assumed to be contiguous over the local set of visible entities (either vertices or elements)
   \param[in]  tag_storage_name_length (int)    The length of the tag_storage_name string
*/
ErrCode GetDoubleTagStorage(iMOAB_AppID pid, iMOAB_String tag_storage_name, int * num_tag_storage_length, int * ent_type, double* tag_storage_data, int tag_storage_name_length)
{
  // exactly the same code, except tag type check
  std::string tag_name(tag_storage_name);
  if (tag_storage_name_length< (int)tag_name.length())
  {
    tag_name = tag_name.substr(0, tag_storage_name_length);
  }
  appData & data = appDatas[*pid];
  if (data.tagMap.find(tag_name)== data.tagMap.end())
    return 1; // tag not defined
  Tag tag =  data.tagMap[tag_name];

  int tagLength =0;
  ErrorCode rval = MBI->tag_get_length(tag, tagLength);
  if (MB_SUCCESS!=rval)
    return 1;

  DataType  dtype;
  rval = MBI->tag_get_data_type(tag, dtype);
  if (MB_SUCCESS!=rval || dtype!=MB_TYPE_DOUBLE)
    return 1;

  // set it on a subset of entities, based on type and length
  Range * ents_to_get;
  if (* ent_type == 0)// vertices
    ents_to_get = &data.all_verts;
  else if (* ent_type == 1)
    ents_to_get = &data.primary_elems;

  int nents_to_get = *num_tag_storage_length /tagLength;

  if (nents_to_get > (int)ents_to_get->size() || nents_to_get<1)
    return 1; // to many entities to get
  // restrict the range; everything is contiguous; or not?

  Range contig_range( *(ents_to_get->begin()), *(ents_to_get->begin()+nents_to_get-1));
  rval = MBI->tag_get_data(tag, contig_range, tag_storage_data);
  if (MB_SUCCESS!=rval)
    return 1;

  return 0; // no error
}

/**
   \fn ErrCode SynchronizeTags(iMOAB_AppID pid,  int * num_tags, int * tag_indices, int * ent_type )
   \brief Exchange tag values for given tags

   <B>Operations:</B> Collective

   \param[in]  pid (iMOAB_AppID)                The unique pointer to the application ID
   \param[in]  num_tags (int*)                  Number of tags to exchange
   \param[in]  tag_indices (int*)               Array with tag indices of interest (size  = *num_tags)
   \param[in]  ent_type (int*)                  type of entity for tag exchange
  */
ErrCode SynchronizeTags(iMOAB_AppID pid, int * num_tag, int * tag_indices, int * ent_type)
{
  appData & data = appDatas[*pid];
  Range ent_exchange;
  std::vector<Tag> tags;
  for (int i = 0; i<* num_tag; i++)
  {
    if (tag_indices[i]<0 || tag_indices[i]>= (int)data.tagList.size())
      return 1 ; // error in tag index
    tags.push_back( data.tagList[tag_indices[i]]);
  }
  if (* ent_type==0)
    ent_exchange = data.all_verts;
  else if (*ent_type ==1 )
    ent_exchange = data.primary_elems;
  else
    return 1; // unexpected type

  ParallelComm * pco = pcomms[*pid];

  ErrorCode rval = pco->exchange_tags(tags, tags, ent_exchange);
  if (rval!=MB_SUCCESS)
    return 1;

  return 0;
}
/**
   \fn ErrCode (iMOAB_AppID pid, iMOAB_LocalID * local_index, int* num_adjacent_elements, iMOAB_LocalID* adjacent_element_IDs)
   \brief retrieve the adjacencies for the element entities

   <B>Operations:</B> Local

   \param[in]  pid (iMOAB_AppID)                      The unique pointer to the application ID
   \param[in]  local_index (iMOAB_LocalID*)           The local element ID for which adjacency information is needed
   \param[out] num_adjacent_elements (int*)           The total number of adjacent elements
   \param[out] adjacent_element_IDs (iMOAB_LocalID*)  The local element IDs of all adjacent elements to the current one (typically, num_total_sides for internal elements or num_total_sides-num_sides_on_boundary for boundary elements)
*/
ErrCode GetNeighborElements(iMOAB_AppID pid, iMOAB_LocalID * local_index, int* num_adjacent_elements, iMOAB_LocalID* adjacent_element_IDs)
{
  //; one neighbor for each subentity of dimension-1
  MeshTopoUtil mtu(MBI);
  appData & data = appDatas[*pid];
  EntityHandle eh = data.primary_elems[*local_index];
  Range adjs;
  ErrorCode rval = mtu.get_bridge_adjacencies(eh, data.dimension-1, data.dimension, adjs);
  if (rval!=MB_SUCCESS)
    return 1;
  if (* num_adjacent_elements<(int)adjs.size())
    return 1; // not dimensioned correctly
  *num_adjacent_elements=(int)adjs.size();
  for (int i=0; i<* num_adjacent_elements; i++)
  {
    adjacent_element_IDs[i] = data.primary_elems.index(adjs[i]);
  }

  return 0;
}
#if 0
/**
   \fn ErrCode GetNeighborVertices(iMOAB_AppID pid, iMOAB_GlobalID global_vertex_ID, int* num_adjacent_vertices, iMOAB_GlobalID* adjacent_vertex_IDs)
   \brief Compute the adjacencies for the vertex entities

   <B>Operations:</B> Collective

   \param[in]  pid (iMOAB_AppID)                      The unique pointer to the application ID
   \param[in]  local_vertex_ID (iMOAB_LocalID*)       The local vertex ID for which adjacency information is needed
   \param[out] num_adjacent_vertices (int*)           The total number of adjacent vertices
   \param[out] adjacent_vertex_IDs (iMOAB_LocalID*)   The local element IDs of all adjacent vertices to the current one (typically, num_total_sides for internal elements or num_total_sides-num_sides_on_boundary for boundary elements)
*/
ErrCode GetNeighborVertices(iMOAB_AppID pid, iMOAB_LocalID* local_vertex_ID, int* num_adjacent_vertices, iMOAB_LocalID* adjacent_vertex_IDs)
{
  return 0;
}
#endif
