#include <iostream>
#include <gtest/gtest.h>
#include "moab/Interface.hpp"
#include "moab/Core.hpp"

#include "DagMC.hpp"

static const char input_file[] = "test_geom.h5m";

class DagmcSimpleTest : public :: testing::Test
{ 
 protected:
  virtual void SetUp() {
  }
};

TEST_F(DagmcSimpleTest,DagmcLoadFile)
{
  moab::DagMC *dagmc = new moab::DagMC();
  moab::ErrorCode rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, moab::MB_SUCCESS);
  delete dagmc;
}

TEST_F(DagmcSimpleTest,DagmcLoadFileExternalMB)
{
  /* 1 - Test with external moab, load file in DAGMC*/
  // make new moab core
  moab::Core *mbi = new moab::Core();
  // make new dagmc into that moab
  moab::DagMC *dagmc = new moab::DagMC(mbi);

  moab::ErrorCode rval;
  // load a file
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval,moab::MB_SUCCESS);
  
  // delete dagmc
  delete dagmc;
  delete mbi;
}

TEST_F(DagmcSimpleTest, DagmcLoadFileUsingMB) 
{
  /* 2 - Test with external moab, load file in MOAB*/
  // load the file into moab rather than dagmc
  moab::ErrorCode rval;

  moab::Core *mbi = new moab::Core();
  rval = mbi->load_file(input_file);
  EXPECT_EQ(rval,moab::MB_SUCCESS);
  moab::DagMC *dagmc = new moab::DagMC(mbi);
  rval = dagmc->load_existing_contents();
  EXPECT_EQ(rval,moab::MB_SUCCESS);
  
  // delete dagmc;
  delete dagmc;
  delete mbi;
}

TEST_F(DagmcSimpleTest,DagmcLoadFileExternMoabBuildOBB)
{
  /* 1 - Test with external moab, load file in DAGMC*/
  // make new moab core
  moab::ErrorCode rval;
    
  moab::Core *mbi = new moab::Core();
  // make new dagmc into that moab
  moab::DagMC *dagmc = new moab::DagMC(mbi);

  // load a file
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval,moab::MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval,moab::MB_SUCCESS);
  // delete dagmc
  delete dagmc;
  delete mbi;
}

TEST_F(DagmcSimpleTest,DagmcLoadExistingExternMoabBuildOBB)
{
  /* 2 - Test with external moab, load file in MOAB*/
  // load the file into moab rather than dagmc
  moab::ErrorCode rval;

  moab::Core *mbi = new moab::Core();
  rval = mbi->load_file(input_file);
  EXPECT_EQ(rval,moab::MB_SUCCESS);
  moab::DagMC *dagmc = new moab::DagMC(mbi);
  rval = dagmc->load_existing_contents();
  EXPECT_EQ(rval,moab::MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval,moab::MB_SUCCESS);
  
  // delete dagmc;
  delete dagmc;
  delete mbi;
}

TEST_F(DagmcSimpleTest,DagmcLoadFileInternalMoabBuildOBB)
{
  /* 3 - Test with internal moab, load file in DAG*/
  // make new dagmc into that moab
  moab::ErrorCode rval;

  moab::DagMC *dagmc = new moab::DagMC();
  // load a file
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval,moab::MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval,moab::MB_SUCCESS);
  delete dagmc;
}

TEST_F(DagmcSimpleTest,DagmcTestOBBRetreval)
{
  // make new dagmc
  moab::DagMC *dagmc = new moab::DagMC();

  moab::ErrorCode rval;
  // load a file
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval,moab::MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval,moab::MB_SUCCESS);

  // write the file
  rval = dagmc->write_mesh("fcad",4);
  EXPECT_EQ(rval,moab::MB_SUCCESS);
  
  // now remove the dagmc instance a
  delete dagmc;

  dagmc = new moab::DagMC();
  rval = dagmc->load_file("fcad");
  EXPECT_EQ(rval,moab::MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval,moab::MB_SUCCESS);

  // delete the fcad file
  remove("fcad");
  delete dagmc;
}

TEST_F(DagmcSimpleTest,DagmcNumVols)
{
  moab::DagMC *dagmc = new moab::DagMC();
  moab::ErrorCode rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, moab::MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, moab::MB_SUCCESS);

  int expect_num_vols = 2;
  int num_vols = dagmc->num_entities(3); 
  EXPECT_EQ(expect_num_vols, num_vols);

  delete dagmc;
}

TEST_F(DagmcSimpleTest,DagmcPointInVol) 
{
  int result = 0;
  int expect_result = 1;
  int vol_idx = 1;
  
  moab::DagMC *dagmc = new moab::DagMC();

  // load a file                      
  moab::ErrorCode rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval,moab::MB_SUCCESS);   
  rval = dagmc->init_OBBTree();       
  EXPECT_EQ(rval,moab::MB_SUCCESS);   

  double xyz[3] = {0.0, 0.0, 0.0};
  moab::EntityHandle vol_h = dagmc->entity_by_index(3, vol_idx);
  rval = dagmc->point_in_volume(vol_h, xyz, result);
  EXPECT_EQ(rval, moab::MB_SUCCESS);
  EXPECT_EQ(expect_result, result);

  delete dagmc;
}

TEST_F(DagmcSimpleTest,DagmcTestObbRetrevalRayfire) 
{
  // make new dagmc  
  moab::DagMC *dagmc = new moab::DagMC();

  moab::ErrorCode rval;
  // load a file
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, moab::MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, moab::MB_SUCCESS);

  // write the file
  rval = dagmc->write_mesh("fcad",4);
  EXPECT_EQ(rval, moab::MB_SUCCESS);
  // now remove the dagmc instance a
  delete dagmc;
  
  // now create new DAGMC
  dagmc = new moab::DagMC();
  rval = dagmc->load_file("fcad");
  EXPECT_EQ(rval, moab::MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, moab::MB_SUCCESS);

  // delete the fcad file
  remove("fcad");

  // now perform full ray fire
  double eps = 1.e-6;
  int vol_idx = 1;
  // note model is cube of side 10, centred at 0,0,0, so ray fire along
  // any unit direction should be exactly 5.0
  double xyz[3] = {0.0, 0.0, 0.0};
  double dir[3] = {0.0, 0.0, 1.0};
  moab::EntityHandle next_surf;
  double next_surf_dist;
  double expect_next_surf_dist = 5.0;
  moab::EntityHandle vol_h = dagmc->entity_by_index(3, vol_idx);

  rval = dagmc->ray_fire(vol_h, xyz, dir, next_surf, next_surf_dist);
  EXPECT_EQ(rval, moab::MB_SUCCESS);
  EXPECT_EQ(expect_next_surf_dist, next_surf_dist);
  delete dagmc;
}

TEST_F(DagmcSimpleTest, DagmcRayfire)
{
  const double eps = 1e-6; // epsilon for test, faceting tol?

  moab::DagMC *dagmc = new moab::DagMC();

  moab::ErrorCode rval;                    
  // load a file                      
  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, moab::MB_SUCCESS);
  rval = dagmc->init_OBBTree(); 
  EXPECT_EQ(rval, moab::MB_SUCCESS); 

  int vol_idx = 1;
  // note model is cube of side 10, centred at 0,0,0, so ray fire along
  // any unit direction should be exactly 5.0
  double xyz[3] = {0.0, 0.0, 0.0};
  double dir[3] = {0.0, 0.0, 1.0};
  moab::EntityHandle next_surf;
  double next_surf_dist;
  double expect_next_surf_dist = 5.0;
  moab::EntityHandle vol_h = dagmc->entity_by_index(3, vol_idx);

  rval = dagmc->ray_fire(vol_h, xyz, dir, next_surf, next_surf_dist);
  EXPECT_EQ(rval, moab::MB_SUCCESS);
  EXPECT_EQ(expect_next_surf_dist, next_surf_dist);
}

TEST_F(DagmcSimpleTest,DagmcClosestTo)
{
  const double eps = 1e-6; // epsilon for test, faceting tol?

  moab::ErrorCode rval;
  // load a file                                                                                                                        
  moab::DagMC *dagmc = new moab::DagMC();

  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, moab::MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, moab::MB_SUCCESS);

  int vol_idx = 1;
  // note model is cube of side 10, centred at 0,0,0, so ray fire along
  // any unit direction should be exactly 5.0
  double xyz[3] = {-6.0, 0.0, 0.0};
  double distance; // distance from point to nearest surface
  double expect_distance = 1.0;
  moab::EntityHandle vol_h = dagmc->entity_by_index(3, vol_idx);

  rval = dagmc->closest_to_location(vol_h, xyz, distance);
  EXPECT_EQ(rval, moab::MB_SUCCESS);
  // distance should be 1.0 cm
  EXPECT_EQ(expect_distance, distance);
}

TEST_F(DagmcSimpleTest,DagmcTestBoundary)
{
  moab::ErrorCode rval;
  // load a file                                                                                                                        
  moab::DagMC *dagmc = new moab::DagMC();

  rval = dagmc->load_file(input_file);
  EXPECT_EQ(rval, moab::MB_SUCCESS);
  rval = dagmc->init_OBBTree();
  EXPECT_EQ(rval, moab::MB_SUCCESS);

  int vol_idx = 1;
  moab::EntityHandle vol_h = dagmc->entity_by_index(3, vol_idx);
  int surf_idx = 1;
  moab::EntityHandle surf_h = dagmc->entity_by_index(2, surf_idx);

  double xyz[3] = {0.0, 0.0, 5.0};
  double dir[3] = {0.0, 0.0, 1.0};
  int result;
  int expect_result = 0;

  rval = dagmc->test_volume_boundary(vol_h, surf_h, xyz, dir, result);
  EXPECT_EQ(rval, moab::MB_SUCCESS);
  // check ray leaving volume
  EXPECT_EQ(expect_result, result);
}
  
