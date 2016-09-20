// FluDAG/src/test/test_FlukaFuncs.cpp

#include <gtest/gtest.h>

#include "DagMC.hpp"
#include "moab/Interface.hpp"
#include "fluka_funcs.h"


#include <cmath>
#include <cassert>

moab::DagMC *DAG = new moab::DagMC();


//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class FluDAGNasaTest : public ::testing::Test
{
 protected:

  // initalize variables for each test
  virtual void SetUp() {
    // Default h5m file for testing
    std::string infile = "nasa_cancer.h5m";

    rloadval = DAG->load_file(infile.c_str());
    assert(rloadval == moab::MB_SUCCESS);

    // DAG call to initialize geometry
    rval = DAG->init_OBBTree();
    assert (rval == moab::MB_SUCCESS);

  }

 protected:

  moab::ErrorCode rloadval;
  moab::ErrorCode rval;

};


//---------------------------------------------------------------------------//
// Test setup outcomes
TEST_F(FluDAGNasaTest, SetUp)
{
  EXPECT_EQ(moab::MB_SUCCESS, rloadval);

  // DAG call to initialize geometry
  EXPECT_EQ(moab::MB_SUCCESS, rval);
  std::vector< std::string > keywords;
  rval = DAG->detect_available_props( keywords );
  EXPECT_EQ(moab::MB_SUCCESS, rval);
  rval = DAG->parse_properties( keywords );
  EXPECT_EQ(moab::MB_SUCCESS, rval);

}
//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Test test to make sure the auxscore card is correct
TEST_F(FluDAGNasaTest, FileCorrectLi6)
{
  // call the fluka function
  std::ostringstream out;
  print_auxscore(out,3,6,"TEST");
  std::string outcome = out.str();
  std::string expected_outcome = "AUXSCORE    USRTRACK   -600300              TEST\n";
  EXPECT_EQ(outcome,expected_outcome);
}
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
// Test test to make sure the auxscore card is correct
TEST_F(FluDAGNasaTest, FileCorrectLi7)
{
  // call the fluka function
  std::ostringstream out;
  print_auxscore(out,3,7,"TEST");
  std::string outcome = out.str();
  std::string expected_outcome = "AUXSCORE    USRTRACK   -700300              TEST\n";
  EXPECT_EQ(outcome,expected_outcome);
}
//---------------------------------------------------------------------------//
// Test test to make sure the auxscore card is correct
TEST_F(FluDAGNasaTest, FileCorrectU238)
{
  // call the fluka function
  std::ostringstream out;
  print_auxscore(out,92,238,"TEST");
  std::string outcome = out.str();
  std::string expected_outcome = "AUXSCORE    USRTRACK -23809200              TEST\n";
  EXPECT_EQ(outcome,expected_outcome);
}
//---------------------------------------------------------------------------//

