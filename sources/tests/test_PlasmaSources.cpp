// MCNP5/dagmc/test/test_Tally.cpp

#include "gtest/gtest.h"
#include "../plasma_source.hpp"

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class PlasmaSourceTestParametricSourceLMode : public ::testing::Test
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
      src = new ParametricSource(1.0e20, //ion_density_ped
						   1.0e19, //ion_density_sep
						   1.0e20, //ion_density_org
						   20.0,  // ion_temp_ped
						   20.0, // ion_temp_sep
						   20.0, // ion_temp_org
						   0.9, // pedistal rad
						   2.0, // ion_density_ppeaking
						   2.0, // ion_temp peaking
						   2.0, // minor rad
						   5.2, // major rad
						   1.2, // elongation
						   0.3, // triangularity
						   0.3, // shfranov shift
						   "", // name
						   1, // id
						   50); // number of bins
    }

    // deallocate memory resources
    virtual void TearDown()
    {
    }

  protected:
    ParametricSource *src;
};

//---------------------------------------------------------------------------//
TEST_F(PlasmaSourceTestParametricSourceLMode,Construct)
{
}

TEST_F(PlasmaSourceTestParametricSourceLMode,XSTest)
{
    double xs = src->dt_xs(20.0);
    EXPECT_EQ(xs,4.28240066943647578715449857714799296919001900593280e-22);
}

TEST_F(PlasmaSourceTestParametricSourceLMode,R2XY_1)
{
    double x,y,wgt;
    src->convert_r_to_xy(3.0,0.0,x,y,wgt);
    EXPECT_EQ(x,0.0);
    EXPECT_EQ(y,3.0);
    EXPECT_EQ(wgt,1.0);
}

TEST_F(PlasmaSourceTestParametricSourceLMode,R2XY_2)
{
    double x,y,wgt;
    src->convert_r_to_xy(3.0,1.0,x,y,wgt);
    EXPECT_LT(x,1.0e-7); // floating point
    EXPECT_EQ(y,3.0);
    EXPECT_EQ(wgt,1.0);
}

TEST_F(PlasmaSourceTestParametricSourceLMode,R2XY_3)
{
    double x,y,wgt;
    src->convert_r_to_xy(3.0,0.5,x,y,wgt);
    EXPECT_LT(x, 1.0e-7); // floating point
    EXPECT_EQ(y,-3.0);
    EXPECT_EQ(wgt,1.0);
}

TEST_F(PlasmaSourceTestParametricSourceLMode,R2XY_4)
{
    double x,y,wgt;
    src->convert_r_to_xy(3.0,0.25,x,y,wgt);
    EXPECT_EQ(x, 3.0); // floating point
    EXPECT_LT(y, 1.0e-7);
    EXPECT_EQ(wgt,1.0);
}

TEST_F(PlasmaSourceTestParametricSourceLMode,R2XY_5)
{
    double x,y,wgt;
    src->convert_r_to_xy(3.0,0.75,x,y,wgt);
    EXPECT_EQ(x, -3.0); // floating point
    EXPECT_LT(y, 1.0e-7);
    EXPECT_EQ(wgt,1.0);
}

TEST_F(PlasmaSourceTestParametricSourceLMode,SampleEnergy)
{
    double energy;
    src->sample_energy(1.0,1.0,0,energy);
    EXPECT_EQ(energy, 14.08); // floating point
}

TEST_F(PlasmaSourceTestParametricSourceLMode,SampleEnergy1)
{
    double energy;
    src->sample_energy(1.0,0.0,1.0,energy);
    EXPECT_EQ(energy, 14.08); // floating point
}

TEST_F(PlasmaSourceTestParametricSourceLMode,SampleEnergy2)
{
    double energy;
    src->sample_energy(0.0,1.0,10.0,energy);
    EXPECT_EQ(energy, 14.08); // floating point
}

TEST_F(PlasmaSourceTestParametricSourceLMode,SampleEnergy3)
{
    double energy;
    src->sample_energy(0,1.0,10.0,energy);
    EXPECT_EQ(energy, 14.08); // floating point
}

/////
class PlasmaSourceTestParametricSourceHMode : public ::testing::Test
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
      src = new ParametricSource(1.0e20, //ion_density_ped
               1.0e19, //ion_density_sep
               1.0e20, //ion_density_org
               20.0,  // ion_temp_ped
               20.0, // ion_temp_sep
               20.0, // ion_temp_org
               1.9, // pedistal rad
               2.0, // ion_density_ppeaking
               2.0, // ion_temp peaking
               2.0, // minor rad
               5.2, // major rad
               1.0, // elongation
               0.0, // triangularity
               0.0, // shfranov shift
               "", // name
               1, // id
               50); // number of bins
    }

    // deallocate memory resources
    virtual void TearDown()
    {
    }

  protected:
    ParametricSource *src;
};

TEST_F(PlasmaSourceTestParametricSourceHMode,Construct)
{
}
// test that we can call the dt_xs function
TEST_F(PlasmaSourceTestParametricSourceHMode,XSTest)
{
    double xs = src->dt_xs(20.0);
    EXPECT_EQ(xs,4.28240066943647578715449857714799296919001900593280e-22);
}

// test the setup function
TEST_F(PlasmaSourceTestParametricSourceHMode,SetupTest)
{
    int ec = src->setup();
    EXPECT_EQ(ec,1);
}

// test the sampler for x,y,z position
TEST_F(PlasmaSourceTestParametricSourceHMode,SampleTest1)
{
    int ec;
    double randoms[6] = {1.0,1.0,1.0,1.0,1.0,1.0};
    double x,y,z,erg,wgt;
    ec = src->setup();
    ec = src->sample(randoms,x,y,z,wgt,erg);
    EXPECT_EQ(ec,1);
    // for rand[0] = 1
    // for rand[1] = 1
    // selects the radial section
    // rand[2] sets rz should have r be major+minor
    // z should be 0.0
    EXPECT_EQ(std::sqrt((x*x)+(y*y)),5.2+2.0);
    // expect basically 0
    EXPECT_LT(z,1.0e-7);
    EXPECT_EQ(wgt,1.0);
}

// test the sampler for x,y,z position
TEST_F(PlasmaSourceTestParametricSourceHMode,SampleTest2)
{
    int ec;
    double randoms[6] = {1.0,1.0,0.25,1.0,1.0,1.0};
    double x,y,z,erg,wgt;
    ec = src->setup();
    ec = src->sample(randoms,x,y,z,wgt,erg);
    EXPECT_EQ(ec,1);
    // for rand[0] = 1
    // for rand[1] = 1
    // selects the radial section
    // rand[2] sets rz should have r be major+minor
    // z should be 0.0
    EXPECT_FLOAT_EQ(std::sqrt((x*x)+(y*y)),5.2);
    // expect basically 0
    EXPECT_DOUBLE_EQ(z,2.0);
    EXPECT_EQ(wgt,1.0);
}

// test the sampler for x,y,z position
TEST_F(PlasmaSourceTestParametricSourceHMode,SampleTest3)
{
    int ec;
    double randoms[6] = {1.0,1.0,0.75,1.0,1.0,1.0};
    double x,y,z,erg,wgt;
    ec = src->setup();
    ec = src->sample(randoms,x,y,z,wgt,erg);
    EXPECT_EQ(ec,1);
    // for rand[0] = 1
    // for rand[1] = 1
    // selects the radial section
    // rand[2] sets rz should have r be major+minor
    // z should be 0.0
    EXPECT_FLOAT_EQ(std::sqrt((x*x)+(y*y)),5.2);
    // expect basically 0
    EXPECT_DOUBLE_EQ(z,-2.0);
    EXPECT_EQ(wgt,1.0);
}

// test the sampler for x,y,z position
TEST_F(PlasmaSourceTestParametricSourceHMode,SampleTest4)
{
    int ec;
    double randoms[6] = {1.0,1.0,0.5,1.0,1.0,1.0};
    double x,y,z,erg,wgt;
    ec = src->setup();
    ec = src->sample(randoms,x,y,z,wgt,erg);
    EXPECT_EQ(ec,1);
    // for rand[0] = 1
    // for rand[1] = 1
    // selects the radial section
    // rand[2] sets rz should have r be major+minor
    EXPECT_NEAR(std::sqrt((x*x)+(y*y)),3.2,1e-8);
    EXPECT_NEAR(z,0.0,1.0e-8);
    // z should be 0.0
    EXPECT_EQ(wgt,1.0);
}

/////
class PlasmaSourceRZMode : public ::testing::Test
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
    }

    // deallocate memory resources
    virtual void TearDown()
    {
    }
};

// instanciate new RZ constructor
TEST_F(PlasmaSourceRZMode,Construct)
{
    RZProfileSource *src = new RZProfileSource(0,360.0,std::string(""));
}

// instanciate setup no file
TEST_F(PlasmaSourceRZMode,SetupNoFile)
{
    RZProfileSource *src = new RZProfileSource(0,360.0,std::string(""));
    int ec = src->setup();
    // expect failure
    EXPECT_EQ(ec,0);
}

// instanciate setup right file
TEST_F(PlasmaSourceRZMode,SetupWithFile)
{
    RZProfileSource *src = new RZProfileSource(0,360.0,std::string("rzsource.txt"));
    int ec = src->setup();
    // expect success
    EXPECT_EQ(ec,1);
}

// instanciate setup right file
TEST_F(PlasmaSourceRZMode,SampleWithFile)
{
    RZProfileSource *src = new RZProfileSource(0,360.0,std::string("rzsource.txt"));
    int ec;
    ec = src->setup();
    // expect success
    EXPECT_EQ(ec,1);
    double randoms[6] = {1.0,1.0,1.0,1.0,1.0,1.0};
    double x,y,z,wgt,erg;
    ec = src->sample(randoms,x,y,z,wgt,erg);
    EXPECT_EQ(wgt,1.0);
    // bin width is 0.1, therefore +/- 0.05 
    EXPECT_NEAR(x,0.0,1.0e-7);
    EXPECT_NEAR(y,12.15,1.0e-7);
    EXPECT_NEAR(z,0.85,1.0e-7);
    EXPECT_NEAR(erg,14.08,1.0e-7);
    // expect success
    EXPECT_EQ(ec,1);
}

// instanciate setup right file
TEST_F(PlasmaSourceRZMode,SampleWithFile2)
{
    RZProfileSource *src = new RZProfileSource(0,360.0,std::string("rzsource.txt"));
    int ec;
    ec = src->setup();
    // expect success
    EXPECT_EQ(ec,1);
    double randoms[6] = {1.0,1.0,1.0,1.0,1.0,1.0};
    double x,y,z,wgt,erg;
    ec = src->sample(randoms,x,y,z,wgt,erg);
    EXPECT_EQ(wgt,1.0);
    // bin width is 0.1, therefore +/- 0.05 
    EXPECT_NEAR(x,0.0,1.0e-7);
    EXPECT_NEAR(y,12.15,1.0e-7);
    EXPECT_NEAR(z,0.85,1.0e-7);
    EXPECT_NEAR(erg,14.08,1.0e-7);
    // expect success
    EXPECT_EQ(ec,1);
}

/////
class PlasmaSourceRProfMode : public ::testing::Test
{
  protected:
    // initialize variables for each test
    virtual void SetUp()
    {
    }

    // deallocate memory resources
    virtual void TearDown()
    {
    }
};

// instanciate new Profile constructor 
TEST_F(PlasmaSourceRProfMode,Construct)
{
  RadialProfileSource *src = new RadialProfileSource(2.0, // minor radius
						     5.2, // major radius
						     1.0, // elongation
						     0.0, // triangularity
						     0.0, // shafranov shift
						     0.0, // start angle
						     360.0, // end angle
						     std::string("")); // filenaame
}

// instanciate setup no file
TEST_F( PlasmaSourceRProfMode,SetupNoFile)
{
  RadialProfileSource *src = new RadialProfileSource(2.0, // minor radius
						     5.2, // major radius
						     1.0, // elongation
						     0.0, // triangularity
						     0.0, // shafranov shift
						     0.0, // start angle
						     360.0, // end angle
						     std::string("")); // filenaame
    int ec = src->setup();
    // expect failure
    EXPECT_EQ(ec,0);
}

// instanciate setup wrong file
TEST_F( PlasmaSourceRProfMode,SetupWrongFile)
{
  RadialProfileSource *src = new RadialProfileSource(2.0, // minor radius
					 5.2, // major radius
					 1.0, // elongation
					 0.0, // triangularity
					 0.0, // shafranov shift
					 0.0, // start angle
					 360.0, // end angle
					 std::string("this_file_doesnt_exist.txt")); // filenaame
    int ec = src->setup();
    // expect failure
    EXPECT_EQ(ec,0);
}


// instanciate setup no file
TEST_F( PlasmaSourceRProfMode,SetupWithFile)
{
  RadialProfileSource *src = new RadialProfileSource(2.0, // minor radius
					 5.2, // major radius
					 1.0, // elongation
					 0.0, // triangularity
					 0.0, // shafranov shift
					 0.0, // start angle
					 360.0, // end angle
					 std::string("r_prof.txt")); // filenaame
    int ec = src->setup();
    // expect success
    EXPECT_EQ(ec,1);
}

// instanciate setup no file
TEST_F( PlasmaSourceRProfMode,SetupWithFileSample)
{
  RadialProfileSource *src = new RadialProfileSource(2.0, // minor radius
						     5.2, // major radius
						     1.0, // elongation
						     0.0, // triangularity
						     0.0, // shafranov shift
						     0.0, // start angle
						     360.0, // end angle
						     std::string("r_prof.txt")); // filenaame
    int ec = src->setup();
    // expect success    
    EXPECT_EQ(ec,1);
    double randoms[6] = {1.0,1.0,1.0,1.0,1.0,1.0};
    double x,y,z,wgt,erg;
    ec = src->sample(randoms,x,y,z,wgt,erg);
    // expect success
    EXPECT_EQ(ec,1);
    // wgt should be 1.0
    EXPECT_EQ(wgt,1.);

    // for rand[0] = 1
    // for rand[1] = 1
    // selects the radial section
    // rand[2] sets rz should have r be major+minor
    // z should be 0.0
    EXPECT_NEAR(std::sqrt((x*x)+(y*y)),7.2,1.0e-7);
    // expect basically 0
    EXPECT_NEAR(z,0.0,1.0e-7);
}

// instanciate setup no file, sample again
TEST_F( PlasmaSourceRProfMode,SetupWithFileSample1)
{
  RadialProfileSource *src = new RadialProfileSource(2.0, // minor radius
						     5.2, // major radius
						     1.0, // elongation
						     0.0, // triangularity
						     0.0, // shafranov shift
						     0.0, // start angle
						     360.0, // end angle
						     std::string("r_prof.txt")); // filenaame
    int ec = src->setup();
    // expect success    
    EXPECT_EQ(ec,1);
    double randoms[6] = {1.0,1.0,0.25,1.0,1.0,1.0};
    double x,y,z,wgt,erg;
    ec = src->sample(randoms,x,y,z,wgt,erg);
    // expect success
    EXPECT_EQ(ec,1);
    // wgt should be 1.0
    EXPECT_EQ(wgt,1.0);

    // for rand[0] = 1
    // for rand[1] = 1
    // selects the radial section
    // rand[2] sets rz should have r be major
    // z should be minor
    EXPECT_NEAR(std::sqrt((x*x)+(y*y)),5.2,1.0e-7);
    // expect basically 2
    EXPECT_NEAR(z,2.0,1.0e-7);
}

// instanciate setup no file, sample again
TEST_F( PlasmaSourceRProfMode,SetupWithFileSample2)
{
  RadialProfileSource *src = new RadialProfileSource(2.0, // minor radius
						     5.2, // major radius
						     1.0, // elongation
						     0.0, // triangularity
						     0.0, // shafranov shift
						     0.0, // start angle
						     360.0, // end angle
						     std::string("r_prof.txt")); // filenaame
    int ec = src->setup();
    // expect success    
    EXPECT_EQ(ec,1);
    double randoms[6] = {1.0,1.0,0.75,1.0,1.0,1.0};
    double x,y,z,wgt,erg;
    ec = src->sample(randoms,x,y,z,wgt,erg);
    // expect success
    EXPECT_EQ(ec,1);
    // wgt should be 1.0
    EXPECT_EQ(wgt,1.0);

    // selects the radial section
    // rand[2] sets rz should have r be major
    // z should be -minor
    EXPECT_NEAR(std::sqrt((x*x)+(y*y)),5.2,1.0e-7);
    // expect basically 0
    EXPECT_NEAR(z,-2.0,1.0e-7);
}

// instanciate setup no file, sample again
TEST_F( PlasmaSourceRProfMode,SetupWithFileSample3)
{
  RadialProfileSource *src = new RadialProfileSource(2.0, // minor radius
						     5.2, // major radius
						     1.0, // elongation
						     0.0, // triangularity
						     0.0, // shafranov shift
						     0.0, // start angle
						     360.0, // end angle
						     std::string("r_prof.txt")); // filenaame
  int ec = src->setup();
  // expect success    
  EXPECT_EQ(ec,1);
  double randoms[6] = {1.0,1.0,0.5,1.0,1.0,1.0};
  double x,y,z,wgt,erg;
  ec = src->sample(randoms,x,y,z,wgt,erg);
  // expect success
  EXPECT_EQ(ec,1);
  // wgt should be 1.0
  EXPECT_EQ(wgt,1.0);
  
  // selects the radial section
  // rand[2] sets rz should have r be major-minor
  // z should be 0.0
  EXPECT_NEAR(std::sqrt((x*x)+(y*y)),3.2,1.0e-7);
  // expect basically 0
  EXPECT_NEAR(z,0.0,1.0e-7);
  std::cout << erg << std::endl;
}
