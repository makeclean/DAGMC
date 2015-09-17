#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

// http://www-ist.cea.fr/publicea/exl-doc/201100003182.pdf

class PlasmaSourceSampler {
  public:
  PlasmaSourceSampler();
  PlasmaSourceSampler(char* filename);
  PlasmaSourceSampler(std::string filename);

  ~PlasmaSourceSampler();

  virtual int setup() =0;

  virtual int sample(double random_numbers[6], double &x, double &y, double &z,
		     double &wgt, double &energy) =0;

  void sample_energy(double random1, double random2, double ion_temp,
    double &neutron_energy);

  double dt_xs(double ion_temp);

  void convert_r_to_xy(const double r, const double random1,
		       double &x, double &y, double &wgt);
  public:
  double min_toroidal_angle;
  double max_toroidal_angle;

};

class ParametricSource : public PlasmaSourceSampler {
  public:
  ParametricSource(double ion_density_ped, double ion_density_sep,
		   double ion_density_origin, double ion_temp_ped,
		   double ion_temp_sep, double ion_temp_origin,
		   double pedistal_rad, double ion_density_peak,
		   double ion_temp_peak, double minor_radius,
		   double major_radius, double elongation,
		   double triangularity, double shafranov,
		   char* plasma_type, int plasma_id,
		   int  number_of_bins, double start_angle = 0.0,
       double end_angle = 360.0 );
  ~ParametricSource();
  int setup();

  int sample(double random_numbers[6], double &x, double &y, double &z,
	     double &wgt, double &energy);

  private:

  double ion_density(const double radius);
  double ion_temperature(const double radius);
  int sample_source_parametric( double rn_store1, double rn_store2,
				double &sampled_radius, int &sampled_bin);
  void convert_rad_to_rz(double normalised_sample_r,
       double rn_store, double &radius, double &height);

  private:
  std::vector<double> ion_kt;
  std::vector<double> source_profile;
  double bin_width;

  double ion_density_pedistal;
  double ion_density_separatrix;
  double ion_density_origin;
  double ion_temp_pedistal;
  double ion_temp_separatrix;
  double ion_temp_origin;
  double pedistal_radius;
  double ion_density_peaking;
  double ion_temp_peaking;
  double minor_radius;
  double major_radius;
  double elongation;
  double triangularity;
  double shafranov_shift;
  char* plasma_type;
  int plasma_id;
  int number_of_bins;
};

class RadialProfileSource : public PlasmaSourceSampler {

  // struct to store the read data from the rprofile data
  struct rprofile_source_line {
  double minor; // normalised minor radius value
  double ion_temp; // ion temp in MeV
  double d_dens;   // deuterium density ions/m3
  double t_dens;   // tritium density ions/m3
  };

  // struct to store the cdf & energy data
  struct rprofile_source_data {
  double probability;
  double ion_temp;
  double minor;
  double bin_width;
  };

  public:
  RadialProfileSource(double minor_rad, double major_rad, double elong, double tri,
                      double shaf, double min_angle = 0.0, double max_angle = 360.0,
                      char* filename = "r_prof.txt");
  RadialProfileSource(double minor_rad, double major_rad, double elong, double tri,
                      double shaf, double min_angle = 0.0, double max_angle = 360.0,
                      std::string filename = "r_prof.txt");

  int setup();

  int sample(double random_numbers[6], double &x, double &y, double &z,
	     double &wgt, double &energy);
  private:
  int read_file();
  int setup_rprofile_source();
  void convert_rad_to_rz(double normalised_sample_r,
                         double rn_store, double &radius, double &height);
  int sample_source_rprofile(double rn_store1, double rn_store2,
                             double &sampled_radius,int &sampled_bin);
  private:
  std::string source_file;
  std::vector<rprofile_source_line> r_source_profile;
  int num_bins;
  std::vector<rprofile_source_data> rprof_cdf_data;
  double minor_radius;
  double major_radius;
  double elongation;
  double triangularity;
  double shafranov_shift;
};

class RZProfileSource : public PlasmaSourceSampler {

  struct rz_source_line {
  double r;
  double z;
  double power;
  double ion_temp;
  };

  struct rz_source_data{
  double r;
  double r_width;
  double z;
  double z_height;
  double probability;
  double ion_temp;
  };

  public:
  RZProfileSource(double min_angle = 0.0, double max_angle = 360.0,
                  char* filename = "rz_source.txt");
  RZProfileSource(double min_angle = 0.0, double max_angle = 360.0,
                  std::string filename = "rz_source.txt");

  int setup();

  int sample(double random_numbers[6], double &x, double &y, double &z,
	     double &wgt, double &energy);

  /// private functions
  private:
  int read_file();
  int setup_rz_source();
  int sample_source_rz(double rn_store1, double rn_store2,
                         double rn_store3, double &sampled_r, double &sampled_z,
                         int &bin_value);


  /// private variables
  private:
  std::string source_file; ///< the file to read
  int num_bins; ///< number of bins in the rz file
  std::vector<rz_source_line> rz_source_profile; ///< the read in data to process
  std::vector<rz_source_data> rz_source_cdf; ///< the cdf to sample from

double minor_radius;
double major_radius;
double elongation;
double triangularity;
double shafranov_shift;
};
