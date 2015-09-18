#include "plasma_source.hpp"
#include <iostream>


#define TWOPI 6.28318530718
#define RADS_PER_DEG 0.01745329251

// Generic Constructor
PlasmaSourceSampler::PlasmaSourceSampler() {
}

// Generic Destructor
PlasmaSourceSampler::~PlasmaSourceSampler() {
}


//Given an Ion Temp in KeV return the DT xs value in barns
double PlasmaSourceSampler::dt_xs(double ion_temp) {
  double dt;
  double c[7] = {2.5663271e-18,19.983026,2.5077133e-2,
                 2.5773408e-3,6.1880463e-5,6.6024089e-2,
                 8.1215505e-3};

  double u = 1.0-ion_temp*(c[2]+ion_temp*(c[3]-c[4]*ion_temp))
    /(1.0+ion_temp*(c[5]+c[6]*ion_temp));

  dt = c[0]/(std::pow(u,5./6.)*std::pow(ion_temp,2./3.0));
  dt *= exp(-1.*c[1]*std::pow(u/ion_temp,1./3.));

  return dt;
}

// turns the radial coordinate into a random xy
void PlasmaSourceSampler::convert_r_to_xy(const double r, const double random,
					     double &x, double &y, double &wgt) {

  double toroidal_extent = this->max_toroidal_angle - this->min_toroidal_angle;
  double toroidal_angle = (toroidal_extent*random) + this->min_toroidal_angle;
  x = r*sin(toroidal_angle*RADS_PER_DEG);
  y = r*cos(toroidal_angle*RADS_PER_DEG);
  wgt = toroidal_extent/360.;
  return;
}

// sample the energy of the neutrons, updates energy neutron in mev
void PlasmaSourceSampler::sample_energy(double random_number1, double
  random_number2, double ion_temp, double &energy_neutron) {
  // generate the normally distributed number
  double sample1;

  if (random_number1 > 0.0 )
    sample1 = std::sqrt(-2.0*std::log(random_number1));
  else
    sample1 = 0.0;

  double sample2 = cos(TWOPI*random_number2);
  energy_neutron = (5.59/2.35)*ion_temp*sample1*sample2;
  energy_neutron += 14.08;
  return;
}


ParametricSource::ParametricSource(double ion_density_ped, double ion_density_sep,
				   double ion_density_org, double ion_temp_ped,
				   double ion_temp_sep, double ion_temp_org,
				   double pedistal_rad, double ion_density_peak,
				   double ion_temp_peak, double minor_rad,
				   double major_rad, double elong,
				   double tri, double shaf,
				   char* p_type, int p_id,
				   int  num, double start_angle,
           double end_angle) {

  ion_density_pedistal = ion_density_ped;
  ion_density_separatrix = ion_density_sep;
  ion_density_origin = ion_density_org;
  ion_temp_pedistal = ion_temp_ped;
  ion_temp_separatrix = ion_temp_sep;
  ion_temp_origin = ion_temp_org;
  pedistal_radius = pedistal_rad;
  ion_density_peaking = ion_density_peak;
  ion_temp_peaking = ion_temp_peak;
  minor_radius = minor_rad;
  major_radius = major_rad;
  elongation = elong;
  triangularity = tri;
  shafranov_shift = shaf;
  plasma_type = p_type;
  plasma_id = p_id;
  number_of_bins = num;
  min_toroidal_angle = start_angle;
  max_toroidal_angle = end_angle;
}

ParametricSource::~ParametricSource()
{
}

double ParametricSource::ion_density(const double radius)
{
  double ion_dens = 0.0;

  if( plasma_id == 0 ) {
    ion_dens = ion_density_origin;
    if ( radius > 0 )
      ion_dens *= (1.0-std::pow(radius/minor_radius,2));
  } else {
    if(radius <= pedistal_radius) {
      ion_dens += ion_density_pedistal;
      double product;
      product = 1.0 - std::pow(radius/pedistal_radius,2);
      product = std::pow(product,ion_density_peaking);
      ion_dens += (ion_density_origin-ion_density_pedistal)*(product);
    } else {
      ion_dens += ion_density_separatrix;
      double product;
      product = ion_density_pedistal-ion_density_separatrix;
      ion_dens += product*(minor_radius-radius)/
        (minor_radius-pedistal_radius);
    }
  }

  return ion_dens;
}

double ParametricSource::ion_temperature(const double radius)
{
  double ion_temp = 0.0;
  if( plasma_id == 0 ) {
    ion_temp = ion_temp_origin;
    ion_temp *= (1.0-std::pow(radius/minor_radius,ion_temp_peaking));
  } else {
    if(radius <= pedistal_radius) {
      ion_temp += ion_temp_pedistal;
      double product;
      if ( radius > 0.0 )
          product = 1.0 - std::pow(radius/pedistal_radius,2);
      else
          product = 1.0;

      product = std::pow(product,ion_temp_peaking);

      ion_temp += (ion_temp_origin-
                   ion_temp_pedistal)*(product);

    } else {
      ion_temp += ion_temp_separatrix;
      double product;
      product = ion_temp_pedistal -
       ion_temp_separatrix;
      ion_temp += product*(minor_radius-radius)/
        (minor_radius-pedistal_radius);
    }
  }

  return ion_temp;
}


int ParametricSource::setup() {
  double ion_d; // ion density
  double ion_t; // ion temp
  double sig_dt; // dt xs

  std::vector<double> src_strength; // the source strength, n/m3
  double r;

  bin_width = minor_radius/float(number_of_bins);
  double total = 0.0; // total source strength

  for (int i = 0 ; i < number_of_bins ; i++) {
    r = bin_width * float(i);
    ion_d = ion_density(r);
    ion_t = ion_temperature(r);
    src_strength.push_back(std::pow(ion_d,2)*dt_xs(ion_t));
    ion_kt.push_back(sqrt(ion_t/1000.0)); // convert to sqrt(MeV)
    total += src_strength[i];
  }


  // normalise the source profile
  double sum = 0.0 ;
  for ( int i = 0 ; i < number_of_bins ; i++) {
    sum += src_strength[i];
    source_profile.push_back(sum/total);
  }
  return 1;
}

/*
 * sample the pdf src_profile, to generate the sampled minor radius
 */
int ParametricSource::sample_source_parametric( double rn_store1,
  double rn_store2,double &sampled_radius, int &sampled_bin) {
  for ( int i = 0 ; i < number_of_bins ; i++ ) {
    if ( rn_store1 <= source_profile[i] ) {
      if ( i > 1 ) {
        sampled_radius = (float(i)*(bin_width)) + (bin_width*rn_store2);
        sampled_bin = i;
        return 1;
      } else {
        sampled_radius = bin_width*rn_store2;
        sampled_bin = i;
        return 1;
      }
    }
  }

  std::cerr << "error" << std::endl;
  std::cerr << "Sample position greater than plasma radius" << std::endl;

  return 0;
}

//  convert the sampled radius to an rz coordinate by using plasma parameters
void ParametricSource::convert_rad_to_rz(double normalised_sample_r,
  double rn_store, double &radius, double &height) {

  double alpha = TWOPI*rn_store;
  //  std::cout << alpha << std::endl;
  double shift = shafranov_shift*(1.0-std::pow((normalised_sample_r/minor_radius),2));
  //  std::cerr << std::pow(*minor_sampled/(*minor_radius),2) << std::endl;

  radius = major_radius + shift + normalised_sample_r*cos(alpha+(triangularity*sin(alpha)));
  height = elongation*normalised_sample_r*sin(alpha);

  //  std::cout << "r=" << radius << " h=" << height << std::endl;
};

int ParametricSource::sample(double random_numbers[6],
  double &x, double &y, double &z, double &wgt, double &energy) {

  double normalised_rad;
  int bin;

  int ec = this->sample_source_parametric(random_numbers[0],random_numbers[1],
    normalised_rad,bin);
  if(ec != 1 ) return 0;
  double radius,height;

  // convert normalised radius to r,z
  this->convert_rad_to_rz(normalised_rad,random_numbers[2],radius,z);

  // convert radius to a random x,y
  this->convert_r_to_xy(radius,random_numbers[3],x,y,wgt);
  // sample the neutron energy
  this->sample_energy(random_numbers[4],random_numbers[5],ion_kt[bin],energy);

  return 1;
}

RZProfileSource::RZProfileSource(double min_angle,
                                 double max_angle, char* filename ) {
  source_file = std::string(filename);
  min_toroidal_angle = min_angle;
  max_toroidal_angle = max_angle;
}

RZProfileSource::RZProfileSource(double min_angle,
                                 double max_angle, std::string filename ) {
    source_file = filename;
    min_toroidal_angle = min_angle;
    max_toroidal_angle = max_angle;
}

int RZProfileSource::setup() {
  int ec = read_file();
  if(ec != 1) return ec;

  ec = setup_rz_source();
  if(ec != 1) return ec;

  return 1;
}

int RZProfileSource::sample(double random_numbers[6], double &x, double &y,
  double &z, double &wgt, double &energy) {
  int bin_sampled;
  double r;
  int ec = sample_source_rz(random_numbers[0],random_numbers[1],random_numbers[2],
    r,z,bin_sampled);

  if(ec != 1) return ec;

  sample_energy(random_numbers[3],random_numbers[4],
                rz_source_cdf[bin_sampled].ion_temp,energy);

  convert_r_to_xy(r,random_numbers[5],x,y,wgt);
  return 1;
}

int RZProfileSource::sample_source_rz(double rn_store1, double rn_store2,
                       double rn_store3, double &sampled_r, double &sampled_z,
                       int &bin_value) {

  for ( int i = 0 ; i < num_bins ; i++ ) {
    if ( rn_store1 <= rz_source_cdf[i].probability ) {
      sampled_r = rz_source_cdf[i].r + ((rn_store2-0.5)*rz_source_cdf[i].r_width);
      sampled_z = rz_source_cdf[i].z + ((rn_store3-0.5)*rz_source_cdf[i].z_height);
      bin_value = i;
      return 1;
    }
  }
  // if we arehere, then somehow we have failed to sample. Stop
  return 0;
}


int RZProfileSource::read_file()
{
  // open up the file
  //char* filename = souce_file.c_str();
  std::ifstream rz_input(source_file);
  std::string line;

  std::string buf;

  if (rz_input.is_open()) {
    while ( getline(rz_input, line) )
      {
        std::vector<double> tokens;
        std::stringstream ss(line);

        while ( ss >> buf)
          tokens.push_back(atof(buf.c_str()));
        rz_source_line data;
        data.r = tokens[0];
        data.z = tokens[1];
        data.power = tokens[2];
        data.ion_temp = tokens[3]/1.0e6; // convert to MeV
        rz_source_profile.push_back(data);
      }

    // now data is read in
    num_bins = rz_source_profile.size();
    rz_input.close();
   } else {
    std::cout << "Failed to open '"<< source_file <<"' is it in the right place?" << std::endl;
    return 0;
  }
  return 1;
}

int RZProfileSource::setup_rz_source(){

  std::vector<rz_source_line>::iterator it;

  // step through data to figure out probabilities
  double total_prob = 0.0; // total prob
  std::vector<double> pdf; // the pdf
  rz_source_line src_line;


  for ( it = rz_source_profile.begin() ; it != rz_source_profile.end() ; ++it ) {
    src_line = *it;
    total_prob += src_line.power;
    pdf.push_back(total_prob); // collect the pdf
  }

  // calculate the bin width, assume symetric and uniform
  double bin_height = rz_source_profile[num_bins-1].z-rz_source_profile[num_bins-2].z;
  double bin_width = bin_height;

  // make the cdf
  rz_source_data data; // actual data stored for use later
  for ( int i = 0 ; i <= num_bins ; i++ ) {
    data.r = rz_source_profile[i].r;
    data.r_width = bin_width;
    data.z = rz_source_profile[i].z;
    data.z_height = bin_height;
    data.probability = pdf[i]/total_prob;
    data.ion_temp = rz_source_profile[i].ion_temp;
    rz_source_cdf.push_back(data);
  }

  return 1;
}

RadialProfileSource::RadialProfileSource(double minor_rad,
                                        double major_rad, double elong, double tri,
                                        double shaf, double min_angle,
                                        double max_angle,char* filename ) {
  source_file = std::string(filename);
  minor_radius = minor_rad;
  major_radius = major_rad;
  elongation = elong;
  triangularity = tri;
  shafranov_shift = shaf;
  min_toroidal_angle = min_angle;
  max_toroidal_angle = max_angle;
}

RadialProfileSource::RadialProfileSource(double minor_rad,
                                        double major_rad, double elong, double tri,
                                        double shaf, double min_angle,
                                        double max_angle, std::string filename ) {
    source_file = filename;
    minor_radius = minor_rad;
    major_radius = major_rad;
    elongation = elong;
    triangularity = tri;
    shafranov_shift = shaf;
    min_toroidal_angle = min_angle;
    max_toroidal_angle = max_angle;
}

int RadialProfileSource::setup() {
  int ec;
  ec = read_file();
  if(ec != 1) return ec;
  ec = setup_rprofile_source();
  if(ec != 1) return ec;
}

int RadialProfileSource::sample(double random_numbers[6], double &x, double &y,
  double &z, double &wgt, double &energy) {
  int bin_sampled;
  double r;
  double normalised_radial;
  int ec = sample_source_rprofile(random_numbers[0],random_numbers[1],
                                  normalised_radial,bin_sampled);

  if(ec != 1) return ec;
  // normalised radial cannot be greater than 1.0
  if(normalised_radial > 1.0) return 0;
  convert_rad_to_rz(normalised_radial, random_numbers[2], r ,z);

  //  std::cout << normalised_radial << " " << r << " " << z << std::endl;
  convert_r_to_xy(r,random_numbers[3],x,y,wgt);
  sample_energy(random_numbers[4],random_numbers[5],
                rprof_cdf_data[bin_sampled].ion_temp,energy);

  return 1;
}

/*
 * sample the cdf src_profile, to generate the sampled minor radius
 */
int RadialProfileSource::sample_source_rprofile(double rn_store1,
                                                double rn_store2,
                                                double &sampled_radius,
                                                int &sampled_bin)
{
  std::cout.precision(10);
  // loop over the bins
  for ( int i = 0 ; i < num_bins-1 ; i++ ) {
    if ( rn_store1 <= rprof_cdf_data[i].probability ) {
      sampled_radius = rprof_cdf_data[i].minor + ((rn_store2-0.5)*rprof_cdf_data[i].bin_width);

      if( i == num_bins-2 ) // temporary hack
	sampled_radius = rprof_cdf_data[i].minor;
      sampled_bin = i;
      return 1;
    }
  }
  std::cerr << "Failed to sample in source_profile!!" << std::endl;
  return 0;
}

//  convert the sampled radius to an rz coordinate by using plasma parameters
void RadialProfileSource::convert_rad_to_rz(double normalised_sample_r,
    double rn_store, double &radius, double &height) {

  double alpha = TWOPI*rn_store;
  //  std::cout << alpha << std::endl;
  double shift = shafranov_shift*(1.0-std::pow(normalised_sample_r,2));
  //  std::cerr << std::pow(*minor_sampled/(*minor_radius),2) << std::endl;
  radius = major_radius + shift + (minor_radius*normalised_sample_r*cos(alpha+(triangularity*sin(alpha))));
  height = elongation*normalised_sample_r*minor_radius*sin(alpha);

  //  std::cout << normalised_sample_r << " " << minor_radius << " " << radius << " " << height << std::endl;
};

// alternative
int RadialProfileSource::read_file()
{
  // open up the file
  std::ifstream rprof_input (source_file);
  std::string line;

  std::vector<rprofile_source_line> r_profile; // data read in
  std::string buf;

  // if the file has been opened
  if (rprof_input.is_open()) {
    // while there is data to read
    while (getline(rprof_input, line)) {
      std::vector<double> tokens;
      std::stringstream ss(line);
      
      while ( ss >> buf)
	tokens.push_back(atof(buf.c_str()));
      
      rprofile_source_line data;
      data.minor = tokens[0];
      data.ion_temp = tokens[1]/1000.; // convert to MeV
      data.d_dens = tokens[2];
      data.t_dens = tokens[3];
      
      // push the data back
      r_source_profile.push_back(data);
    }
    
    // now data is read in
    num_bins = r_source_profile.size();
    rprof_input.close();
    
   } else {
    std::cout << "Failed to open '"<< source_file <<"' is it in the right place?" << std::endl;
    return 0;
  }
  return 1;
}

// setup the rprofile source
int RadialProfileSource::setup_rprofile_source()
{
  // iterator for the radial source profile
  std::vector<rprofile_source_line>::iterator it;

  double total = 0.0;
  std::vector<double> total_rr; // running total of the reaction rate

  // step through the data to figure out the reaction rates
  for ( int i = 0 ; i < num_bins-1 ; i++ ) {
	  //  for ( it = r_source_profile.begin() ; it != r_source_profile.end() ; ++it ) {
    // take the ion temp and calc the cross section
    //    rprofile_source_line source_data = *it;
    rprofile_source_line source_data = r_source_profile[i];
    double nd_d = source_data.d_dens; // deuterium number density
    double nd_t = source_data.t_dens; // tritium number density
    double rr = nd_d*nd_t*dt_xs(source_data.ion_temp); // number of dt reactions
    total += rr; // take sum of total reaction rate
    total_rr.push_back(rr);
  }
  
  // step through the data once more to figre out the cdf
  for ( int i = 0 ; i < num_bins-1 ; i++ ) {
    double cdf = 0.0;
    rprofile_source_data data;
    for ( int j = 0 ; j <= i ; j++ )
      cdf += total_rr[j]/total;

    // tmp var
    int k = i+1;

    // store the cdf
    data.probability = cdf;
    data.ion_temp = r_source_profile[k].ion_temp;
    data.minor = r_source_profile[k].minor;
    data.bin_width = r_source_profile[k].minor - r_source_profile[k-1].minor;
    rprof_cdf_data.push_back(data);
    //    std::cout << r_source_profile[k].minor << " " << data.bin_width << " " <<  cdf << std::endl;

  }

  return 1;
}
