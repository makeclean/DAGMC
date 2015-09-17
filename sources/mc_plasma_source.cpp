#include "plasma_source.hpp"
#include "mc_plasma_source.hpp"
#include <iostream>

RadialProfileSource *rd_src;
ParametricSource *p_src;
RZProfileSource *rz_src;

int setting = -1;

// setup the appropriate source
void setup_(int *set, double rdums[50], int *erc) {
  int ec = 0;
  if (*set < 0 || *set > 3 ) {
    *erc = 0;
  }
  setting = *set;
  if(setting == 1) {
    p_src = new ParametricSource(1.0e20, //ion_density_ped
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
    ec = p_src->setup();
  }
  if(setting == 2) {
    rd_src = new RadialProfileSource(2.0, // minor radius
				     5.2, // major radius
				     1.0, // elongation
				     0.0, // triangularity
				     0.0, // shafranov shift
				     0.0, // start angle
				     360.0, // end angle
				     std::string("r_prof.txt"));
    ec = rd_src->setup();
  }
  if(setting == 3) {
    rz_src = new RZProfileSource(0,360.0,std::string("rzsource.txt"));
    ec = rz_src->setup();
  }
  *erc = ec;
  return;
}

// sample a valid position, wgt and direction
void sample_(double randoms[6], double *x, double *y, double *z, 
	     double *wgt, double *erg, int *ec) {
  int erc = 0;
  double x_c,y_c,z_c,w,e;
  if(setting == 1)
    erc = p_src->sample(randoms,x_c,y_c,z_c,w,e);
  if(setting == 2)
    erc = rd_src->sample(randoms,x_c,y_c,z_c,w,e);
  if(setting == 3)
    erc = rz_src->sample(randoms,x_c,y_c,z_c,w,e);

  *x = x_c;
  *y = y_c;
  *z = z_c;
  *wgt = w;
  *erg = e;
  *ec = erc;
}
