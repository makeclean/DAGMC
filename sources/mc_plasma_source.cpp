#include "plasma_source.hpp"
#include "mc_plasma_source.hpp"

RadialProfileSource *rd_src;
ParametricSource *p_src;
RZProfileSource *rz_src;

int setting = -1;

// setup the appropriate source
void setup_(int idums[50], double rdums[50], int *erc) {
  int ec = 0;
  if (idums[0] < 0 || idums[0] > 3 ) {
    *erc = 0;
  }
  setting = idums[0];
  if(setting == 1 || setting == 0 ) {
    p_src = new ParametricSource(rdums[0], // ion density ped
				 rdums[1], // ion density sep
				 rdums[2], // ion density org
				 rdums[3], // ion temp ped
				 rdums[4], // ion temp sep
				 rdums[5], // ion temp org
				 rdums[6], // pedistal radius
				 rdums[7], // ion density peak
				 rdums[8], // ion temp peaking
				 rdums[9], // minor radius
				 rdums[10], // major radius
				 rdums[11], // elongation
				 rdums[12], // triangularity
				 rdums[13], // shafranov shift
				 "", 
				 idums[0], // plasma selector
				 idums[1], // number of bins
				 rdums[14], // start angle
				 rdums[15]); //end angle
    ec = p_src->setup();
  }
  if(setting == 3) {
    rd_src = new RadialProfileSource(rdums[0], // minor radius
				     rdums[1], // major radius
				     rdums[2], // elongation
				     rdums[3], // triangularity
				     rdums[4], // shafranov shif
				     rdums[5], // start angle
				     rdums[6], // end angle
				     std::string("r_prof.txt"));
    ec = rd_src->setup();
  }
  if(setting == 2) {
    rz_src = new RZProfileSource(rdums[0],
				 rdums[1],
				 std::string("rzsource.txt"));
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
  if(setting == 1 || setting == 0)
    erc = p_src->sample(randoms,x_c,y_c,z_c,w,e);
  if(setting == 3)
    erc = rd_src->sample(randoms,x_c,y_c,z_c,w,e);
  if(setting == 2)
    erc = rz_src->sample(randoms,x_c,y_c,z_c,w,e);

  *x = x_c;
  *y = y_c;
  *z = z_c;
  *wgt = w;
  *erg = e;
  *ec = erc;
}
