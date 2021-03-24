#include <iostream>
#include <fstream>
#include <set>
#include <algorithm>
#include "bfield.hpp"
#include <cmath>

MagneticField field = MagneticField();
/*
int main() {
  setup_();
  double bxn,byn,bzn,b;
  get_field_(0.113,0.,-2.08,
	     bxn,byn,bzn,b);
}
*/

void get_field_(double &x, double &y, double &z,
		double &bxn, double &byn, double &bzn, double &b){

  double bfield[3];
  field.interpolate_field(x,y,z,bfield[0],bfield[1],bfield[2]);
  double bnorm = std::sqrt(bfield[0]*bfield[0] +
			   bfield[1]*bfield[1] +
			   bfield[2]*bfield[2]);
  //  std::cout << x << " " << y << " " << z << std::endl;
  //std::cout << bfield[0] << " " << bfield[1] << " " << bfield[2] << std::endl;
  // fluka expects this
  if(bnorm == 0.0) {
    bxn = 0;
    byn = 0;
    bzn = 1.0;
    b = 0;
  } else {
    bxn = bfield[0]/bnorm;
    byn = bfield[1]/bnorm;
    bzn = bfield[2]/bnorm;
    b = bnorm;
  }
  return;
} 

// Constructor
MagneticField::MagneticField() {
  mode = NEAREST;
  load_file();
}

// Desctructor
MagneticField::~MagneticField() {
}

// load data from file
void MagneticField::load_file(std::string filename) {
  std::ifstream input;
  input.open(filename, std::ifstream::in);

  if(!input.good()) {
    std::cout << "Failed to open file" << std::endl;
    exit(10);
  }
  
  double r,z,br,bz,btheta;
  // while we can read
  while ( input.good() ) {
    input >> r >> z >> br >> bz >> btheta;
    // push_back to the storage array * cm conversion
    radial_pos.push_back(r*100.0);
    vertical_pos.push_back(z*100.0);
    // note the order br bz btheta
    b_field.push_back({br,bz,btheta});
  }
  
  return process_field();
}

// process the magnetic field
void MagneticField::process_field() {
  std::set<double> unique_r(radial_pos.begin(), radial_pos.end());
  std::set<double> unique_z(vertical_pos.begin(), vertical_pos.end());
 
  // set the size of the data
  n_r = unique_r.size();
  n_z = unique_z.size();

  //  std::cout << "bins: " << n_r << " " << n_z << std::endl;
  
  // copy in the set data
  std::copy(unique_r.begin(), unique_r.end(),std::back_inserter(r_point));
  std::copy(unique_z.begin(), unique_z.end(),std::back_inserter(z_point));

  // sort the vectors
  std::sort(r_point.begin(), r_point.end());
  std::sort(z_point.begin(), z_point.end());
  
  return;
}

// find out which radial and vertical bins are nearest
void MagneticField::find_indices(const double r, const double z, int &r_ind, int &z_ind) {
  std::vector<double>::iterator it_r = std::lower_bound(r_point.begin(),
							r_point.end(), r);
  r_ind = it_r - r_point.begin();
  std::vector<double>::iterator it_z = std::lower_bound(z_point.begin(),
							z_point.end(), z);
  z_ind = it_z - z_point.begin();
  return;
};

double MagneticField::bilinear_interpolation(const std::array<double,2> x_stencil,
					     const std::array<double,2> y_stencil,
					     const std::array<double,4> values,
					     const double x,
					     const double y) {
  // perform binlinear interpolation using f(x,y)  = a0 + a1*x + a2*y + a3*x*y

  // x1 = x_stencil[0]
  // x2 = x_stencil[1]
  // y1 = y_stencil[0]
  // y2 = y_stencil[1]
  // fQ11 = values[0]
  // fQ12 = values[1]
  // fQ21 = values[2]
  // fQ22 = values[3]
  
  double x1x2y1y2 = (x_stencil[0] - x_stencil[1])*(y_stencil[0] - y_stencil[1]);
  double x1x2y2y1 = (x_stencil[0] - x_stencil[1])*(y_stencil[1] - y_stencil[0]);

  double a0 = (values[0]*x_stencil[1]*y_stencil[1]/x1x2y1y2) +
              (values[1]*x_stencil[1]*y_stencil[0]/x1x2y2y1) +
              (values[2]*x_stencil[0]*y_stencil[1]/x1x2y2y1) +
              (values[3]*x_stencil[0]*y_stencil[0]/x1x2y1y2);

  double a1 = (values[0]*y_stencil[1]/x1x2y1y2) +
              (values[1]*y_stencil[0]/x1x2y2y1) +
              (values[2]*y_stencil[1]/x1x2y2y1) +
              (values[3]*y_stencil[0]/x1x2y1y2);

  double a2 = (values[0]*x_stencil[1]/x1x2y1y2) +
              (values[1]*x_stencil[1]/x1x2y2y1) +
              (values[2]*x_stencil[0]/x1x2y2y1) +
              (values[3]*x_stencil[0]/x1x2y1y2);

  double a3 = (values[0]/x1x2y1y2) +
              (values[1]/x1x2y2y1) +
              (values[2]/x1x2y2y1) +
              (values[3]/x1x2y1y2);

  double value = a0 + a1*x + a2*y + a3*x*y;
  return value;
}

// get the magnetic field by lookup to nearest point
void MagneticField::field_by_nearest(const int index,
				     double &field_r,
				     double &field_z,
				     double &field_t) {

    field_r = b_field[index][0];
    field_z = b_field[index][1];
    field_t = b_field[index][2];
    return;
}

void MagneticField::field_by_interpolation(const double r,
					   const double z,
					   const int r_index,
					   const int z_index,
					   double &field_x,
					   double &field_y,
					   double &field_z) {
    
    std::array<double,2> x_pts;
    std::array<double,2> y_pts;
    
    // set the r coordinate values
    // if radial value returned is the 1st bin
    // resort to 2d interpolation
    if ( r_index == 0 ) {
      x_pts[0] = r_point[r_index];
      x_pts[1] = r_point[r_index];
    } else if ( r_index > n_r ) {
      // if the r_index is greater than the
      // number of indices
      x_pts[0] = r_point[r_index];
      x_pts[1] = r_point[r_index];
    } else {
      // standard case
      x_pts[0] = r_point[r_index-1];
      x_pts[1] = r_point[r_index];
    }
    
    // set the z coordinate values
    // if radial value returned is the 1st bin
    // resort to 2d interpolation
    if ( z_index == 0 ) {
      y_pts[0] = z_point[z_index];
      y_pts[1] = z_point[z_index];
    } else if ( z_index > n_z ) {
      // if the r_index is greater than the
      // number of indices
      y_pts[0] = z_point[z_index];
      y_pts[1] = z_point[z_index];
    } else {
      // standard case
      y_pts[0] = z_point[z_index-1];
      y_pts[1] = z_point[z_index];
    }
    
    std::array<double,4> values;
    int index = (r_index-1*n_z) + z_index-1;
    values[0] = b_field[index][0];
    index = (r_index-1*n_z) + z_index;
    values[1] = b_field[index][0];
    index = (r_index*n_z) + z_index-1;
    values[2] = b_field[index][0];
    index = (r_index-1*n_z) + z_index;
    values[3] = b_field[index][0];
    
    // b field x
    double value = bilinear_interpolation(x_pts,y_pts,values,r,z);
    
    //    std::cout << values[0] << " " << values[1] << " " << values[2] << " " << values[3] << " " << value << std::endl;

}

// array index given radial and vertical index
int MagneticField::index_by_rz(const int r_index, const int z_index) {
  int index = (r_index)*n_z + (z_index) ;
  return index;
}

// given some position, iterpolate the field from the nearest points
void MagneticField::interpolate_field(const double x,
				      const double y,
				      const double z,
				      double &field_x,
				      double &field_y,
				      double &field_z) {
  // determine the radial value
  double r = std::sqrt(x*x + y*y);

  //  std::cout << "interpolat_field: r ";
  //std::cout << r << " " << r_point.front() << " " << r_point.back() << std::endl;
  
  // early return ensure r and z are within bounds else return 0
  if ( r < r_point.front() || r > r_point.back() ) {
    field_x = 0.;
    field_y = 0.;
    field_z = 0.;
    return;
  }
  //  std::cout << "interpolat_field: z ";
  // std::cout << z << " " << z_point.front() << " " << z_point.back() << std::endl;
  
  // early return ensure z is in bounds
  if ( z < z_point.front() || z > z_point.back() ) {
    field_x = 0.;
    field_y = 0.;
    field_z = 0.;
    return;
  }

  int r_index,z_index;
  // find the index
  find_indices(r,z,r_index,z_index);

  double field_r, field_t;
  // note passing field_r and field_t below so we can rotate the field
  if( mode == NEAREST ) {
    field_by_nearest(index_by_rz(r_index,z_index),field_r,field_z,field_t);
  }
  if( mode == BILINEAR_INTERPOLATION ) {
    field_by_interpolation(r, z, r_index, z_index, field_r, field_z, field_t);
  }

  // need to rotate the b-field appropraitely
  // determine theta angle
  double theta = std::atan2(y,x);
  //
  field_x = field_r*std::cos(theta) - field_t*std::sin(theta);
  field_y = field_r*std::sin(theta) + field_t*std::cos(theta);
  
  return;
}
