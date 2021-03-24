#include <vector>
#include <array>

enum Sampling { NEAREST, BILINEAR_INTERPOLATION };

#ifdef __cplusplus
extern "C" {
#endif

void setup_();

void get_field_(double &x, double &y, double &z,
		double &bx, double &by, double &bz, double &b);
  
#ifdef __cplusplus
} // extern "C"                                                                            
#endif


  
class MagneticField {
  public:
  MagneticField();
  ~MagneticField();

  public:
  void load_file(std::string filename = "../mast-u-fields.txt");
  
  void process_field();

  double bilinear_interpolation(const std::array<double,2> x_stencil,
			       const std::array<double,2> y_stencil,
			       const std::array<double,4> values,
			       const double x,
			       const double y);
  
  void interpolate_field(const double x,
			 const double y,
			 const double z,
			 double &field_x,
			 double &field_y,
			 double &field_z);
  
  private:
  void find_indices(const double r,
		    const double z,
		    int &r_ind,
		    int &z_ind);
 
  void field_by_nearest(const int index,
			double &field_x,
			double &field_y,
			double &field_z);

  void field_by_interpolation(const double r,
			      const double z,
			      const int r_index,
			      const int z_index,
			      double &field_x,
			      double &field_y,
			      double &field_z);

  int index_by_rz(const int r_index,
		  const int z_index);
  
  private:
  Sampling mode;
  int n_r; // number of r points
  int n_z; // number of z points
  std::vector<double> r_point; // the unique grid points in r
  std::vector<double> z_point; // the unique grid points in z
  
  public:
  std::vector<double> radial_pos; /// temp
  std::vector<double> vertical_pos; /// temp
  std::vector<std::array<double,3> > b_field;
  
};
