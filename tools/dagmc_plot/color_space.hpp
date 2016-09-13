#include <cmath>
#include <vector>

struct rgb_f {
    float r;
    float g;
    float b; 
  };

struct rgb_i {
    int r;
    int g;
    int b; 
  };

class colorspace {

  public:
  colorspace();
  ~colorspace();
  
  std::vector<rgb_i> getRGBtable(unsigned int num_divisions);

  private:
  double c_value(double v, double s);
  double x_value(double c, double h);
  double m_value(double v, double c);

  void get_hsv();

  rgb_f getRGBfromHSV(double h, double s, double v);
  rgb_f getRGBprime(double c, double x, double h);
  rgb_f getRGBvalue(rgb_f rgbprimed, double m);
  
  rgb_i rgbFloatToInt(rgb_f float_rgb);  
};

  
