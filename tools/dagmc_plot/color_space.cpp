#include "color_space.hpp"
#include <iostream> 

// constructor
colorspace::colorspace() {
}

// destructor
colorspace::~colorspace() {
}

// Get the vector of RGB colors
std::vector<rgb_i> colorspace::getRGBtable(unsigned int num_divisions) {
   std::vector<rgb_f> rgb_table_f;
   std::vector<rgb_i> rgb_table;
  // loop over h and s, keep v constant
  double bin_width = 360.0/float(num_divisions);
  // build floating point list of rgb colors
  for ( unsigned int i = 0 ; i < num_divisions ; i++ ) {
      double h = bin_width*float(i);
      rgb_f value = getRGBfromHSV(h,1.0,1.0);
      rgb_table.push_back(rgbFloatToInt(value));
      std::cout << rgb_table[i].r << " " << rgb_table[i].g << " " << rgb_table[i].b << std::endl;
  } 
  return rgb_table;
}
/// End of public functions ///
/// Private ///

// given value and saturation get c
double colorspace::c_value(double v, double s) {
    double c = v*s;
    return c;
}

// given c and hue get x value
double colorspace::x_value(double c, double h) {
    double value = std::abs(fmod(h/60.,2.)-1.0);
    double x = c*(1-value);
    return x;
}

// given value and c determine m
double colorspace::m_value(double v, double c) {
    double m = v-c;
    return m;
}

// for a given hsv value get 
rgb_f colorspace::getRGBfromHSV(double h, double s, double v) {
  double c = c_value(v,s);
  double x = x_value(c,h);
  double m = m_value(v,c);

  rgb_f rgb_primed = getRGBprime(c,x,h);
  rgb_f rgb_value = getRGBvalue(rgb_primed,m);
  return rgb_value; 
}

// given c,x,h get r' g' b' values
rgb_f colorspace::getRGBprime(double c, double x, double h) {
    rgb_f rgb_primed;
    if ( h < 60.) {
        rgb_primed.r = c;
        rgb_primed.g = x;
        rgb_primed.b = 0.0;
    } else if ( h >= 60. && h < 120. ) {
        rgb_primed.r = x;
        rgb_primed.g = c;
        rgb_primed.b = 0.0;
    } else if ( h >= 120. && h < 180. ) {
        rgb_primed.r = 0.0;
        rgb_primed.g = c;
        rgb_primed.b = x;
    } else if ( h >= 180. && h < 240. ) {
        rgb_primed.r = 0.0;
        rgb_primed.g = x;
        rgb_primed.b = c;
    }  else if ( h >= 240. && h < 300. ) {
        rgb_primed.r = x;
        rgb_primed.g = 0.0;
        rgb_primed.b = c;
    }  else if ( h >= 300. && h < 360. ) {
        rgb_primed.r = c;
        rgb_primed.g = 0.0;
        rgb_primed.b = x;
    }
    return rgb_primed;
}


// given r'g'b' and m get the true rgb
rgb_f colorspace::getRGBvalue(rgb_f rgb_primed, double m) {
    rgb_f rgb_data;
    rgb_data.r = (rgb_primed.r + m)*255.0;
    rgb_data.g = (rgb_primed.g + m)*255.0;
    rgb_data.b = (rgb_primed.b + m)*255.0;
    return rgb_data;
}

// convert a floating point rgb to int rgb
rgb_i colorspace::rgbFloatToInt(rgb_f float_rgb) {
  rgb_i int_rgb;
  int_rgb.r = int(float_rgb.r);
  int_rgb.g = int(float_rgb.g);
  int_rgb.b = int(float_rgb.b);
  return int_rgb;
}