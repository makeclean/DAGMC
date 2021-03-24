#include <iostream>
#include <array>
#include <vector>

#include "plasma_data.hpp"

#ifndef PLASMA_SOURCE_HPP
#define PLASMA_SOURCE_HPP 1

#ifdef __cplusplus
extern "C" {
#endif

void setup_(const double &ionDensityOrigin,
	    const double &ionTempOrigin,
	    const double &ionTempPeak,
	    const double &majorRadius,
	    const double &minorRadius,
	    const double &elongation,
	    const double &triangularity,
	    const double &shafranov);

void sample_(double *randoms,
	     double &x,
	     double &y,
	     double &z,
	     double &u,
	     double &v,
	     double &w,
	     double &e);

#ifdef __cplusplus
}
#endif

namespace plasma_source {
  
struct xs_params {
  double c[7];
};

enum reaction_type {DD,DT};

  /*
// Derived Classes should define their 
class PlasmaSourceDTHMode : public PlasmaSource {
 
};
  */
  
class PlasmaSource {
  public:
    // constructor
    PlasmaSource();

    // destructor
    ~PlasmaSource();

    // large constructor
    PlasmaSource(const double ion_density_ped,
                 const double ion_density_sep,
		 const double ion_density_origin,
                 const double ion_temp_ped,
		 const double ion_temp_sep,
                 const double ion_temp_origin, 
		 const double pedestal_rad,
                 const double ion_density_peak,
		 const double ion_temp_peak,
                 const double ion_temp_beta,
                 const double minor_radius,
                 const double major_radius,
                 const double elongation,
                 const double triangularity,
                 const double shafranov,
                 const std::string plasma_type = "",
                 const int plasma_id = 0,
                 const int number_of_bins = 100,
		 const double min_toroidal_angle = 0.0,
		 const double max_toridal_angle = 360.);

    // main sample fucnction
    void sample(std::array<double,8> randoms,
                double &x,
                double &y,
                double &z,
                double &u,
                double &v,
                double &w,
                double &E);

    /*
     * Function to setup the plasma source in the first case.
     */
    void setup_plasma_source();

    /*
     * function to calculate the ion density at a specific minor 
     * radius
     */
    double ion_density(const double sample_radius);

    /*
     * function to calculate the ion temperature at a specific minor 
     * radius
     */
    double ion_temperature(const double sample_radius);

    /*
     * function to determine the value of the dt xs cross sections at 
     * a specific ion temp
     */
    double dt_xs(double ion_temp);

    /*
     * function to determine the value of the dt xs cross sections at 
     * a specific ion temp
     */
    double dd_xs(double ion_temp);

    protected:
    /*
     * sample the source, returns the minor radius sampled
     * expects new rn_store every call
     */
    void sample_radial(double rn_store1,
                       double rn_store2,
                       double &sampled_radius,
                       int &sampled_bin);

    /*
     * sample the neutron energy  in MeV
     */
    void sample_energy(const int bin_number, double random_number1, double random_number2,
		                   double &energy_neutron);

    /*
     * take the sampled minor radius and convert to cylindrical coordinates
     */
    void convert_rad_to_rz(const double minor_sampled,
                           const double rn_store, 
                           double &radius,
                           double &height);

    /*
     * convert partial cylindrical coords to xyz
     */
    void convert_r_to_xy(const double r, const double rn_store, double &x, double &y);

    /*
     * get an isotropically direction vector
     */
    void isotropic_direction(const double random1,
                             const double random2,
		                         double &u,
                             double &v,
                             double &w);

    public:
    /*
     * get a key-value pair string representation of this instance of the source
     */
    std::string to_string();

    /*
     * create a new source from the provided key-value pair string representation
     * of the source
     */
    static PlasmaSource from_string(std::string parameters);


    void setIonDensityPedestal(const double value) {
      ionDensityPedestal = value;
    }

    /*
     * Getter for pedestal ion density
     */
    const double &getIonDensityPedestal() const { return ionDensityPedestal; }

    void setIonDensitySeparatrix(const double value) {
      ionDensitySeparatrix = value;
    }
  
    /*
     * Getter for separatrix ion density
     */
    const double &getIonDensitySeparatrix() const { return ionDensitySeparatrix; }

  
    void setIonDensityOrigin(const double value) {
      ionDensityOrigin = value;
    }
  
    /*
     * Getter for origin ion density
     */
    const double &getIonDensityOrigin() const { return ionDensityOrigin; }

    void setIonTemperaturePedestal(const double value) {
      ionTemperaturePedestal = value;
    }
  
    /*
     * Getter for pedestal ion temperature
     */
    const double &getIonTemperaturePedestal() const { return ionTemperaturePedestal; }

    void setIonTemperatureSeparatrix(const double value) {
      ionTemperaturePedestal = value;
    }
  
    /*
     * Getter for separatrix ion temperature
     */
    const double &getIonTemperatureSeparatrix() const { return ionTemperatureSeparatrix; }

    void setIonTemperatureOrigin(const double value) {
      ionTemperatureOrigin = value;
    }

    /*
     * Getter for origin ion temperature
     */
    const double &getIonTemperatureOrigin() const { return ionTemperatureOrigin; }

    void setSetPedestalRadius(const double value) {
      pedestalRadius = value;
    }
  
    /*
     * Getter for pedestal radius
     */
    const double &getPedestalRadius() const { return pedestalRadius; }

    void setIonDensityPeaking(const double value) {
      ionDensityPeaking = value;
    }
  
    /*
     * Getter for ion density peaking factor
     */
    const double &getIonDensityPeaking() const { return ionDensityPeaking; }

    void setIonTemperaturePeaking(const double value) {
      ionTemperaturePeaking = value;
    }
  
    /*
     * Getter for ion temperature peaking factor
     */
    const double &getIonTemperaturePeaking() const { return ionTemperaturePeaking; }

    void setIonTemperatureBeta(const double value) {
      ionTemperatureBeta = value;
    }
 
    /*
     * Getter for ion temperature beta factor
     */
    const double &getIonTemperatureBeta() const { return ionTemperatureBeta; }

    void setMinorRadius(const double value) {
      minorRadius = value;
    }
 
    /*
     * Getter for minor radius
     */
    const double &getMinorRadius() const { return minorRadius; }

    void setMajorRadius(const double value) {
      majorRadius = value;
    }
  
    /*
     * Getter for major radius
     */
    const double &getMajorRadius() const { return majorRadius; }

    void setElongation(const double value) {
      elongation = value;
    }
  
    /*
     * Getter for elongation
     */
    const double &getElongation() const { return elongation; }

    void setTriangularity(const double value) {
      triangularity = value;
    }
  
    /*
     * Getter for triangularity
     */
    const double &getTriangularity() const { return triangularity; }

    void setShafranovShift(const double value) {
      shafranov = value;
    }
   
    /*
     * Getter for shafranov shift
     */
    const double &getShafranov() const { return shafranov; }

    void setMinToroidalAngle(const double value) {
      minToroidalAngle = value;
    }
  
    /*
     * Getter for minimum toroidal angle
     */
    const double &getMinToroidalAngle() const { return minToroidalAngle; }

    void setMaxToroidalAngle(const double value) {
      maxToroidalAngle = value;
    }
  
    /*
     * Getter for maximum toroidal angle
     */
    const double &getMaxToroidalAngle() const { return maxToroidalAngle; }

    void setPlasmaType(const double value) {
      plasmaType = value;
    }
  
    /*
     * Getter for plasma type
     */
    const std::string &getPlasmaType() const { return plasmaType; }

    /*
     * set the fusion reaction type (dd or dt)
     */
    void setReactionType(const int value) {
      if (value == 0) reactionType = DD;
      else if (value == 1) reactionType = DT;
      else reactionType = DT;
    }
  
    /*
     * Getter for reaction type (dd or dt)
     */
    const int &getReactionType() const { return reactionType; }

    /*
     * Setter for plasma identifier
     */
    void setPlasmaID(const int value ) { plasmaId = value; }
  
    /*
     * Getter for plasma identifier
     */
    const int &getPlasmaId() const { return plasmaId; }

    void setNumberOfBins(int value) { numberOfBins = value; }
  
    /*
     * Getter for number of bins
     */
    const int &getNumberOfBins() const { return numberOfBins; }

  protected:
    std::vector<double> source_profile;
    std::vector<double> ion_kt;

    double ionDensityPedestal;
    double ionDensitySeparatrix;
    double ionDensityOrigin;
    double ionTemperaturePedestal;
    double ionTemperatureSeparatrix;
    double ionTemperatureOrigin;
    double pedestalRadius;
    double ionDensityPeaking;
    double ionTemperaturePeaking;
    double ionTemperatureBeta;
    double minorRadius;
    double majorRadius;
    double elongation;
    double triangularity;
    double shafranov;
    double minToroidalAngle;
    double maxToroidalAngle;

    std::string plasmaType;
    reaction_type reactionType;
  
    int plasmaId;
    double binWidth;
    int numberOfBins;
};

}// end of namespace

#endif
