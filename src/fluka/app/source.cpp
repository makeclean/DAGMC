#include "plasma_source.hpp"

plasma_source::PlasmaSource* PlasmaSource;

void setup_source_(double &major_radius,
		   double &minor_radius,
		   double &elongation,
		   double &shafranov_shift,
		   double &triangularity,
		   double &ion_density_core,
		   double &ion_density_ped,
		   double &ion_density_sep,
		   double &ion_temp_core,
		   double &ion_temp_ped,
		   double &ion_temp_sep,
		   double &ion_density_peak,
		   double &ion_temp_peak,
		   double &ion_temp_beta) {

  if ( PlasmaSource == NULL ) {
    PlasmaSource = new plasma_source::PlasmaSource(ion_density_ped,
				    ion_density_sep,
				    ion_density_core,
				    ion_temp_ped,
				    ion_temp_sep,
				    ion_temp_core,
				    ion_density_peak,
				    ion_temp_peak,
				    ion_temp_beta,
				    minor_radius,
				    major_radius,
				    elongation,
				    triangularity,
				    shafranov_shift);
       	    			  
  }
}
