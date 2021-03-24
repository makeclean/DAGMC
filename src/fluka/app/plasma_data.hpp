#include <cmath>

#define U 1.66053906660e-27
#define NEUTRON_MASS 1.00866491588*U
#define DEUTERON_MASS 2.01410177811*U
#define TRITON_MASS 3.0160492*U
#define DT_MASS DEUTERON_MASS+TRITON_MASS
#define DD_MASS 2*DEUTERON_MASS
#define DT_ENERGY 14.08
#define DD_ENERGY 2.54
#define DD_FACTOR 4*NEUTRON_MASS*DD_ENERGY/DD_MASS
#define DT_FACTOR 4*NEUTRON_MASS*DT_ENERGY/DT_MASS

#ifndef PLASMA_DATA_HPP
#define PLASMA_DATA_HPP 1

struct nrl_formulary_xs {
  double a1;
  double a2;
  double a3;
  double a4;
  double a5;
};

// nrl formulary data from https://en.wikipedia.org/wiki/Nuclear_fusion#Formulas_of_fusion_cross_sections
  
nrl_formulary_xs dtxs = {45.95,50200,1.368e-2,1.076,409};
nrl_formulary_xs dd1xs = {46.097,372,4.36e-4,1.22,0};
nrl_formulary_xs dd2xs = {47.88,482,3.08e-4,1.177,0};
nrl_formulary_xs ttxs = {38.39,448,1.02e-3,2.09,0};

// plasma type 
struct plasma_type {
  nrl_formulary_xs xs;
  double factor;
};

plasma_type dt_plasma = {dtxs,std::sqrt(DT_FACTOR)};
plasma_type dd_plasma = {dd2xs,std::sqrt(DD_FACTOR)};

  
// given an energy in kev sample the xs
double sample_nrl_formulary(nrl_formulary_xs xs, double e) {
  double numerator = xs.a5 + xs.a2*std::pow(std::pow(xs.a4-xs.a3*e,2)+1,-1);
  double denominator = e*(std::exp(xs.a1/std::sqrt(e))-1);
  return numerator/denominator;
}

#endif
