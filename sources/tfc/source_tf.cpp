#include <iostream>

#include "source_tf.h"

double ebins[47] = {1.0E-7,4.14e-07,8.76e-07,1.86e-06,3.93e-06,8.32e-06,1.76e-05,3.73e-05,
		    7.89e-05,0.000167,0.000354,0.000749,0.00158,0.00335,0.0071,0.015,
		    0.0318,0.0674,0.123,0.166,0.224,0.302,0.408,0.55,0.743,1,1.35,1.83,
		    2.47, 2.73,3.01,3.33,3.68,4.07,4.49,4.97,5.49,6.07,6.7,7.41,8.19,9.05,
		    10,11.1,12.2,13.8,14.2};

double fluxes[46] = {1.11182e+07,7.2947e+06,1.15913e+07,1.69177e+07,2.31575e+07,2.91347e+07,
		     3.41075e+07,4.68462e+07,5.93648e+07,5.75665e+07,7.3275e+07,9.36479e+07,
		     8.76419e+07,8.02542e+07,9.10929e+07,1.51276e+08,1.85725e+08,2.06504e+08,
		     1.35166e+08,1.29167e+08,1.63143e+08,1.76031e+08,1.40119e+08,1.7264e+08,
		     9.82925e+07,6.76057e+07,3.99465e+07,2.45524e+07,5.70741e+06,4.94299e+06,
		     4.32173e+06,3.59835e+06,2.92893e+06,2.52765e+06,2.34186e+06,2.024e+06,
		     1.92309e+06,1.91356e+06,1.98723e+06,1.98182e+06,2.10613e+06,2.24288e+06,
		     2.43055e+06,2.80265e+06,1.2578e+07,6.25074e+06};

double weights[46];
double pdf[46];

int bins = 46;
double total;
double tot_weight = 0.0;

int i;

void setup_()
{
  // calculate normalisation
  for (int i = 0 ; i < bins ; i++ )
    {
      total += fluxes[i];
      weights[i] = 0.0;
    }

  // calculate weighting factor
  for (int i = 0 ; i < bins ; i++ )
    {
      weights[i] = fluxes[i]/total;
      tot_weight += weights[i];
      pdf[i] = tot_weight;
      //      std::cout << i << " " << weights[i] << " " << pdf[i] << std::endl;
    }

  return;
}

/* sample a distribution linearly*/
void sample_linear_(double max, double min, double random, double &sample)
{
  sample = ((max-min)*random) + min;
  return;
}

/* sample the particle energy given the spectra */
void sample_(double rand1, double rand2, double &energy, double &weight)
{
  for ( i = 1 ; i < bins+1 ; i++ )
    {
      if ( i == 1 )
	{
	  if(rand1 <= pdf[i-1])
	    {
	      sample_linear_(ebins[i], ebins[i-1], rand2, energy);
	      weight = weights[i];
	      return;
	    }
	}
      else
	{
	  if(rand1 > pdf[i-1] && rand1 <= pdf[i] )
	    {
	      sample_linear_(ebins[i], ebins[i-1], rand2, energy);
	      weight = weights[i];
	      return;
	    }
	}
    }
}
