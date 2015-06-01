
#ifdef __cplusplus
extern "C" {
#endif

/**
 *  setup the problem arrays
 */
void setup_();

/**
 * sample between max and min linearly 
 * double max - maximum bin bound
 * double min - minimum bin bound
 * double random - random number between 0 and 1
 * double sample - the returned value
 */
void linear_sample_(double &max, double &min, 
		    double &random, double &sample); 

/**
 * samples the build in source distributions
 * double rand1 - a random number between 0 and 1
 * double rand2 - a random number between 0 and 1
 * double energy - the returned energy 
 * double weight - the returned statistical weight
 */
void sample_(double &rand1, double &rand2, double &energy,
	     double &weight);

#ifdef __cplusplus
}
#endif
