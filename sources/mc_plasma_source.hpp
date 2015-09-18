#ifdef __cplusplus
extern "C" {
#endif

// setup the appropriate source
void setup_(int idums[50], double rdums[50], int *ec);

// sample from the source
void sample_(double randoms[6], double *x, double *y, double *z, 
	     double *wgt, double *erg, int *ec);

#ifdef __cplusplus
} // extern C
#endif
