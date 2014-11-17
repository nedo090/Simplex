#include<stdlib.h>


double **inversion(double **base)
{
	double**inverse = (double**) malloc(2*sizeof(double *));
	double det=base[0][0]*base[1][1]-base[1][0]*base[0][1];
	inverse[0]=(double *)malloc(2*sizeof(double));
	inverse[1]=(double *)malloc(2*sizeof(double));
	inverse[0][0]= base[1][1]/det;
	inverse[0][1]= -1*base[1][0]/det;
	inverse[1][0]= -1*base[0][1]/det;
	inverse[1][1]= base[0][0]/det;
	return inverse;
}


