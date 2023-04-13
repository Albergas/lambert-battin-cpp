
#include <cmath>
#include "vector.h"

// Returns the norm of a given vector v of dimension n
double norm(double v[], int n) {
	
	double suma = 0;
	int i;
	
	if(n <= 0)
		throw "Empty vector";

	for(i = 0; i < n; i++)
		suma += v[i]*v[i];

	return(sqrt(suma));
}

// Returns the dot product of given vectors v, w; both of dimension n
double dot(double v[], double w[], int n) {

	double suma = 0;
	int i;

	if(n <= 0)
		throw "Empty vector";

	for(i = 0; i < n; i++)
		suma += v[i]*w[i];

	return suma;
}

// Computes the cross product of given vectors v, w (both of dimension 3) and stores it in vector u
void cross(double v[], double w[], double u[]) {
	u[0] = v[1] * w[2] - v[2] * w[1];
	u[1] = v[2] * w[0] - v[0] * w[2];
	u[2] = v[0] * w[1] - v[1] * w[0];
}