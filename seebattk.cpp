
#include <cmath>
#include "seebattk.h"

double seebattk(double v){

	double d[22];
	d[1] =     1.0e0 /    3.0e0;
	d[2] =     4.0e0 /   27.0e0;
	d[3] =     8.0e0 /   27.0e0;
	d[4] =     2.0e0 /    9.0e0;
	d[5] =    22.0e0 /   81.0e0;
	d[6] =   208.0e0 /  891.0e0;
	d[7] =   340.0e0 / 1287.0e0;
	d[8] =   418.0e0 / 1755.0e0;
	d[9] =   598.0e0 / 2295.0e0;
	d[10] =  700.0e0 / 2907.0e0;
	d[11]=   928.0e0 / 3591.0e0;
	d[12]=  1054.0e0 / 4347.0e0;
	d[13]=  1330.0e0 / 5175.0e0;
	d[14]=  1480.0e0 / 6075.0e0;
	d[15]=  1804.0e0 / 7047.0e0;
	d[16]=  1978.0e0 / 8091.0e0;
	d[17]=  2350.0e0 / 9207.0e0;
	d[18]=  2548.0e0 /10395.0e0;
	d[19]=  2968.0e0 /11655.0e0;
	d[20]=  3190.0e0 /12987.0e0;
	d[21]=  3658.0e0 /14391.0e0;

	double sum1 = d[1];
	double delold = 1.0e0;
	double termold = d[1];

	int i = 1;
	double del, term;

	while(true){

		del = 1.0e0 / (1.0e0 + d[i+1]*v*delold);
		term = termold * (del - 1.0e0);
		sum1 += term;
		i = i + 1;
		delold = del;
		termold = term;
		if ( (i > 20) || (fabs(termold) > 0.000001e0))
			break;
	}

	return sum1;
}
