#include <iostream>
#include <cmath>
#include "lambert.h"

using namespace std;

int main() {

	double ro[3] = {20.0e6, 20.0e6, 0};
	double r[3] = {-20.0e6, 10.0e6, 0};
	double tof = 1.0 * 86400;
	char dm[] = "retro";

	double v1[3] = {0, 0, 0};
	double v2[3] = {0, 0, 0};

	lambert(ro, r, dm, tof, v1, v2);

	cout.precision(15);

	cout << v1[0] << ", " << v1[1] << ", " << v1[2] << endl;
	cout << v2[0] << ", " << v2[1] << ", " << v2[2] << endl;

	double v1_test[3] = {4144.30717367665, -1571.15318557575, 0};
	double v2_test[3] = {3223.39508300486, 4103.76281774997, 0};

	double max_err = pow(10, -10);

	double e[6] = { fabs(v1[0] - v1_test[0]), fabs(v1[1] - v1_test[1]), fabs(v1[2] - v1_test[2]), fabs(v2[0] - v2_test[0]), fabs(v2[1] - v2_test[1]), fabs(v2[2] - v2_test[2]) };

	cout << e[0] << ", " << e[1] << ", " << e[2] << ", " << e[3] << ", " << e[4] << ", " << e[5] << endl;

	if( e[0] < max_err && e[1] < max_err && e[2] < max_err && e[3] < max_err && e[4] < max_err && e[5] < max_err )
		cout << "Test lambert(): passed" << endl;
	else
		cout << "Test lambert(): failed" << endl;

	return 0;
}
