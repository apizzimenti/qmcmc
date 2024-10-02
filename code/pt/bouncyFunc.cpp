#include <random>
#include <math.h>

double bouncyFunctionCost(double x, double imstupidanddontdoanything) {
	if (x < 0 || x > 3.14) {
		return 1000;
	}
	double y = 0.0;
	for (int n = 0; n < 5; n++) {
		y += (10.0-n)*sin(20.0*n*(1.0/5.0));
		y += (10.0-n)*cos(20.0*n*(1.0/5.0));
	}
	return y;
}

double calculate_expression(double x, double dummy) {
    double result;
    double term1 = 0.5 * pow((x - 2), 4);
    double term2 = 2 * pow((x - 2), 3);
    double term3 = 3 * (x - 2);
    result = term1 - term2 + term3;
    return result;
}


double update(double x) {
    static std::mt19937 rng(std::random_device{}()); // Random number generator
    std::uniform_real_distribution<double> uni(0.0, 1.0);
	if (uni(rng) < .5) {
		return x += .1;
	}
	else {
		return x -= .1;
	}
}
