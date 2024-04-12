#ifndef PARALLEL_TEMPERING_H
#define PARALLEL_TEMPERING_H

#include <vector>
#include <cmath>
#include <random>

// Assuming delta's type is std::vector<double>. Adjust if necessary.
void ParallelTempering(std::vector<double>& delta, std::vector<double>& T, int L, std::vector<int>& R, FILE* log);

// Assuming the return type is void for generateNeighbor and the arguments are as needed. Adjust if necessary.
// Also, assuming generateNeighbor modifies its argument by reference.
void generateNeighbor(double* value);

// Assuming temperedCost takes a double and returns a double. Adjust as necessary.
double temperedCost(double value);

#endif // PARALLEL_TEMPERING_H
