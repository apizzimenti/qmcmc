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

// Toggle outputting debug statements to a 'log.txt' file. 1 is assert, 0 is deassert.
#define DEBUG 1

#if DEBUG
	#define DEBUG_CMD(...) do {__VA_ARGS__;} while(0)
#else
	#define DEBUG_CMD(...)
#endif

#endif // PARALLEL_TEMPERING_H
