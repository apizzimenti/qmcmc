#ifndef PARALLEL_TEMPERING_H
#define PARALLEL_TEMPERING_H

#include <vector>
#include <cmath>
#include <random>

// Assuming delta's type is std::vector<double>. Adjust if necessary.
template<class V>
std::vector<V> parallelTempering(int n, int m, std::vector<double>& temperatures, std::vector<V>& vertices, double (*temperedCost)(V, double), V (*getCandidate)(V), void (*analysis)(const std::vector<const V>&, const std::vector<double>&)=nullptr);

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


#define RECORD 1


//void record(int num, char* msg, FILE* log) {
//	static std::vector<int> line_loc
//	if (vector.size() < num) {
//		line_loc.resize(num+1)
//	}
//}

#endif // PARALLEL_TEMPERING_H
