#include "ParallelTempering.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream> // For std::coutging
#include <stdexcept> // For runtime exception
#include <omp.h> // For parallel computing
#include "bouncyFunc.h"

// n is the current stopping condition
template<class V>
V parallelTempering(int n, int m, std::vector<double>& T, V& s, double (*temperedCost)(V, double), V (*getCandidate)(V)) {
	// Random Number Generator
	auto gen = std::mt19937(std::random_device{}());
    
    //Define and initialize walkers (states of the system)
    struct walker {
        double t; // Temperature for walker
        V v; // State of the system
    };
    std::vector<walker> W;
	
    // Initialize walkers with each temperature
    for (auto &&t : T) {
        W.push_back({t, s});
    }
	// Output the initial tempered costs for each walker
	for (int i = 0; i < W.size(); i++) {
		double tmp_cost = temperedCost(W[i].v, 0);
		printf("%d:%.3lf\n", i, tmp_cost);
	}

     // Loop until the stopping condition is met
    for (int i = 0; i < n; i++) {
         // Parallel section using OpenMP
        #pragma omp parallel
        {
	// For each walker in the parallelized loop
        for (auto &&w: W) {
            // Propose m updates for the current walker
            for (int j = 0; j < m; j++) {
		// Generate a candidate new state
                V u = getCandidate(w.v);
		// Perform Metropolis acceptance check using the tempered cost
                if (std::uniform_real_distribution<double>(0.0,1.0)(gen)
                    < temperedCost(u, w.t)/temperedCost(w.v, w.t)) {
                    w.v = u; // Accept the new state if the condition is met
                }
		                // Output temporary cost (debug feature)
				double tmp_cost = temperedCost(W[i].v, 0);
				printf("%d:%.3lf\n", j, tmp_cost);
            }
        }
        }
        // odd-even style sorting temperatures
        for (int j = 0; j < W.size() - 1; j++) {
	   // Check if both indices i and j are even or odd
            if (i % 2 == j % 2) {
		 // Generate a uniform random number in the range [0,1]
                if (std::uniform_real_distribution<double>(0.0,1.0)(gen)
                    < temperedCost(W[i+1].v,W[i+1].t)/temperedCost(W[i].v,W[i].t)) {
                    std::swap(W[i], W[i+1]); // Swap adjacent walkers based on their temperatures
                }
                                // Output temporary cost (debug feature)
				double tmp_cost = temperedCost(W[i].v, 0);
				printf("%d:%.3lf\n", j, tmp_cost);
            }
        }
    }
	
// Calculate and return the current vertex (state) with the minimum cost
    V minCostV = W.at(0).v; // Initialize with the first vertex in the walkers list
    for (auto &&w : W) {
	// If the cost of the current walker's vertex is less than the current minimum cost vertex
        if (temperedCost(w.v, 0) < temperedCost(minCostV, 0)) {
	    // Update minCostV to be the vertex with the lower cost
            minCostV = w.v;
        }
	// If the costs are equal, randomly choose between the current and minCostV
        else if (temperedCost(w.v, 0) == minCostV) {
	 // Generate a random number in the range [0,1] and pick the current vertex with 50% probability
            if (std::uniform_real_distribution<double>(0.0,1.0)(gen) < .5) {
                minCostV = w.v;
            }
        }
    }
    
    return minCostV;
}


int main() {
	std::vector<double> T = {.25, .5, .75, 1.0, 1.25}; // Define a vector of temperatures for the tempering process
	double startingVal = .6; // Initialize the starting value for the state

	// Call the parallel tempering function with:
	// - 50 iterations
	// - 5 updates per walker
	// - Temperature vector T
	// - Initial state startingVal
	// - Function calculate_expression as the cost function
	// - Function update to generate candidate states
	auto p = parallelTempering<double>(50, 5, T, startingVal, calculate_expression, update);
	std::cout << p << std::endl; // Output the result of the parallel tempering process
	return 0;
}


// int DEBUG = 1;

// // Adjust according to actual std::coutic
// double temperedCost(double value) {
//     // Placeholder
//     return value;
// }

// // Adjust later
// void generateNeighbor(double* value) {
//     static std::mt19937 generator(std::random_device{}());
//     std::uniform_real_distribution<double> distribution(0.0, 5.0);
    
//     // Placeholder. Modify according to your actual neighbor generation std::coutic.
//     *value = distribution(generator);

// }

// void ParallelTempering(std::vector<double>& delta, std::vector<double>& T, 
// 		int L, std::vector<int>& R) {

// 	// Make sure T, R, and delta have same size
// 	if ( (delta.size() != T.size()) ||
// 			(T.size() != R.size()) ||
// 			(delta.size() != R.size())) {
// 		throw std::runtime_error("Parameters must have same length");
// 	}


// 	if (DEBUG) {
// 		std::cout << "I'm in debug mode!\n";
// 	}


//     static std::mt19937 rng(std::random_device{}()); // Random number generator
//     std::uniform_real_distribution<double> uni(0.0, 1.0);

//     for (int k = 0; k < 1000; ++k) {
//         for (size_t i = 0; i < T.size(); ++i) {
//             for (int j = 0; j < 5; ++j) {
//                 auto delta_prime = delta[i];
// 				if (DEBUG) {
// 					std::cout << "Working with " << delta_prime << "...\n";
// 				}

//                 generateNeighbor(&delta_prime);
// 				if (DEBUG) {
// 					std::cout << "Proposed update is " << delta_prime << "\n";
// 					std::cout << "Cost of original is " << temperedCost(delta[i]) << '\n';
// 					std::cout << "Cost of proposed is " << temperedCost(delta_prime) << '\n';
// 				}


//                 double costDiff = temperedCost(delta[i]) - temperedCost(delta_prime);
// 				if (DEBUG) {
// 					std::cout << "Cost difference is " << costDiff << '\n';
// 				}
// 				double p = std::exp((1.0 / T[i]) * 
// 						(temperedCost(delta[R[i]]) - temperedCost(delta_prime)));
//                 if (uni(rng) < p) {
//                     delta[i] = delta_prime;
// 					if (DEBUG) {
// 						std::cout << "Accepted!\n";
// 					}
//                 }
// 				else {
// 					if (DEBUG) {
// 						std::cout << "Not Accepted....\n";
// 					}
// 				}

// 				if (DEBUG) {
// 					std::cout << std::endl;
// 				}
//             }

// 			for (size_t i = 0; i < R.size(); ++i) {
// 				if (i+1 >= R.size()) {break;}

// 				double p = std::exp((1.0 / T[R[i]]) - (1.0 / T[R[i] + 1]) * 
// 						(temperedCost(delta[R[i]]) - temperedCost(delta[R[i] + 1])));
// 				if (DEBUG) {
// 					std::cout << "Proposed swap " << delta[R[i]] 
// 						<< " and " << delta[R[i] +1] << " with probability "
// 						<< p << '\n';
// 				}
// 				if (uni(rng) < p) {
// 					std::swap(delta[R[i]], delta[R[i] + 1]);
// 					if (DEBUG) {
// 						std::cout << "Accepted!\n";
// 					}
// 				}
// 				else {
// 					if (DEBUG) {
// 						std::cout << "Not Accepted....\n";
// 					}
// 				}
// 				std::cout << std::endl;
// 			}
//         }
//     }
// }

// int main() {
//     std::vector<double> delta = {5, 5, 5, 5, 5}; // Input actual data
//     std::vector<double> T = {1, 2, 3, 4, 5}; // Temperatures
//     std::vector<int> R = {1, 2, 3, 4, 5}; // Indices for swapping
// 	std::ofstream log("log.txt");
// 	log << "poop\n";
// 	log.close();

//     ParallelTempering(delta, T, 0, R);
//     // Example output. Also placeholder. Replace later
//     std::cout << "Final states of delta after Parallel Tempering:" << std::endl;
//     for (auto val : delta) {
//         std::cout << val << " ";
//     }
//     std::cout << std::endl;

//     return 0;
// }
