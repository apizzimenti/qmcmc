#include "ParallelTempering.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <omp.h>
#include "bouncyFunc.h"

/**
 * @brief Runs a parallel tempering optimization algorithm to search for low-cost configurations of vertices.
 * 
 * @tparam V The type of vertices being optimized.
 * @param n The number of iterations to run the algorithm for.
 * @param m The number of proposals or candidate configurations to generate for each vertex per iteration.
 * @param temperatures A vector containing the temperature values for each walker (used to adjust acceptance probabilities).
 * @param vertices A vector of initial vertex configurations, one for each walker.
 * @param temperedCost A function pointer that takes a vertex and a temperature, returning a cost value.
 * @param getCandidate A function pointer to generate candidate configurations from a given vertex.
 * @param analysis (Optional) A function pointer for any analysis to perform at the end of each iteration, taking vertices and temperatures as inputs.
 * @return A vector of vertices representing the final configuration of each walker after the specified number of iterations.
 * 
 * @details This function implements the parallel tempering (replica exchange) algorithm. Each "walker" searches
 * in parallel for a low-cost configuration of a vertex. The temperature values are used to modify the
 * acceptance probability for each proposed candidate configuration, allowing the algorithm to explore a range
 * of states that balance exploration and exploitation.
 */
template<class V>
std::vector<V> parallelTempering(int n, int m, std::vector<double>& temperatures, std::vector<V>& vertices, double (*temperedCost)(V, double), V (*getCandidate)(V), void (*analysis)(const std::vector<const V>&, const std::vector<double>&)=nullptr) {
    if (temperatures.size() != vertices.size()) {
        std::cout << "Temperature vector size and Vertices vector size do not match" << std::endl;
        return nullptr;
    }
    int num_walkers = temperatures.size();

    // Initialize vertices array
    // std::vector<V> vertices(num_walkers, s);
    

    // Debug: Print initial cost of each vertex
    // for (int i = 0; i < num_walkers; i++) {
    //     double tmp_cost = temperedCost(vertices[i], 0);
    //     printf("%d:%.3lf\n", i, tmp_cost);
    // }

    // Loop until stopping condition is met
    for (int i = 0; i < n; i++) {
        #pragma omp parallel for
        for (int j = 0; j < num_walkers; j++) {
            V& vertex = vertices[j];
            double temperature = temperatures[j];

            // Thread-specific random generator
            std::mt19937 gen(std::random_device{}() + omp_get_thread_num());

            // Propose m updates for each vertex
            for (int k = 0; k < m; k++) {
                V candidate = getCandidate(vertex);
                //TODO
                if (std::uniform_real_distribution<double>(0.0, 1.0)(gen)
                    < temperedCost(candidate, temperature) / temperedCost(vertex, temperature)) {
                    vertex = candidate;
                }
            }
        }

        // Odd-even style sorting of vertices only
        int parity = i % 2;
        #pragma omp parallel for
        for (int j = parity; j < num_walkers - 1; j+=2) {
            std::mt19937 gen(std::random_device{}() + omp_get_thread_num());

            double cost1 = temperedCost(vertices[j], temperatures[j]);
            double cost2 = temperedCost(vertices[j + 1], temperatures[j + 1]);
            if (std::uniform_real_distribution<double>(0.0, 1.0)(gen) < cost2 / cost1) {
                std::swap(vertices[j], vertices[j + 1]);
            }
        }
        if (analysis != nullptr) {
            analysis(vertices, temperatures);
        }
    }

    return vertices;

    // Find and return the vertex with the minimum cost
    // V minCostVertex = vertices[0];
    // double minCost = temperedCost(minCostVertex, 0);

    // for (int i = 1; i < num_walkers; i++) {
    //     double cost = temperedCost(vertices[i], 0);
    //     if (cost < minCost) {
    //         minCostVertex = vertices[i];
    //         minCost = cost;
    //     }
    // }

    // return minCostVertex;
}

int main() {
    std::vector<double> T = {0.25, 0.5, 0.75, 1.0, 1.25};
    double startingVal = 0.6;
    // auto p = parallelTempering<double>(50, 5, T, startingVal, calculate_expression, update);
    // std::cout << p << std::endl;
    return 0;
}




// #include "ParallelTempering.hpp"
// #include <iostream>
// #include <vector>
// #include <cmath>
// #include <random>
// #include <fstream> // For std::coutging
// #include <stdexcept> // For runtime exception
// #include <omp.h>
// #include "bouncyFunc.h"

// // n is current stopping condition
// template<class V>
// V parallelTempering(int n, int m, std::vector<double>& T, V& s, double (*temperedCost)(V, double), V (*getCandidate)(V)) {
// 	// Make an rng
// 	auto gen = std::mt19937(std::random_device{}());

//     edge cases
//     if (T.size() == 0) {
//         return nullptr;
//     }
    
//     // define and initialize walkers
//     struct walker {
//         double t;
//         V v;
//     };
//     std::vector<walker> W;
//     for (auto &&t : T) {
//         W.push_back({t, s});
//     }
// 	// for (int i = 0; i < W.size(); i++) {
// 	// 	double tmp_cost = temperedCost(W[i].v, 0);
// 	// 	printf("%d:%.3lf\n", i, tmp_cost);
// 	// }

    
//     // loop until stopping condition is met

//     for (int i = 0; i < n; i++) {
//         // for each walker
//         #pragma omp parallel for
//         for (auto &&w: W) {
//             // propose m updates
//             for (int j = 0; j < m; j++) {
//                 V u = getCandidate(w.v);
//                 if (std::uniform_real_distribution<double>(0.0,1.0)(gen)
//                     < temperedCost(u, w.t)/temperedCost(w.v, w.t)) {
//                     w.v = u;
//                 }
// 				double tmp_cost = temperedCost(W[i].v, 0);
// 				printf("%d:%.3lf\n", j, tmp_cost);
//             }
//         }
//         // odd even style sorting temperatures
//         #pragma omp parallel for
//         for (int j = 0; j < W.size() - 1; j++) {
//             if (i % 2 == j % 2) {
//                 if (std::uniform_real_distribution<double>(0.0,1.0)(gen)
//                     < temperedCost(W[j+1].v,W[j+1].t)/temperedCost(W[j].v,W[j].t)) {
//                     std::swap(W[j], W[j+1]);
//                 }

// 				double tmp_cost = temperedCost(W[i].v, 0);
// 				printf("%d:%.3lf\n", j, tmp_cost);
//             }
//         }
//     }
	
// 	// calculate and return current min cost vertex
//     V minCostV = W.at(0).v;
//     for (auto &&w : W) {
//         if (temperedCost(w.v, 0) < temperedCost(minCostV, 0)) {
//             minCostV = w.v;
//         }
//         else if (temperedCost(w.v, 0) == minCostV) {
//             if (std::uniform_real_distribution<double>(0.0,1.0)(gen) < .5) {
//                 minCostV = w.v;
//             }
//         }
//     }
    
//     return minCostV;
// }


// int main() {
// 	std::vector<double> T = {.25, .5, .75, 1.0, 1.25};
// 	double startingVal = .6;
// 	auto p = parallelTempering<double>(50, 5, T, startingVal, calculate_expression, update);
// 	std::cout << p << std::endl;
// 	return 0;
// }


// // int DEBUG = 1;

// // // Adjust according to actual std::coutic
// // double temperedCost(double value) {
// //     // Placeholder
// //     return value;
// // }

// // // Adjust later
// // void generateNeighbor(double* value) {
// //     static std::mt19937 generator(std::random_device{}());
// //     std::uniform_real_distribution<double> distribution(0.0, 5.0);
    
// //     // Placeholder. Modify according to your actual neighbor generation std::coutic.
// //     *value = distribution(generator);

// // }

// // void ParallelTempering(std::vector<double>& delta, std::vector<double>& T, 
// // 		int L, std::vector<int>& R) {

// // 	// Make sure T, R, and delta have same size
// // 	if ( (delta.size() != T.size()) ||
// // 			(T.size() != R.size()) ||
// // 			(delta.size() != R.size())) {
// // 		throw std::runtime_error("Parameters must have same length");
// // 	}


// // 	if (DEBUG) {
// // 		std::cout << "I'm in debug mode!\n";
// // 	}


// //     static std::mt19937 rng(std::random_device{}()); // Random number generator
// //     std::uniform_real_distribution<double> uni(0.0, 1.0);

// //     for (int k = 0; k < 1000; ++k) {
// //         for (size_t i = 0; i < T.size(); ++i) {
// //             for (int j = 0; j < 5; ++j) {
// //                 auto delta_prime = delta[i];
// // 				if (DEBUG) {
// // 					std::cout << "Working with " << delta_prime << "...\n";
// // 				}

// //                 generateNeighbor(&delta_prime);
// // 				if (DEBUG) {
// // 					std::cout << "Proposed update is " << delta_prime << "\n";
// // 					std::cout << "Cost of original is " << temperedCost(delta[i]) << '\n';
// // 					std::cout << "Cost of proposed is " << temperedCost(delta_prime) << '\n';
// // 				}


// //                 double costDiff = temperedCost(delta[i]) - temperedCost(delta_prime);
// // 				if (DEBUG) {
// // 					std::cout << "Cost difference is " << costDiff << '\n';
// // 				}
// // 				double p = std::exp((1.0 / T[i]) * 
// // 						(temperedCost(delta[R[i]]) - temperedCost(delta_prime)));
// //                 if (uni(rng) < p) {
// //                     delta[i] = delta_prime;
// // 					if (DEBUG) {
// // 						std::cout << "Accepted!\n";
// // 					}
// //                 }
// // 				else {
// // 					if (DEBUG) {
// // 						std::cout << "Not Accepted....\n";
// // 					}
// // 				}

// // 				if (DEBUG) {
// // 					std::cout << std::endl;
// // 				}
// //             }

// // 			for (size_t i = 0; i < R.size(); ++i) {
// // 				if (i+1 >= R.size()) {break;}

// // 				double p = std::exp((1.0 / T[R[i]]) - (1.0 / T[R[i] + 1]) * 
// // 						(temperedCost(delta[R[i]]) - temperedCost(delta[R[i] + 1])));
// // 				if (DEBUG) {
// // 					std::cout << "Proposed swap " << delta[R[i]] 
// // 						<< " and " << delta[R[i] +1] << " with probability "
// // 						<< p << '\n';
// // 				}
// // 				if (uni(rng) < p) {
// // 					std::swap(delta[R[i]], delta[R[i] + 1]);
// // 					if (DEBUG) {
// // 						std::cout << "Accepted!\n";
// // 					}
// // 				}
// // 				else {
// // 					if (DEBUG) {
// // 						std::cout << "Not Accepted....\n";
// // 					}
// // 				}
// // 				std::cout << std::endl;
// // 			}
// //         }
// //     }
// // }

// // int main() {
// //     std::vector<double> delta = {5, 5, 5, 5, 5}; // Input actual data
// //     std::vector<double> T = {1, 2, 3, 4, 5}; // Temperatures
// //     std::vector<int> R = {1, 2, 3, 4, 5}; // Indices for swapping
// // 	std::ofstream log("log.txt");
// // 	log << "poop\n";
// // 	log.close();

// //     ParallelTempering(delta, T, 0, R);
// //     // Example output. Also placeholder. Replace later
// //     std::cout << "Final states of delta after Parallel Tempering:" << std::endl;
// //     for (auto val : delta) {
// //         std::cout << val << " ";
// //     }
// //     std::cout << std::endl;

// //     return 0;
// // }
