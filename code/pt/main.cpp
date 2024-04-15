#include "ParallelTempering.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream> // For std::coutging
#include <stdexcept> // For runtime exception

double bouncyFunctionCost(double x) {
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

// n is current stopping condition
template<class V>
V parallelTempering(int n, int m, std::vector<double>& T, V& s, double (*temperedCost)(V, double), V (*getCandidate)(V)) {
    // edge cases
    if (T.length() == 0) {
        return NULL;
    }
    
    // define and initialize walkers
    struct walker {
        double t;
        V v;
    };
    std::vector<walker> W;
    for (auto &&t : T) {
        walkers.push_back({t, s});
    }
    
    // loop until stopping condition is met
    for (int i = 0; i < n; i++) {
        // for each walker
        for (auto &&w: W) {
            // propose m updates
            for (int j = 0; j < m; j++) {
                V u = getCandidate(w.v)
                if (std::uniform_real_distribution<double>(0.0,1.0)(std::random_device{}())
                    < temperedCost(u.v, w.t)/temperedCost(w.v, w.t)) {
                    w.v = u;
                }
            }
        }
        // odd even style sorting temperatures
        for (int j = 0; j < W.length() - 1; j++) {
            if (i % 2 == j % 2) {
                if (std::uniform_real_distribution<double>(0.0,1.0)(std::random_device{}())
                    < temperedCost(W[i+1].v,W[i+1].t)/temperedCost(W[i].v,W[i].t)) {
                    std::swap(W[i], W[i+1])
                }
            }
        }
    }
	
	// calculate and return current min cost vertex
    V minCostV = W.at(0).v;
    for (auto &&w : W) {
        if (temperedCost(w.v) < temperedCost(minCostV)) {
            minCostV = w.v
        }
        else if (temperedCost(w.v) == minCostV) {
            if (std::uniform_real_distribution<double>(0.0,1.0)(std::random_device{}()) < .5) {
                minCostV = w.v
            }
        }
    }
    
    return minCostV;
}


int main() {
	std::vector<double> T = {.25, .5, .75, 1.0, 1.25};
	double startingVal = .6;
	parallelTempering(50, 5, T, startingVal, bouncyFunctionCost, update);
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

