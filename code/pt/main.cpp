#include "ParallelTempering.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <time.h>
#include <stdexcept> // For runtime exception



// Adjust according to actual std::coutic
double temperedCost(double value) {
    // Placeholder
    return value;
}

// Adjust later
void generateNeighbor(double* value) {
    static std::mt19937 generator(std::random_device{}());
    std::uniform_real_distribution<double> distribution(0.0, 5.0);
    
    // Placeholder. Modify according to your actual neighbor generation std::coutic.
    *value = distribution(generator);

}

void ParallelTempering(std::vector<double>& delta, std::vector<double>& T, 
		int L, std::vector<int>& R, std::ofstream& log) {


	DEBUG_CMD(log << "I'm in debug mode!\n";);


    static std::mt19937 rng(std::random_device{}()); // Random number generator
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    for (int k = 0; k < 1000; ++k) {
        for (size_t i = 0; i < T.size(); ++i) {
            for (int j = 0; j < 5; ++j) {
                auto delta_prime = delta[i];
				DEBUG_CMD(log << "Working with " << delta_prime << " ...\n";);

                generateNeighbor(&delta_prime);
				DEBUG_CMD(log << "Proposed update is " << delta_prime << '\n';
					log << "Cost of original is " << temperedCost(delta[i]) << '\n';
					log << "Cost of proposed is " << temperedCost(delta_prime) << '\n';);


                double costDiff = temperedCost(delta[i]) - temperedCost(delta_prime);
				double p = std::exp((1.0 / T[i]) * 
						(temperedCost(delta[i]) - temperedCost(delta_prime)));

				DEBUG_CMD(log << "Probability is " << fmin(1, p) << '\n';);

                if (uni(rng) < p) {
                    delta[i] = delta_prime;
					DEBUG_CMD(log << "Accepted!\n";);
                }
				else {
					DEBUG_CMD(log << "Not Accepted....\n";);
				}

				DEBUG_CMD(log << std::endl;);
            }

			for (size_t i = 0; i < R.size(); ++i) {
				if (i+1 >= R.size()) {break;}

				double p = std::exp((1.0 / T[R[i]]) - (1.0 / T[R[i] + 1]) * 
						(temperedCost(delta[R[i]]) - temperedCost(delta[R[i] + 1])));
				DEBUG_CMD(log << "Proposed swap " << delta[R[i]] 
						<< " and " << delta[R[i] +1] << " with probability " 
						<< fmin(1, p) << '\n';);
				if (uni(rng) < p) {
					std::swap(delta[R[i]], delta[R[i] + 1]);
					DEBUG_CMD(log << "Accepted!\n";);
				}
				else {
					DEBUG_CMD(log << "Not Accepted....\n";);
				}
				DEBUG_CMD(log << std::endl;);
			}
        }
    }
}

int main(int argc, char *argv[]) {
    std::vector<double> delta = {5, 5, 5, 5, 5}; // Input actual data
    std::vector<double> T = {1, 2, 3, 4, 5}; // Temperatures
    std::vector<int> R = {1, 2, 3, 4, 5}; // Indices for swapping

	std::ofstream log;
	log.open("log.txt");
	time_t cur_time = time(NULL);

	log << "Starting logging at " << asctime(localtime(&cur_time)) << '\n';

    ParallelTempering(delta, T, 0, R, log);
    // Example output. Also placeholder. Replace later
    std::cout << "Final states of delta after Parallel Tempering:" << std::endl;
    for (auto val : delta) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    return 0;
}

