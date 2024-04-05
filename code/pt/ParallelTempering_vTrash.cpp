#include <iostream>
#include <vector>
#include <cmath>
#include <random>

double temperedCost() {
}

generateNeighbor(/* arguments */) {
}

void ParallelTempering(std::vector<>& delta, std::vector<double>& T, int L, std::vector<int>& R) {
    for (int k = 0; k < 100; ++i){
        for (int i = 0; i < T.size(); ++i) {
            for (int j = 0; j < 5; ++j) {
                auto delta_prime = generateNeighbor(delta[i]);
                double p = std::min(1.0, temperedCost(delta[i]) - temperedCost(delta_prime));
                if (/* random number between 0 and 1 < p */) {
                    delta[i] = delta_prime;
                }
            }
        }
        
        for (int i : R) {
            double p = std::min(1.0, exp((1.0 / T[i] - 1.0 / T[i + 1]) * (omega(delta[i]) - omega(delta[i + 1]))));
            if (/* random number between 0 and 1 < p */) {
                std::swap(delta[i], delta[i + 1]);
            }
        }
    }
}

int main() {
    std::vector<> delta;
    std::vector<double> T;
    int L =;
    std::vector<int> R;

    ParallelTempering(delta, T, L, R);
    return 0;
}
