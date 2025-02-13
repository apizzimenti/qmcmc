#include <vector>
#include <random>
#include <cstdlib>
#include "main.cpp"
#include <omp.h>
// n=7 g=7
// 99    - 605
// 60    - 11
// 1     - .65
// .0001 - .25

// n=9 g=13
// 99    - 1600
// 60    - 26
// 1     - 1?
// .0001 - 

int main(int argc, char const *argv[]) {
    int acceptedNeighbors = 0;
    int energyIncreasingUpdates = 0;
    float temp = .8;
    std::cout << "Temperature: " << temp << std::endl;
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // float randNums[50000];
    // for (int i = 0; i < sizeof(randNums); i++) {
    //     randNums[i] = dist(rng);
    // }

    size_t numPrimary = 9;
    size_t numGates = 13;
    MajTree vertices[50] = {{numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary},
                            {numGates+numPrimary, (int)numPrimary}}; 
    
    // #pragma omp parallel for
    for (int i = 0; i < 50; i++) {
        // printf("I am thread: %d", omp_get_thread_num());
        vertices[i].assignRandomValues(rng);
        vertices[i].outputsInit();
    }

    int sum_boolshit = 0;
    for (int i = 0; i < 200000*2; i++) {
        #pragma omp parallel for reduction (+:energyIncreasingUpdates, acceptedNeighbors, sum_boolshit)
        for (int j = 0; j < 50; j++) {
            // #pragma omp critical
            // {
                sum_boolshit++;
            // }

            MajTree candidate = getNeighbor(vertices[j], rng);
            float oldCost, newCost;
            oldCost = temperedMajCost(vertices[j], temp);
            newCost = temperedMajCost(candidate, temp);
            // #pragma omp critical
            // {
                if (newCost > oldCost) {
                    energyIncreasingUpdates++;
                    if (dist(rng) * newCost < oldCost) {
                        acceptedNeighbors++;
                    }
                }
            // }
            vertices[j] = candidate;
        }
        
    }

    std::cout << "Counter: " << sum_boolshit << std::endl;
    std::cout << "Accepted Neighbors: " << acceptedNeighbors << std::endl;
    std::cout << "Energy Incresing Updates: " << energyIncreasingUpdates << std::endl;
    std::cout << "Acceptance rate of all encountered energy increasing updates: " << 100.0 * acceptedNeighbors/energyIncreasingUpdates << std:: endl;

    return 0;
}
