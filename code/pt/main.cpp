#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <random>
#include <omp.h>
#include <cstdlib>
#include <ctime>
#include <queue>

#define ROWSIZE 3
#define LEFT 0
#define MIDDLE 1
#define RIGHT 2
#define OUTPUT 3
#define UNKNOWN -3

#define CONST_0 -2
#define CONST_1 -1

#define MAJ3(a, b, c) (((a) + (b) + (c)) > 1 ? 1 : 0)

#define N 9
#define G 13
#define II 10000
#define END_EARLY 0

int memoized_majorities[1 << N];
int memoinit = 0;

int actualMaj(int num, int bits) {
    if (!memoinit) {
        for (int i = 0; i < 1 << N; i++) {
            memoized_majorities[i] = -1;
        }
        memoinit = 1;
    }
    if (memoized_majorities[num] != -1) {
        return memoized_majorities[num];
    } else {
        int counter1s = 0;
        int numCopy = num;
        for (int i = 0; i < bits; i++) {
            if (numCopy & 1) {
                counter1s++;
            }
            numCopy >>= 1;
        }
        if (counter1s > bits/2) {
            return (memoized_majorities[num] = 1);
        } else {
            return (memoized_majorities[num] = 0);
        }
    }
}

class MajTree {
public:
    std::shared_ptr<int>* data;  // Public 2D array
    size_t gates; // Number of gates
    int n;       // Additional integer member
    std::shared_ptr<char>* outputs;
    int cost;

    // Constructor to initialize the array and set the value of n
    MajTree(size_t rowCount, int nValue) : gates(rowCount), n(nValue) {
        data = new std::shared_ptr<int>[gates];
        for (size_t i = 0; i < gates; ++i) {
            data[i] = std::shared_ptr<int> (new int[ROWSIZE], [](int* p) { delete[] p; });
        }

        outputs = new std::shared_ptr<char>[gates+1];
        for (size_t i = 0; i < gates+1; ++i) {
            outputs[i] = std::shared_ptr<char> (new char[1 << n], [](char* p) { delete[] p; });
        }
    }

    // Default constructor
    MajTree() : MajTree(0, 0) {}

    // Copy Constructor
    MajTree(const MajTree& other) {
        n = other.n;
        gates = other.gates;
        cost = other.cost;
        data = new std::shared_ptr<int>[gates];
        for (size_t i = 0; i < gates; ++i) {
            data[i] = other.data[i];
            // for (int j = 0; j < ROWSIZE; j++) {
            //     data[i].get()[j] = other.data[i].get()[j];
            // }
        }

        outputs = new std::shared_ptr<char>[gates+1];
        for (size_t i = 0; i < gates+1; ++i) {
            outputs[i] = other.outputs[i];
            // for (size_t j = 0; j < (1 << n); j++) {
            //     outputs[i][j] = other.outputs[i][j];
            // }
        }
    }

    MajTree& operator=(const MajTree& other) {
        if (this == &other) {
            return *this;
        }

        delete[] data;

        delete[] outputs;
        
        n = other.n;
        gates = other.gates;
        cost = other.cost;
        data = new std::shared_ptr<int>[gates];
        for (size_t i = 0; i < gates; ++i) {
            data[i] = other.data[i];
            // for (int j = 0; j < ROWSIZE; j++) {
            //     data[i].get()[j] = other.data[i].get()[j];
            // }
        }

        outputs = new std::shared_ptr<char>[gates+1];
        for (size_t i = 0; i < gates+1; ++i) {
            outputs[i] = other.outputs[i];
            // for (size_t j = 0; j < (1 << n); j++) {
            //     outputs[i][j] = other.outputs[i][j];
            // }
        }

        return *this;
    }

    // Destructor to clean up dynamically allocated memory
    ~MajTree() {
        delete[] data;

        delete[] outputs;
    }

    // Method to randomly assign values to columns 0-2
    // void assignRandomValues() {
    //     std::srand(static_cast<unsigned>(std::time(0)));
    //     for (int gate = n; gate < gates; ++gate) {
    //         for (int child = 0; child < (ROWSIZE - 1); ++child) {
    //             data[gate].get()[child] = (std::rand() % (gate)); // Random integer [-2, gate)
    //         }
    //         while ((data[gate].get()[0] == data[gate].get()[1]) || (data[gate].get()[0] == data[gate].get()[2])) {
    //             data[gate].get()[0] = (std::rand() % (gate));
    //         }
    //         // worst case: <x,y,y>
    //         while ((data[gate].get()[1] == data[gate].get()[0]) || (data[gate].get()[1] == data[gate].get()[2])) {
    //             data[gate].get()[1] = (std::rand() % (gate));
    //         }
    //         // worst case: <x,z,y>
    //     }
    //     // TODO
    //     for (int i = 0; i < gates; i++) {
    //         data[i].get()[3]=UNKNOWN;
    //     }
        
    // }
    void assignRandomValues(std::mt19937& rng) {
        std::uniform_int_distribution<int> dist;
        for (int gate = n; gate < gates; ++gate) {
            for (int child = 0; child < (ROWSIZE); ++child) {
                data[gate].get()[child] = (dist(rng) % (gate + 2)) - 2; // Random integer [-2, gate)
            }
            while ((data[gate].get()[0] == data[gate].get()[1]) || (data[gate].get()[0] == data[gate].get()[2])) {
                data[gate].get()[0] = (dist(rng) % (gate + 2)) -2;
            }
            // worst case: <x,y,y>
            while ((data[gate].get()[1] == data[gate].get()[0]) || (data[gate].get()[1] == data[gate].get()[2])) {
                data[gate].get()[1] = (dist(rng) % (gate + 2)) -2;
            }
            // worst case: <x,z,y>
        }
    }

    void outputsInit() {
        cost = 0;
        for (int i = 0; i < (1 << n); i++) {
            int iCopy = i;
            for (int j = 0; j < n; j++) {
                outputs[j].get()[i] = iCopy & 1;
                iCopy >>= 1;
            }
            for (int j = n; j < gates; j++) {
                char a = data[j].get()[LEFT];
                char b = data[j].get()[MIDDLE];
                char c = data[j].get()[RIGHT];
                a = (a<0) ? a + 2 : outputs[a].get()[i];
                b = (b<0) ? b + 2 : outputs[b].get()[i];
                c = (c<0) ? c + 2 : outputs[c].get()[i];
                outputs[j].get()[i] = MAJ3(a,b,c);
            }
            if (actualMaj(i, n) != outputs[gates-1].get()[i]) {
                outputs[gates].get()[i] = 0;
                cost++;
            } else {
                outputs[gates].get()[i] = 1;
            }
        }
    }

    // Method to print the array
    void print() const {
        std::cout << "Array contents:" << std::endl;
        for (size_t i = n; i < gates; ++i) {
            std::cout << "x_" << i << " = < ";
            for (size_t j = 0; j < ROWSIZE; ++j) {
                std::cout << data[i].get()[j] << " ";
            }
            std::cout << ">" << std::endl;
        }
        std::cout << "Value of n: " << n << std::endl;
    }
};

// POOP
int majCost2(MajTree& array) {
    int cost = 0;
    int outs[1 << array.n];
    int exTree[array.gates];
    for (int i = 0; i < 1 << array.n; i++) {
        int iCopy = i;
        for (int j = 0; j < array.n; j++) {
            exTree[j] = iCopy & 1;
            iCopy >>= 1;
        }
        for (int j = array.n; j < array.gates; j++) {
            char a = array.data[j].get()[LEFT];
            char b = array.data[j].get()[MIDDLE];
            char c = array.data[j].get()[RIGHT];
            a = (a<0) ? a + 2 : exTree[a];
            b = (b<0) ? b + 2 : exTree[b];
            c = (c<0) ? c + 2 : exTree[c];
            exTree[j] = MAJ3(a,b,c);
        }
        outs[i] = exTree[array.gates-1];
    }

    for (int i = 0; i < (1 << array.n); i++) {
        if (outs[i] != actualMaj(i, array.n)) {
            cost++;
        }
    }
    return cost;
}
// // END POOP

int majCost(MajTree& array) {
    // if (array.cost != majCost2(array)) {
    //     std::cout << "I'M GONNA DIE!!" << std::endl;
    // }
    return array.cost;
}

float temperedMajCost(MajTree& array, float temp) {
    int untemperedCost = majCost(array);
    float beta = 1.0/temp;
    if (untemperedCost == 0 && END_EARLY) {
        std::cout << "Done:" << std::endl;
        array.print();
        std::cout << "DOUBLE CHECKING COST: " << majCost2(array) << std::endl;
        exit(0);
        return 0;
    } else {
        return exp(beta * untemperedCost);
    }
}

MajTree getNeighbor(MajTree& array, std::mt19937& rng) {
    // MajTree array2(array.getSize(), array.getN());
    MajTree array2(array.gates, array.n);
    // for (int i = 0; i < array2.gates; i++) {
    //     for (int j = 0; j < 4; j++) {
    //         array2.data[i][j] = array.data[i][j];
    //     }
    // }
    array2 = array;
    
    std::uniform_int_distribution<int> distNode(-2, array2.gates - 2);
    std::uniform_int_distribution<int> distGate(array2.n, array2.gates - 1);
    std::uniform_int_distribution<int> distInput(0, 2);

    int randNode, randGate, randInput;
    do {
        randNode = distNode(rng);
        randGate = distGate(rng);
        randInput = distInput(rng);
    } while ((randNode == array2.data[randGate].get()[(randInput + 1) % 3]) ||
             (randNode == array2.data[randGate].get()[(randInput + 2) % 3]) ||
              randNode >= randGate ||
              randNode == array2.data[randGate].get()[randInput]);
    

    // int size2 = array2.gates;
    // int randomNode = (std::rand() % (size2 + 2)) - 2;
    // int randGate = (std::rand() % (size2-array2.n)) + array2.n;
    // int randInput = std::rand() % 3;

    // while ((randomNode == array2.data[randGate][(randInput + 1) % 3]) || (randomNode == array2.data[randGate][(randInput + 2) % 3]) || randomNode >= randGate) {
    //     randomNode = (std::rand() % (size2 + 2)) - 2;
    //     randGate = (std::rand() % (size2-array2.n)) + array2.n;
    //     randInput = std::rand() % 3;
    // }

    array2.data[randGate] = std::shared_ptr<int> (new int[ROWSIZE], [](int* p) { delete[] p; });
    array2.data[randGate].get()[randInput] = randNode;
    array2.data[randGate].get()[(randInput + 1) % 3] = array.data[randGate].get()[(randInput + 1) % 3];
    array2.data[randGate].get()[(randInput + 2) % 3] = array.data[randGate].get()[(randInput + 2) % 3];

    char flags[array2.gates] = {0};
    flags[randGate] = 1;
    char flagsDone = 0;
    char minGate = randGate;
    while (!flagsDone) {
        for (char gate = minGate+1; gate < array2.gates; gate++) {
            if (flags[gate] != 1 &&
                (array2.data[gate].get()[LEFT] == minGate ||
                array2.data[gate].get()[MIDDLE] == minGate ||
                array2.data[gate].get()[RIGHT] == minGate)) {
                    flags[gate] = 1;
                }
        }
        flagsDone = 1;
        for (char gate = minGate+1; gate < array2.gates; gate++) {
            if (flags[gate] == 1) {
                minGate = gate;
                flagsDone = 0;
                break;
            }
        }
    }

    for (char j = randGate; j < array2.gates; j++) {
        if (flags[j] == 1) {
            array2.outputs[j] = std::shared_ptr<char> (new char[1 << array2.n], [](char* p) { delete[] p; });
            for (int i = 0; i < (1 << array2.n); i++) {
                char a = array2.data[j].get()[LEFT];
                char b = array2.data[j].get()[MIDDLE];
                char c = array2.data[j].get()[RIGHT];
                a = (a<0) ? a + 2 : array2.outputs[a].get()[i];
                b = (b<0) ? b + 2 : array2.outputs[b].get()[i];
                c = (c<0) ? c + 2 : array2.outputs[c].get()[i];
                array2.outputs[j].get()[i] = MAJ3(a,b,c);
            }
        }
    }
    array2.outputs[array2.gates] = std::shared_ptr<char> (new char[1 << array2.n], [](char* p) { delete[] p; });
    for (int i = 0; i < (1 << array2.n); i++) {
        if (actualMaj(i, array2.n) != array2.outputs[array2.gates-1].get()[i]) {
            array2.outputs[array2.gates].get()[i] = 0;
            if (array.outputs[array2.gates].get()[i] != 0) {
                array2.cost++;
            }
        } else {
            array2.outputs[array2.gates].get()[i] = 1;
            if (array.outputs[array2.gates].get()[i] != 1) {
                array2.cost--;
            }
        }
    }
    
    // for (int i = 0; i < (1 << array2.n); i++) {
    //     for (int j = randGate; j < array2.gates; j++) {
    //         char a = array2.data[j].get()[LEFT];
    //         char b = array2.data[j].get()[MIDDLE];
    //         char c = array2.data[j].get()[RIGHT];
    //         a = (a<0) ? a + 2 : array2.outputs[a].get()[i];
    //         b = (b<0) ? b + 2 : array2.outputs[b].get()[i];
    //         c = (c<0) ? c + 2 : array2.outputs[c].get()[i];
    //         array2.outputs[j].get()[i] = MAJ3(a,b,c);
    //     }
    //     if (actualMaj(i, array2.n) != array2.outputs[array2.gates-1].get()[i]) {
    //         if (array2.outputs[array2.gates].get()[i] != 0) {
    //             array2.outputs[array2.gates].get()[i] = 0;
    //             array2.cost++;
    //         }
    //     } else {
    //         if (array2.outputs[array2.gates].get()[i] != 1) {
    //             array2.outputs[array2.gates].get()[i] = 1;
    //             array2.cost--;
    //         }
    //     }
    // }
    

    return array2;
}

/**
 * @brief Runs a parallel tempering optimization algorithm to search for low-cost configurations of vertices.
 * 
 * @tparam V The type of vertices being optimized.
 * @param numSwaps The number of iterations to run the algorithm for.
 * @param numUpdates The number of proposals or candidate configurations to generate for each vertex per iteration.
 * @param arrSize The number of vertices/temperatures.
 * @param temperatures An array containing the temperature values for each walker (used to adjust acceptance probabilities).
 * @param vertices An array of initial vertex configurations, one for each walker.
 * @param temperedCost A function pointer that takes a vertex and a temperature, returning a cost value.
 * @param getCandidate A function pointer to generate candidate configurations from a given vertex.
 * @param analysis (Optional) A function pointer for any analysis to perform at the end of each iteration, taking vertices and temperatures as inputs.
 * @return An array of vertices representing the final configuration of each walker after the specified number of iterations.
 * 
 * @details This function implements the parallel tempering (replica exchange) algorithm. Each "walker" searches
 * in parallel for a low-cost configuration of a vertex. The temperature values are used to modify the
 * acceptance probability for each proposed candidate configuration, allowing the algorithm to explore a range
 * of states that balance exploration and exploitation.
 */
template<class V>
V* parallelTempering(int numSwaps, int numUpdates, int arrSize, float temperatures[], V vertices[], float (*temperedCost)(V&, float), V (*getNeighbor)(V&, std::mt19937&), void (*analysis)(V[], float[])) {
    std::random_device rd;
    std::mt19937 rng(rd());
    
    int rand1size,rand2size;
    float randNums1[rand1size = (arrSize * numUpdates)];
    float randNums2[rand2size = (arrSize / 2)];
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < numSwaps; i++) {
    // for (int i = 0; i == i; i++) {
        for (int j = 0; j < rand1size; j++) {
            randNums1[j] = dist(gen);
        }
        
        #pragma omp parallel for
        for (int j = 0; j < arrSize; j++) {
            V& vertex = vertices[j];
            float temperature = temperatures[j];
            // Propose m updates for each vertex
            for (int k = 0; k < numUpdates; k++) {
                V candidate = getNeighbor(vertex, rng);
                if (randNums1[j] * temperedCost(candidate, temperature) < temperedCost(vertex, temperature)) {
                    vertex = candidate;
                }
            }
        }

        // Odd-even style sorting of vertices only
        for (int j = 0; j < rand2size; j++) {
            randNums2[j] = dist(gen);
        }
        int parity = i % 2;
        for (int j = parity, k = 0; j < arrSize - 1; j+=2, k++) {
            if ((temperedCost(vertices[j], temperatures[j]) * randNums2[k]) < temperedCost(vertices[j+1], temperatures[j+1])) {
                std::swap(vertices[j], vertices[j + 1]);
            }
        }
        if (analysis != nullptr) {
            analysis(vertices, temperatures);
        }
        // if (i % 10 == 0) {
        //     std::cout << i << ":" << std::endl;
        //     for (int i = 0; i < 5; i++) {
        //         vertices[i].print();
        //     }
        // }
        
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

// int main() {
//     std::random_device rd;
//     std::mt19937 rng(rd());


//     // float temps[] = {0.25, 0.5, 0.75, 1.0, 1.25};
//     float temps[50];
//     float t1 = 605, t2 = 11, tk = .65, tl=.25;
//     // double startingVal = 0.6;
//     // n=7 g=7
//     // 99    - 605
//     // 60    - 11
//     // 1     - .65
//     // .0001 - .25
//     temps[0] = t1;
//     temps[1] = t2;
//     float t2k = (t2-tk) / (10-2);
//     for (int i = 2; i < 10; i++) {
//         temps[i] = temps[i-1] - t2k;
//     }
//     float tkl = (tk-tl) / (50-10);
//     for (int i = 10; i < 50; i++) {
//         temps[i] = temps[i-1] - tkl;

//     }

//     // for (int i = 0; i < 50; i++) {
//     //     std::cout << i << ": " << temps[i] << std::endl;
//     // }
    
//     size_t numPrimary = N;
//     size_t numGates = G;
//     MajTree vertices[50] = {{numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary},
//                             {numGates+numPrimary, (int)numPrimary}}; 
//     for (int i = 0; i < 50; i++) {
//         // vertices[i] = new MajTree(numGates+numPrimary, numPrimary);
//         vertices[i].assignRandomValues(rng);
//         vertices[i].outputsInit();
//         // if (vertices[i].cost != majCost2(vertices[i])) {
//         //     std::cout << "STRAIGHT UP KYS BRUH" << std::endl;
//         // }
//     }

//     int cost;
//     vertices[0].print();
//     cost = majCost(vertices[0]);
//     std::cout << "Calculated Original Cost: " << cost << std::endl;

//     // // BEGIN getNeighbor test
//     // int i = 0;
//     // int oldCost = cost, newCost = cost;
//     // vertices[1] = vertices[0];
//     // int j = 0;
    
//     // while (8 > j++) {
//     //     while (newCost >= oldCost) {
//     //         vertices[1] = getNeighbor(vertices[1], rng);
//     //         newCost = majCost(vertices[1]);
//     //         i++;
//     //     }
//     //     oldCost = newCost;
//     //     vertices[1].print();
//     //     std::cout << "Calculated Updated Cost: " << newCost << std::endl;
//     //     std::cout << "Number of steps taken: " << i << std::endl;
//     // }
//     // // END getNeighbor test
    
//     // vertices[1].print();
//     // std::cout << "Calculated Updated Cost: " << newCost << std::endl;
//     // std::cout << "Number of steps taken: " << i << std::endl;
    


//     MajTree* newVertices = parallelTempering<MajTree>(II, 5, 50, temps, vertices, temperedMajCost, getNeighbor, nullptr);
    
//     int minCost = majCost(newVertices[0]);
//     int mini = 0;
//     for (int i = 1; i < 5; i++) {
//         int currCost = majCost(newVertices[i]);
//         if (currCost < minCost) {
//             minCost = currCost;
//             mini = i;
//         }
//     }
//     std::cout << "bruh?" << std::endl;
//     (newVertices[mini]).print();
//     std::cout << "Min Calculated Cost: " << majCost(newVertices[mini]) << std::endl;

//     return 0;
// }
