#include <vector>
#include <iostream>
#include <ctime>
// #include <cstdlib>

struct IsingModel {
    int size;
    std::vector<std::vector<int>> lattice;
    int cost;

    // initializes for the struct function similar to classes
    IsingModel(int n) {
        this->size = n;
        this->lattice = std::vector<std::vector<int>>(n, std::vector<int>(n));
        std::srand(std::time(0));
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                lattice[i][j] = 0;
                if (std::rand() % 2 == 0) {
                    lattice[i][j] = 1;
                }
            }
        }
        calculateCost();
    }

    // ungodly design, but changes the spin of an element and accurately updates the cost
    void updateSpin(int i, int j) {
        lattice[i][j] ^= 1;
        if (((i-1) > -1)) {
            if (lattice[i][j] != lattice[i-1][j]) {
                cost++;
            } else {
                cost--;
            }
        }
        if (((j-1) > -1)) {
            if (lattice[i][j] != lattice[i][j-1]) {
                cost++;
            } else {
                cost--;
            }
        }
        if (((i+1) < size)) {
            if (lattice[i][j] != lattice[i+1][j]) {
                cost++;
            } else {
                cost--;
            }
        }
        if (((j+1) < size)) {
            if (lattice[i][j] != lattice[i][j+1]) {
                cost++;
            } else {
                cost--;
            }
        }
    }

    // exists to calculate cost at the outset, but would work whenever called
    void calculateCost() {
        int cost = 0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (((i+1)<size) && (lattice[i][j] != lattice[i+1][j])) {
                    cost++;
                }
                if (((j+1)<size) && lattice[i][j] != lattice[i][j+1]) {
                    cost++;
                }
            }
        }
        this->cost = cost;
    }

    int getCost() {
        return cost;
    }

    // exists for debugging purposes to be able to visually see the ising model
    void IsingModelPrinter() {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                std::cout << (lattice[i][j] == 1 ? "1 " : "0 ");
            }
            std::cout << std::endl;
        }
    }
};

int main() {
    IsingModel model(3);
    model.IsingModelPrinter();
    std::cout << model.getCost() << std::endl;
    std::cout << std::endl;
    model.updateSpin(0, 0);
    model.updateSpin(1, 0);
    model.updateSpin(0, 2);
    // model.updateSpin(2, 0);
    // model.updateSpin(3, 0);
    // model.updateSpin(4, 0);
    // model.updateSpin(5, 0);
    // model.updateSpin(6, 0);
    // model.updateSpin(7, 0);
    // model.updateSpin(8, 0);
    // model.updateSpin(9, 0);
    model.IsingModelPrinter();
    std::cout << model.getCost() << std::endl;
    return 0;
}
