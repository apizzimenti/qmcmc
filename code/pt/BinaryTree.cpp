#include "BinaryTree.hpp"
#include <cstdlib>
#include <ctime>
#include <stack>

BinaryTree::BinaryTree(int depth) : maxDepth(depth) {
    root = new TreeNode(0, 1.0); // Root node at depth 0 with initial cost
    AssignCosts(root, 0);        // Assign costs to all nodes
    PerturbPath(root);           // Perturb the path to the target leaf
}

BinaryTree::~BinaryTree() {
    DeleteTree(root);
}

void BinaryTree::AssignCosts(TreeNode* node, int currentDepth) {
    if (currentDepth >= maxDepth) return;

    // Create left and right children
    node->left = new TreeNode(currentDepth + 1, std::pow(2, -(currentDepth + 1)));
    node->right = new TreeNode(currentDepth + 1, std::pow(2, -(currentDepth + 1)));

    // Recursively assign costs to the child nodes
    AssignCosts(node->left, currentDepth + 1);
    AssignCosts(node->right, currentDepth + 1);
}

void BinaryTree::PerturbPath(TreeNode* node) {
    TreeNode* current = node;
    srand(static_cast<unsigned>(time(0)));

    // Randomly select a path to the leaf and perturb the costs
    while (current->left && current->right) {
        perturbedPathNodes.push_back(current);
        current->cost *= 0.5; // Perturb the cost
        if (rand() % 2 == 0) {
            current = current->left;
        } else {
            current = current->right;
        }
    }
    perturbedPathNodes.push_back(current); // Add the final leaf
    current->cost *= 0.5; // Perturb the leaf cost
}

void BinaryTree::SimulateWalkers(int numWalkers, int steps, const std::string& filename) {
    std::ofstream outfile(filename);

    // Each walker starts at the root
    std::vector<TreeNode*> walkers(numWalkers, root);
    
    // Simulate walker movements
    for (int step = 0; step < steps; ++step) {
        for (int i = 0; i < numWalkers; ++i) {
            if (walkers[i]->left && walkers[i]->right) {
                // Randomly move to the left or right child
                if (rand() % 2 == 0) {
                    walkers[i] = walkers[i]->left;
                } else {
                    walkers[i] = walkers[i]->right;
                }
            }
            // Output walker's current position (depth and cost)
            outfile << step << " " << i << " " << walkers[i]->depth << " " << walkers[i]->cost << "\n";
        }
    }

    outfile.close();
}

void BinaryTree::DeleteTree(TreeNode* node) {
    if (!node) return;
    DeleteTree(node->left);
    DeleteTree(node->right);
    delete node;
}

void BinaryTree::PrintTree() {
    PrintTree(root, 0);
}

void BinaryTree::PrintTree(TreeNode* node, int level) {
    if (!node) return;
    
    PrintTree(node->right, level + 1);
    
    for (int i = 0; i < level; ++i) {
        std::cout << "    ";
    }
    
    std::cout << "(" << node->cost << ")\n";
    
    PrintTree(node->left, level + 1);
}
