#ifndef BINARY_TREE_HPP
#define BINARY_TREE_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

// Node structure for the binary tree
struct TreeNode {
    int depth;
    double cost;
    TreeNode* left;
    TreeNode* right;

    TreeNode(int d, double c) : depth(d), cost(c), left(nullptr), right(nullptr) {}
};

// Binary tree class
class BinaryTree {
public:
    BinaryTree(int depth);
    ~BinaryTree();
    void PrintTree();
    void SimulateWalkers(int numWalkers, int steps, const std::string& filename);
    TreeNode* GetRoot() const { return root; }
    
private:
    TreeNode* root;
    int maxDepth;
    void AssignCosts(TreeNode* node, int currentDepth);
    void PerturbPath(TreeNode* node);
    void DeleteTree(TreeNode* node);
    void PrintTree(TreeNode* node, int level);
    std::vector<TreeNode*> perturbedPathNodes;
};
