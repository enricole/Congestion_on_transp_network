#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <algorithm>
#include <iostream>
#include <random>
#include <stack>
#include <stdexcept>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

class Graph {
  int nNode;
  std::vector<std::vector<int>> adj;

public:
  Graph(int node) : nNode(node) { adj.resize(node); }

  void addEdge(int u, int v) {
    adj[u].push_back(v);
    adj[v].push_back(u);
  }

  void clearGraph() {
    for (int i = 0; i < nNode; ++i) {
      adj[i].clear();
    }
  }
  // Prints out the adjacency matrix
  void printAdjMatrix();
  // Returns yhe neighbours of a specific node
  std::vector<int> getNeighbours(int myNode);
  // Find nodes available for new links from specified connected component
  std::vector<int> findAvailableNeighs(int MaxLink, std::vector<int> component);
  // Creates a random graph in which nodes have between MinLink and MaxLInk links
  void randGraph(int MinLink, int MaxLink);
  // Algorithm implemented to map out a component of the graph
  void DFSUtil(int v, std::vector<bool> &visited, std::vector<int> &component);
  // Returns the largest component of the graph
  std::vector<int> findLargestComponent();
  // Return all the different components of the graph
  std::vector<std::vector<int>> findComponents();
  // Method that randomly connects the samller components with the largest one
  void reconnectComponents(int MaxLink);
  // Returns the transition matrix
  std::vector<std::vector<double>> createTransitionMatrix();
};

#endif