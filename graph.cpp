#include "graph.hpp"

using size = long unsigned int;

void Graph::printAdjMatrix() {
  for (int i = 0; i < nNode; ++i) {
    std::sort(adj[i].begin(), adj[i].end());
    std::cout << "Node " << i << ": ";
    for (auto &neigh : adj[i]) {
      std::cout << neigh << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}

std::vector<int> Graph::getNeighbours(int myNode) { return adj[myNode]; }

std::vector<int> Graph::findAvailableNeighs(int MaxLink,
                                            std::vector<int> component) {
  size max = MaxLink;
  std::vector<int> availableNeighs;
  for (size i = 0; i < component.size(); ++i) {
    int node_label = component[i];
    if (adj[node_label].size() != max) {
      availableNeighs.push_back(node_label);
    }
  }
  return availableNeighs;
}

void Graph::randGraph(int MinLink, int MaxLink) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> real_distr(0.0, 1.0);

  for (int i = 0; i < nNode; ++i) {
    std::vector<int> Current_neigh = adj[i];
    int max = MaxLink - Current_neigh.size();
    int min = MinLink - Current_neigh.size();

    size treshold = MaxLink;
    if (Current_neigh.size() > treshold) {
      throw std::runtime_error{"The current node has too many links."};
    }
    std::vector<int> future_neigh;
    int SumOfLinks;
    std::vector<int> AllOfThem(nNode);
    std::iota(std::begin(AllOfThem), std::end(AllOfThem), 0);
    AllOfThem.erase(AllOfThem.begin() + i);

    std::sort(AllOfThem.begin(), AllOfThem.end());
    std::sort(adj[i].begin(), adj[i].end());
    std::vector<int> difference;
    std::set_difference(AllOfThem.begin(), AllOfThem.end(), adj[i].begin(),
                        adj[i].end(), std::back_inserter(difference));

    std::vector<int> availableNeighs =
        this->findAvailableNeighs(MaxLink, difference);
    double Plink = 0.9 / availableNeighs.size();

    do {
      future_neigh.clear();

      for (size j = 0; j < availableNeighs.size(); ++j) {
        int possibleNeigh = availableNeighs[j];
        if (real_distr(gen) <= Plink) {
          future_neigh.push_back(possibleNeigh);
        }
      }
      SumOfLinks = future_neigh.size();
    } while (SumOfLinks > max || SumOfLinks < min);

    for (int neigh : future_neigh) {
      // Eccezione sul selflinking
      if (i == neigh) {
        throw std::runtime_error{"Selflinking is not permitted."};
      }
      // Eccezione su link doppioni.
      if (adj[i].end() != std::find(adj[i].begin(), adj[i].end(), neigh)) {
        throw std::runtime_error{
            "The selected node is already linked to the i-th one."};
      }
      if (adj[neigh].end() !=
          std::find(adj[neigh].begin(), adj[neigh].end(), i)) {
        throw std::runtime_error{
            "The i-th node is already linked to the selected one."};
      }

      this->addEdge(i, neigh);
    }
  }

  for (int i = 0; i < nNode; ++i) {
    std::sort(adj[i].begin(), adj[i].end());

    // Eccezione sul numero di links
    if (adj[i].size() < 2 || adj[i].size() > 5) {
      throw std::runtime_error{"There's an incorrect number of links."};
    }

    for (auto &neigh : adj[i]) {
      // Eccezione sulla simmetria dei links.
      if (adj[neigh].end() ==
          std::find(adj[neigh].begin(), adj[neigh].end(), i)) {
        throw std::runtime_error{"The links are not symmetric."};
      }
    }
  }
}

void Graph::DFSUtil(int v, std::vector<bool> &visited,
                    std::vector<int> &component) {
  visited[v] = true;
  component.push_back(v);

  for (int neighbor : adj[v]) {
    if (!visited[neighbor]) {
      DFSUtil(neighbor, visited, component);
    }
  }
}

std::vector<int> Graph::findLargestComponent() {
  std::vector<bool> visited(nNode, false);
  std::vector<std::vector<int>> components;

  for (int v = 0; v < nNode; ++v) {
    if (!visited[v]) {
      std::vector<int> component;
      DFSUtil(v, visited, component);
      components.push_back(component);
    }
  }

  auto maxComponent =
      max_element(components.begin(), components.end(),
                  [](const std::vector<int> &a, const std::vector<int> &b) {
                    return a.size() < b.size();
                  });

  return *maxComponent;
}

std::vector<std::vector<int>> Graph::findComponents() {
  std::vector<bool> visited(nNode, false);
  std::vector<std::vector<int>> components;

  for (int v = 0; v < nNode; ++v) {
    if (!visited[v]) {
      std::vector<int> component;
      DFSUtil(v, visited, component);
      components.push_back(component);
    }
  }
  for (size i = 0; i < components.size(); ++i) {
    sort(components[i].begin(), components[i].end(),
         [](int a, int b) { return a < b; });
  }
  return components;
}

void Graph::reconnectComponents(int MaxLink) {
  std::random_device rd;
  std::mt19937 gen(rd());
  // Find ALL the connected components
  std::vector<std::vector<int>> components = this->findComponents();

  // Sort the components
  sort(components.begin(), components.end(),
       [](const std::vector<int> &a, const std::vector<int> &b) {
         return a.size() > b.size();
       });
  // Connect the smaller components to the biggest one
  std::uniform_real_distribution<double> real_distr(0.0, 1.0);
  // Repeat for all smaller components
  for (size smallComp = 1; smallComp < components.size(); ++smallComp) {
    // Store the set of new links for the smaller component
    std::vector<std::vector<int>> newLinks;
    // Store the number of new links
    int countNewLinks;
    do {
      newLinks.clear();
      for (int smallComp_Node :
           components[smallComp]) { // {smallComp[0], smallComp[1]}
        // Find nodes in the biggest components available for new links
        std::vector<int> availableNeighs =
            this->findAvailableNeighs(MaxLink, components[0]);
        // We'll loop over the biggest component's available nodes to build
        // new links, so we use their number. Set 0.5 to have on average half
        // of the smaller components' nodes to be connected to the biggest
        // component; could be set to be higher, especially if the components
        // are very small (usually they have just 3 nodes, we may want to have
        // at least 2 new links)
        double Plink = 0.1 / (availableNeighs.size());
        std::vector<int> Current_neigh = adj[smallComp_Node];
        int max = MaxLink - Current_neigh.size();
        std::vector<int> lbl_neigh;
        int SumOfLinks;
        do {
          lbl_neigh.clear();
          for (int bigComp_Node : availableNeighs) {
            // Randomly generate links between the node in the smaller comp
            // and the node in the biggest one
            // Until a good set of links is generated
            if (real_distr(gen) <= Plink) {
              lbl_neigh.push_back(bigComp_Node);
            }
          }
          SumOfLinks = lbl_neigh.size();
        } while (SumOfLinks > max);
        // Add the set of new links
        newLinks.push_back(
            lbl_neigh); // maybe can't do it, need to use vector of pointers
                        // and set null pointer if there are no neighs, ore
                        // use array of vector and put manually
      }
      // Count the number of new links for the smaller component
      countNewLinks = 0;
      for (size i = 0; i < newLinks.size(); ++i) {
        countNewLinks += newLinks[i].size();
      }
      // Until the smaller component has at least 1 link with the biggest one
    } while (countNewLinks < 1);

    // Add newLinks to graph
    for (size i = 0; i < newLinks.size(); ++i) {
      for (int neigh : newLinks[i]) {
        int smallComp_Node = components[smallComp][i];
        this->addEdge(smallComp_Node, neigh);
      }
    }
  }

  components = this->findComponents();
  if (components.size() > 1) {
    throw std::runtime_error{"The reconnection of the components didn't work."};
  }
}

void printVector(const Eigen::VectorXd &Vec) {
  for (const auto &element : Vec) {
    std::cout << element << " ";
  }
  std::cout << std::endl;
}

void printMatrix(const Eigen::MatrixXd &matrix) {
  double sum = 0.;
  for (int i = 0; i < matrix.rows(); ++i) {
    std::cout << "Node " << i << ": ";
    for (int j = 0; j < matrix.cols(); ++j) {
      std::cout << matrix(i, j) << " ";
      sum += matrix(i, j);
    }
    std::cout << "\n";
  }
  std::cout << "\n";
  // std::cout << sum << "\n";
}

void printMatrix(std::vector<std::vector<double>> matrix) {

  for (int i = 0; i < matrix.size(); ++i) {
    std::cout << "Node " << i << ": ";
    for (int j = 0; j < matrix[i].size(); ++j) {
      std::cout << matrix[i][j] << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
  // std::cout << sum << "\n";
}

std::vector<std::vector<double>> Graph::createTransitionMatrix() {
  for (int i = 0; i < nNode; ++i) {std::sort(adj[i].begin(), adj[i].end());}
  Eigen::MatrixXd transMatrix(nNode, nNode);
  for (int j = 0; j < nNode; ++j) {
    int degree = adj[j].size();
    size count = 0;
    for (int i = 0; i < nNode; ++i) {     
      if (i == adj[j][count]) {
        if (count < adj[j].size() - 1) {++count;}
        transMatrix(i, j) = double(1) / degree;
      } else {
        transMatrix(i, j) = 0.;
      }
    }
  }

  // Eccezioni
  Eigen::VectorXd colSums = transMatrix.colwise().sum();
  std::vector<size> vCount(nNode, 0);
  
  for (int j = 0; j < nNode; ++j) {
    if (std::abs(colSums[j] - 1.) > 0.000001) {
      std::cerr << "Node: " << j << "\n";
      printMatrix(transMatrix);
      throw std::runtime_error{"Vanilla probabilities are not normalized to 1 for this node."};
    }
    for (int i = 0; i < nNode; ++i) {
      if (transMatrix(i, j) != 0) {++vCount[j];}
    }
  }
  for (int i = 0; i < nNode; ++i) {
    if (vCount[i] != adj[i].size()) {
      throw std::runtime_error{"The transition matrix row has more elements than the adjiacency matrix row."};
    }
  }

  // Compute eigenvalues and eigenvectors
  Eigen::EigenSolver<Eigen::MatrixXd> solver(transMatrix);
  // Find the index of the maximum eigenvalue
  int maxEigenvalueIndex; //Eccezione che verifichi che l'autovalore sia 1
  solver.eigenvalues().real().maxCoeff(&maxEigenvalueIndex);
  // Extract the eigenvector corresponding to the maximum eigenvalue
  Eigen::VectorXd maxEigenvector =
      solver.eigenvectors().real().col(maxEigenvalueIndex);


  // Normalizza sul massimo elemento dell'autovettore
  double maxCoeff = maxEigenvector.maxCoeff();
  for (int i = 0; i < maxEigenvector.size(); ++i) {
    maxEigenvector[i] /= maxCoeff;
  }
  // Moltiplica la colonna j per il j-esimo elemento dell'autovettore
  for (int j = 0; j < maxEigenvector.size(); ++j) {
    transMatrix.col(j) *= maxEigenvector[j];
  }
  
  // Eccezione
  // Somma sugli elementi di colonna per ottenere l'autovettore 
  Eigen::VectorXd rowSum = transMatrix.rowwise().sum();
  Eigen::VectorXd colSum = transMatrix.colwise().sum();
  for (int j = 0; j < nNode; ++j) {
    if (std::abs(rowSum[j] - maxEigenvector[j]) > 0.000001) {
      throw std::runtime_error{"The sum over the j-th row doesn't return the j-th element of the eigenvector."};
    }

    if (std::abs(colSum[j] - maxEigenvector[j]) > 0.000001) {
      throw std::runtime_error{"The sum over the j-th column doesn't return the j-th element of the eigenvector."};
    }
  }

  std::vector<std::vector<double>> LILtransMatrix(nNode);
  for (int j = 0; j < nNode; ++j) {
    LILtransMatrix[j].push_back(0.0);
    for (int i = 0; i < nNode; ++i) {
      if (transMatrix(i, j) != 0) {
        LILtransMatrix[j].push_back(transMatrix(i, j));
      }
    }
  } 
  return LILtransMatrix;
}
