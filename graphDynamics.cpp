#include "graphDynamics.hpp"
using size = long unsigned int;
template <typename Container> int lastIndex(const Container &container) {
  if (container.empty()) {
    throw std::runtime_error{"The container is empty."};
  }
  return container.size() - 1;
}

template <typename T> void Normalize(std::vector<T> &vector) {
  T sum = 0;
  for (size i = 0; i < vector.size(); ++i) {
    sum += vector[i];
  }

  for (auto &elem : vector) {
    elem /= sum;
  }
}

template <class T> void printVector(const std::vector<T> &Vec) {
  for (const auto &element : Vec) {
    std::cout << element << " ";
  }
  std::cout << std::endl;
}

template<typename T>
void printLIL(const std::vector<std::vector<T>>& vec) {
    for (int i = 0; i < vec.size(); ++i) {
        std::cout << i << ": ";
        for (const auto& element : vec[i]) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
}

void GraphDynamics::createNewState() {
  int index = lastIndex(population);
  population.push_back(population[index]);
}

int GraphDynamics::getNodePopulation(int node) {
  int index = lastIndex(population);
  return population[index - 1][node];
}

int GraphDynamics::sumPopulation() {
  int sum = 0;
  int index = lastIndex(population);
  for (auto &node : population[index]) {
    sum += node;
  }
  return sum;
}

void GraphDynamics::clearPopulation() {
  // Elimina tutte le iterazioni del vettore population tranne l'ultima.
  population.erase(population.begin(), population.end() - 1);
}

std::vector<int> GraphDynamics::findEmptyableNodes(int load) {
  int index = lastIndex(population);
  std::vector<int> emptyableNodes;
  for (int i = 0; i < population[index].size(); ++i) {
    if (population[index][i] >= load) {
      emptyableNodes.push_back(i);
    }
  }
  for (int i = 0; i < emptyableNodes.size(); ++i) {
    int node = emptyableNodes[i];
    if (population[index][node] < load) {
      throw std::runtime_error{"findEmptyableNodes didn't work."};
    }
  }
  return emptyableNodes;
}

bool GraphDynamics::nodeIsNegative() {
  int index = lastIndex(population);
  for (auto &node : population[index]) {
    if (node < 0) {
      return true;
    }
  }
  return false;
}

void GraphDynamics::fillPopulation(const int i, const int nNode, int nPart,
                                   int nIter_Load, int load) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> int_distr(0, nNode - 1);
  int node = 0;
  int index = lastIndex(population);
  // std::cout << "fillPopulation called on the index: " << index << "\n";

  if (i == 0) {
    for (int k = 0; k < nPart; ++k) {
      do {
        node = int_distr(gen);
      } while (population[index][node] >= capacity);
      population[index][node] += 1;
    }
    int sum = sumPopulation();
    if (sum != nPart) {
      throw std::runtime_error{
          "Population wasn't filled correctly, in the stage 0."};
    }
  } else if (i > 0 && i <= nIter_Load / 2) {
    for (int k = 0; k < load; ++k) {
      // if (i == 9) {std::cout << "\niterazione: " << k << "\n";}
      do {
        node = int_distr(gen);

      } while (population[index][node] >= capacity);

      population[index][node] += 1;      
    }
    int sum = sumPopulation();
    if (sum != nPart + i * load) {
      printVector(population[index - 1]);
      printVector(population[index]);
      std::cerr << "sum: " << sum << "\n";
      std::cerr << "nPart + i * load: " << nPart + i * load << "\n";
      throw std::runtime_error{"Network wasn't filled correctly."};
    }
  } else {
    for (int k = 0; k < load; ++k) {
      do {
        node = int_distr(gen);
      } while (population[index][node] == 0);
      population[index][node] -= 1;  
    }
    
    int sum = sumPopulation();
    int j = nIter_Load - i - 1;
    if (sum != nPart + j * load) {
      printVector(population[index - 1]);
      printVector(population[index]);
      std::cerr << "sum: " << sum << "\n";
      std::cerr << "nPart + j * load: " << nPart + j * load << "\n";
      throw std::runtime_error{"Network wasn't emptied correctly."};
    }
    if (nodeIsNegative()) {
      throw std::runtime_error{"a node is negative."};
    }
  }
}

// Slow unloading
void GraphDynamics::fillPopulation(const int i, const int nNode, int nPart,
                                   int start_unloading, int fill_load, int empty_load) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> int_distr(0, nNode - 1);
  int node = 0;
  int index = lastIndex(population);
  // std::cout << "fillPopulation called on the index: " << index << "\n";

  if (i == 0) {
    for (int k = 0; k < nPart; ++k) {
      do {
        node = int_distr(gen);
      } while (population[index][node] >= capacity);
      population[index][node] += 1;
    }
    int sum = sumPopulation();
    if (sum != nPart) {
      throw std::runtime_error{
          "Population wasn't filled correctly, in the stage 0."};
    }
  } else if (i > 0 && i <= start_unloading) {
    for (int k = 0; k < fill_load; ++k) {
      do {
        node = int_distr(gen);
      } while (population[index][node] >= capacity);
      population[index][node] += 1;      
    }
    int sum = sumPopulation();
    if (sum != nPart + i * fill_load) {
      printVector(population[index - 1]);
      printVector(population[index]);
      std::cerr << "sum: " << sum << "\n";
      std::cerr << "nPart + i * fill_load: " << nPart + i * fill_load << "\n";
      throw std::runtime_error{"Network wasn't filled correctly."};
    }
  } else { // 
    std::cout << "sta scaricando.\n";
    int sum1 = sumPopulation();
    for (int k = 0; k < empty_load; ++k) {
      do {
        node = int_distr(gen);
      } while (population[index][node] == 0);
      population[index][node] -= 1;  
    }
    int sum2 = sumPopulation();
    if (sum1 != sum2 + empty_load) {
      printVector(population[index - 1]);
      printVector(population[index]);
      std::cerr << "sum1: " << sum1 << "\n";
      std::cerr << "sum2 + empty_load: " << sum2 + empty_load << "\n";
      throw std::runtime_error{"Network wasn't emptied correctly."};
    }
    if (nodeIsNegative()) {
      throw std::runtime_error{"a node is negative."};
    }
  }
}

// Version of the method to measure prob distr of Rho
void GraphDynamics::fillPopulation(int l, const int nNode) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> int_distr(0, nNode - 1);
  int node = 0;
  int index = lastIndex(population);


  for (int k = 0; k < l * nNode; ++k) {
    do {
      node = int_distr(gen);
    } while (population[index][node] >= capacity);

    population[index][node] += 1;
  }
  int sum = sumPopulation();
  if (sum != l * nNode) {
    throw std::runtime_error{
        "Population wasn't filled correctly, in the stage 0."};
  }
}
// Fai un eccezione in cui conti gli elementi tolti e verifichi che la taglia
// finale concordi
std::vector<double> GraphDynamics::createProbVector(
    const std::vector<std::vector<double>> &transMatrix, int myNode) {
  // Create the prob vector
  std::vector<double> trans_probs = {0.0};
  for (size i = 0; i < transMatrix[myNode].size(); ++i) {
    double prob = transMatrix[myNode].at(i);
    trans_probs.push_back(prob);
  }
  return trans_probs;
}

// sum elements of a vector within specified range //You can use an algorithm
// for this
double sumPartial(const std::vector<double> &V, int start, int end) {
  double sum = 0.0;
  for (int i = start; i <= end; ++i)
    sum += V[i];
  return sum;
}

// move particle from myNode to a randomly chosen neighbour
// make two version based on the transport capacity
void GraphDynamics::updatePopulation(const std::vector<double> &probs,
                                    double rand, int myNode,
                                    const std::vector<int> &neigh_labels) {
  int index = lastIndex(population);
  // int fluxInd = lastIndex(fluxHistory);
  int ini_pop = sumPopulation();
  for (size i = 0; i < probs.size() - 1; ++i) {
    if (rand > sumPartial(probs, 0, i) && rand <= sumPartial(probs, 0, i + 1) &&
        population[index - 1][neigh_labels[i]] < capacity) {

        population[index][neigh_labels[i]]++;
        population[index][myNode]--;
        // Conteggio per tenere conto del numero di particelle mosse ad ogni
        // iterazione
        ++fluxHistory[index - 1];
      int fin_pop = sumPopulation();

      if (nodeIsNegative()) {
        throw std::runtime_error{"A node has negative population."};
      }

      if (ini_pop != fin_pop) {
        throw std::runtime_error{
            "In an exchange the total number of vehicles changed."};
      }

    }
  }

}

// Svuota il vettore storia e lo reinizializza
void GraphDynamics::clearFluxHistory(int nIter) {
  fluxHistory.clear();
  std::vector<int> initializer(nIter, 0);
  fluxHistory = initializer;
}

int GraphDynamics::getAverageFlux(int nIter) {
  int sum = 0;
  for (auto &hist : fluxHistory) {
    sum += hist;
  }
  return int(sum / nIter);
}
// Svuota il vettore storia e lo reinizializza
void GraphDynamics::clearRhoHistory() {
  rhoHistory.clear();
  std::vector<long int> initializer(capacity, 0);
  // rhoHistory.reserve(capacity);
  std::copy(initializer.begin(), initializer.end(),
            std::back_inserter(rhoHistory));
}

// Trova i nodi con una popolazione oltre la capacità.
// std::vector<int> GraphDynamics::nodesOverCapacity() {
//   int index = lastIndex(population);
//   int nNode = population[index].size();
//   std::vector<int> nodesOver;
//   for (int i = 0; i < nNode; ++i) {
//     if (population[index][i] > capacity) {
//       nodesOver.push_back(i);
//     }
//   }
//   return nodesOver;
// }

void GraphDynamics::updateRhoHistory() {
  int index = lastIndex(population);
  int nNode = population[index].size();
  // There are nodes that are over the capacity because of fillPopulation, and
  // we have to ignore them
  // std::vector<int> nodesOver = nodesOverCapacity();
  for (int i = 0; i < nNode; ++i) {
      int j = population[index][i];
      rhoHistory[j]++;
  }
}


std::vector<long double> GraphDynamics::getMeanDistribution(int start, int end) {
  std::vector<long double> meanRhoDistribution(capacity + 5, 0);
  int nNode = population[start].size();
  for (int k = start; k < end; ++k) {
    for (int i = 0; i < nNode; ++i) {
        int j = population[k][i];
        meanRhoDistribution[j]++;
    }
  }
  Normalize(meanRhoDistribution);
  return meanRhoDistribution;
}

int removeUnits(int number) {
    return (number / 10) * 10; // Remove units digit by integer division and then multiply by 10
}

std::vector<long double> GraphDynamics::getDeviations(int coefficient) {
  int n = population.size() / 3;
  std::cout << "n: " << n << "\n";
  int numb = removeUnits(n) / (10 * coefficient);
  std::cout << "numb: " << numb << "\n";
  std::vector<std::vector<long double>> distributionVector;
  for (int i = 0; i < numb + 1; ++i) {
    int start = i * 20 * coefficient;
    int end = (i * 20 + 10) * coefficient;
    
    std::vector<long double> meanDist = getMeanDistribution(start, end);
    distributionVector.push_back(meanDist);
  }

  std::vector<long double> deviations;
  int distSize = distributionVector[0].size();
  for (int i = 1; i < distributionVector.size(); ++i) {
    long double sum = 0;
    for (int j = 0; j < distSize; ++j) {
      long double difference = distributionVector[i][j] - distributionVector[i - 1][j];
      long double squared = pow(difference, 2);
      sum += squared;
    }
    deviations.push_back(sum);
  }
  
  return deviations;
}

// Prova a fare le medie di flusso sugli stati stazionari per vedere che il fenomeno di isteresi si neutralizza.

// Trasla le x dei grafici

// Nel capitolo 3 parli della misura del flusso in stati transistivi e dell'importanza
// di non raggiungere lo stato stazionario. A seguito di uno studio preliminare sul raggiungimento 
// dello stato stazionario abbiamo detereminato un intervallo adeguato d iterazioni. 

// Nel capitolo 4 nell'introduzione mettici i grafici. 
