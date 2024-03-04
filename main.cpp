#include "graph.hpp"
#include "graphDynamics.hpp"

#include <experimental/iterator>
#include <fstream>

// To compile: g++ -I ./Eigen -g -o simu main.cpp graph.cpp graphDynamics.cpp

// Porre sulle ascisse il tempo e non il carico, così da
// poter graficare gli andamenti delle singole iterazioni.

using size = long unsigned int;

template <typename Container> int lastIndex(const Container &container) {
  if (container.empty()) {
    // Return -1 or any value that indicates an empty container
    return -1;
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

template <class T> void printVector(const std::vector<T> &Vec) { //??
  for (const auto &element : Vec) {
    std::cout << element << " ";
  }
  std::cout << std::endl;
}

// To print in a specific file the elements of the population matrix
void print_out(std::vector<std::vector<int>> const &population,
               std::string fname) {
  std::ofstream output_file(fname);
  if (output_file.is_open()) {
    for (size i = 0; i < population.size(); ++i) {
      int row_size = population[0].size();
      for (int j = 0; j < row_size;
           ++j) { // the rows should have all the same size
        output_file << " " << population[i][j] << " ";
      }
      output_file << "\n"; // put the endline outside to get the values of a
                           // single iteration in a row
    }
  }
}

// To use overload of << operator for vectors
template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &Vec) {
  // Printing all the elements using <<, or could use copy function without the
  // loop
  if (!Vec.empty()) {
    for (const auto &element : Vec) {
      os << element << " ";
    }
  }
  return os;
}

template <typename T>
void printVectorOut(const std::vector<T> &vector, std::string fname,
                    std::string gname) { // std::to_string(something you wnat
                                         // ot write at run time)
  std::ofstream output_file(fname);
  if (output_file.is_open()) {
    output_file << gname << "\n";
    for (size i = 0; i < vector.size(); ++i) {
      output_file << i << ", " << vector[i] << "\n";
    }
  }
}

template <typename T>
void printFluxOut(const std::vector<T> &vector, std::string fname,
                  std::string gname) { // std::to_string(something you wnat
                                       // ot write at run time)
  std::ofstream output_file(fname);
  if (output_file.is_open()) {
    output_file << gname << "\n";
    for (size i = 0; i < vector.size(); ++i) {
      if (i <= vector.size() / 2) {
        output_file << i << ", " << vector[i] << "\n";
      } else {
        int j = vector.size() - i - 1;
        output_file << j << ", " << vector[i] << "\n";
      }
    }
  }
}

template <typename T>
void slowPrintFluxOut(const std::vector<T> &vector, int start_unloading, std::string fname,
                  std::string gname) { // std::to_string(something you wnat
                                       // ot write at run time)
  std::ofstream output_file(fname);
  int count = 1;
  if (output_file.is_open()) {
    output_file << gname << "\n";
    for (size i = 0; i < vector.size(); ++i) {
      if (i <= start_unloading) {
        output_file << i << ", " << vector[i] << "\n";
      } else {
        // stampa un punto di scarico ogni 4, cosicché andat e ritorno abbiano gli stessi punti, 
        // corrispondenti allo stesso livello di carico.
        // è i % 4 perché start_unloading è 40 e noi vogliamo prendere le misure di scarico ogni volta che 
        // la popolazione del grafo corrisponde con quella delle misure prese durante il carico.
        // i % 4 Sarebbe da ricalcolare per valori diversi di start_unloading.
        //
        if (i % 4 == 0) {
          int j = start_unloading - count;
          ++count;
          output_file << j << ", " << vector[i] << "\n";
        }
      }
    }
  }
}

int main() {
  // Define the parameters for the network
  int nNode = 500; // Number of vertices in the graph
  int MinLink = 2;
  int MaxLink = 5;
  int nPart = nNode * 0;
  int nPart_fin = nNode * 9;
  int load = 250;
  int fill_load = 100;
  int empty_load = 25;
  int capacity = 10;
  int transport_Treshold = 5; 
  // caso normale
  int nIter_Load = (2 * nNode / load) * (nPart_fin - nPart) / nNode + 1; 
  // caso in cui la rete si scarica più lentamente
  // int nIter_Load = (1 + fill_load / empty_load) * (nNode / fill_load) * (nPart_fin - nPart) / nNode + 1;
  int start_unloading = (nNode / fill_load) * (nPart_fin - nPart) / nNode;
  std::cout << "nIter_load: " << nIter_Load << "\n";
  std::cout << "start_unloading: " << start_unloading << "\n";

  // Graph construction
  Graph graph(nNode);
  graph.randGraph(MinLink, MaxLink);
  graph.reconnectComponents(MaxLink);
  std::cout << "Created the graph.\n";
  std::vector<std::vector<double>> transitionMatrix =
      graph.createTransitionMatrix();
  std::cout << "Created the transition matrix.\n";

  // Synchronous initialization
  int nIter_synchronous = 10; // Time of evolution
  std::vector<std::vector<int>> population_synchronous;
  population_synchronous.reserve(nIter_synchronous + 1);
  std::vector<int> ini_state(nNode, 0);
  population_synchronous.push_back(ini_state);

  std::vector<int> synchronous_fluxHistory(nIter_synchronous, 0);
  // synchronous_fluxHistory.reserve(nIter_synchronous + 1);
  std::vector<long double> synchronous_RhoHistory(capacity + 1, 0);
  GraphDynamics synchronous_dynamics(capacity, population_synchronous,
                                     synchronous_fluxHistory,
                                     synchronous_RhoHistory);


  std::vector<int> synchronous_totalFluxHistory; /*(nIter_Load, 0);*/
  synchronous_totalFluxHistory.reserve(nIter_Load + 1);
  std::cout << "Finished synchronous initilization\n";
  // Asynchronous initialization
  int nIter_asynchronous = nIter_synchronous * nNode;
  std::vector<std::vector<int>> population_asynchronous;
  population_asynchronous.reserve(nIter_asynchronous + 1);
  population_asynchronous.push_back(ini_state);

  std::vector<int> asynchronous_fluxHistory(nIter_asynchronous + 1, 0);
  asynchronous_fluxHistory.reserve(nIter_asynchronous + 1);
  std::vector<long double> asynchronous_RhoHistory(capacity + 1, 0);
  GraphDynamics asynchronous_dynamics(capacity, population_asynchronous,
                                      asynchronous_fluxHistory,
                                      asynchronous_RhoHistory);

  // std::vector<std::vector<int>> asynchronous_totalFluxHistory;
  std::vector<int> asynchronous_totalFluxHistory;
  asynchronous_totalFluxHistory.reserve(nIter_Load + 1);
  std::cout << "Finished asynchronous initilization\n\n";

  for (int i = 0; i < nIter_Load; ++i) {
    /////////////////////////////////////////////////////////
    // Synchronous dynamics
    /////////////////////////////////////////////////////////
    std::cout << "Load phase: " << i << "\n";
    synchronous_dynamics.clearPopulation();
    // Normal case
    synchronous_dynamics.fillPopulation(i, nNode, nPart, nIter_Load, load);
    // Slow unloading
    // synchronous_dynamics.fillPopulation(i, nNode, nPart, start_unloading, fill_load, empty_load);
    // Redefine the number of iteration when there's slow unloading
    // if (i > start_unloading) {
    //   nIter_synchronous = 25;
    // }
    // std::cout << "iterazioni sincrone: " << Iter_synchronous << "\n";
    for (int ind = 0; ind < nIter_synchronous; ++ind) {

      std::random_device rd;
      std::mt19937 gen_synchronous(rd());
      std::uniform_int_distribution<int> int_distr(0, nNode - 1);
      std::uniform_real_distribution<double> real_distr(0.0, 1.0);
      synchronous_dynamics.createNewState();

      // Loop over nodes to move (at most) one particle per node
      for (int myNode = 0; myNode < nNode; ++myNode) {
        int pop = synchronous_dynamics.getNodePopulation(myNode);
        if (pop > 0) {
          std::vector<int> neigh_labels = graph.getNeighbours(myNode);
          std::vector<double> trans_probs = transitionMatrix[myNode];

          if (pop > transport_Treshold) {
            pop = transport_Treshold;
            for (int p = 0; p < pop; ++p) {
              double rand_real = real_distr(gen_synchronous);
              synchronous_dynamics.updatePopulation(trans_probs, rand_real,
                                                    myNode, neigh_labels);
            }
          } else {
            for (int p = 0; p < pop; ++p) {
              double rand_real = real_distr(gen_synchronous);
              synchronous_dynamics.updatePopulation(trans_probs, rand_real,
                                                    myNode, neigh_labels);
            }
          }
        }
      }
    }
    int average_flux = synchronous_dynamics.getAverageFlux(nIter_synchronous);
    synchronous_totalFluxHistory.push_back(average_flux);
    synchronous_dynamics.clearFluxHistory(nIter_synchronous);
    std::cout << "Synchornous Dynamics finished. \n";

    /////////////////////////////////////////////////////////
    // Asynchronous dynamics
    /////////////////////////////////////////////////////////
    asynchronous_dynamics.clearPopulation();
    // Normal case
    asynchronous_dynamics.fillPopulation(i, nNode, nPart, nIter_Load, load);
    // Slow unloading
    // asynchronous_dynamics.fillPopulation(i, nNode, nPart, start_unloading, fill_load, empty_load);
    // nIter_asynchronous = nIter_synchronous * nNode;

    for (int ind = 0; ind < nIter_asynchronous; ++ind) {

      asynchronous_dynamics.createNewState();
      std::random_device rd;
      std::mt19937 gen_asynchronous(rd());
      std::uniform_int_distribution<int> int_distr(0, nNode - 1);
      std::uniform_real_distribution<double> real_distr(0.0, 1.0);
      std::vector<int> non_empty_nodes;
      non_empty_nodes.reserve(nNode);
      for (int node = 0; node < nNode; ++node) {
        int pop = asynchronous_dynamics.getNodePopulation(node);
        if (pop != 0) {
          non_empty_nodes.push_back(node);
        }
      }
      // Apply Markov dynamics if there are particles
      if (!non_empty_nodes.empty()) {
        std::uniform_int_distribution<int> non_empty_distr(
            0, non_empty_nodes.size() - 1);
        int myNode = non_empty_nodes[non_empty_distr(gen_asynchronous)];

        std::vector<int> neigh_labels =
            graph.getNeighbours(myNode); // extract vector of neighboring rows
                                         // of myNode
        double rand_real = real_distr(gen_asynchronous);
        std::vector<double> trans_probs = transitionMatrix[myNode];

        // move particle from myNode
        int myNodePop = asynchronous_dynamics.getNodePopulation(myNode);
        if (myNodePop > transport_Treshold) {
          myNodePop = transport_Treshold;
          for (int p = 0; p < myNodePop; ++p) {
            double rand_real = real_distr(gen_asynchronous);
            asynchronous_dynamics.updatePopulation(trans_probs, rand_real, myNode,
                                                   neigh_labels);
          }
        } else {
          for (int p = 0; p < myNodePop; ++p) {
            double rand_real = real_distr(gen_asynchronous);
            asynchronous_dynamics.updatePopulation(trans_probs, rand_real, myNode,
                                                   neigh_labels);
          }
        }
      }

    }
    int average_flux2 = asynchronous_dynamics.getAverageFlux(nIter_synchronous);
    // std::this_thread::sleep_for(std::chrono::milliseconds(100));
    asynchronous_totalFluxHistory.push_back(average_flux2);
    asynchronous_dynamics.clearFluxHistory(nIter_asynchronous);
    std::cout << "Asynchornous Dynamics finished. \n\n";
  }
  // Synchronous measurement
  std::string strAdd = std::to_string(capacity) + "_" + std::to_string(nIter_Load);
  std::string str1 = "SyncFlux_.csv";
  str1.insert(9, strAdd);
  printFluxOut(synchronous_totalFluxHistory, str1, "Load, Flux");
  // slowPrintFluxOut(synchronous_totalFluxHistory, start_unloading , str1, "Load, Flux");
  std::cout << "synchronous flux: \n";
  printVector(synchronous_totalFluxHistory);

  // Asynchronous measurement
  std::string str2 = "AsyncFlux_.csv";
  str2.insert(10, strAdd);
  printFluxOut(asynchronous_totalFluxHistory, str2, "Load, Flux");
  // slowPrintFluxOut(asynchronous_totalFluxHistory, start_unloading, str2, "Load, Flux");
  std::cout << "asynchronous flux: \n";
  printVector(asynchronous_totalFluxHistory);

  return 0;
}

// To print out the population hsitory
// std::vector<std::vector<int>> popRecord =
// synchronous_dynamics.getPopulation(); int count = 0; for (size s = 0; s <
// popRecord.size(); ++s) {
//   if (s % 2 == 0) {
//     std::cout << s << ": \n";
//     std::cout << "f: ";
//     for (size j = 0; j < nNode; ++j) {
//     std::cout << popRecord[s].at(j) << " ";
//     }
//     std::cout << "\n";
//     ++count;
//   } else {
//     std::cout << "d: ";
//     for (size j = 0; j < nNode; ++j) {
//     std::cout << popRecord[s].at(j) << " ";
//     }
//     std::cout << "\n";
//     ++count;
//   }
// }
// std::cout << "Number of recorded states of the network: " << count << "\n";
