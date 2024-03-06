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
void printFluxOut(const std::vector<T> &vector1, const std::vector<T> &vector2,
                  const std::vector<T> &vector3, std::string fname,
                  std::string gname) { // std::to_string(something you wnat
                                       // ot write at run time)
  std::ofstream output_file(fname);
  if (output_file.is_open()) {
    output_file << gname << "\n";
    for (size i = 0; i < vector1.size(); ++i) {
      if (i <= vector1.size() / 2) {
        output_file << i << ", " << vector1[i] << ", " << vector2[i] << ", " << vector3[i] << "\n";
      } else {
        int j = vector1.size() - i - 1;
        output_file << j << ", " << vector1[i] << "\n";
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

template <typename T>
void printDevHistOut(const std::vector<T> &stdDeviations, std::string fname,
                  std::string gname) { // std::to_string(something you wnat
                                       // ot write at run time)
  std::vector<int> devHist(20, 0);
  for (size i = 0; i < stdDeviations.size(); ++i) {
    if (stdDeviations[i] > 0. && stdDeviations[i] <= 0.5) {
      devHist[0]++;
    } else if (stdDeviations[i] > 0.5 && stdDeviations[i] <= 1.) {
      devHist[1]++;
    } else if (stdDeviations[i] > 1. && stdDeviations[i] <= 1.5) {
      devHist[2]++;
    } else if (stdDeviations[i] > 1.5 && stdDeviations[i] <= 2) {
      devHist[3]++;
    } else if (stdDeviations[i] > 2 && stdDeviations[i] <= 2.5) {
      devHist[4]++;
    } else if (stdDeviations[i] > 2.5 && stdDeviations[i] <= 3.) {
      devHist[5]++;
    } else if (stdDeviations[i] > 3. && stdDeviations[i] <= 3.5) {
      devHist[6]++;
    } else if (stdDeviations[i] > 3.5 && stdDeviations[i] <= 4.) {
      devHist[7]++;
    } else if (stdDeviations[i] > 4. && stdDeviations[i] <= 4.5) {
      devHist[8]++;
    } else if (stdDeviations[i] > 4.5 && stdDeviations[i] <= 5.) {
      devHist[9]++;
    } else if (stdDeviations[i] > 5. && stdDeviations[i] <= 5.5) {
      devHist[10]++;
    } else if (stdDeviations[i] > 5.5 && stdDeviations[i] <= 6.){
      devHist[11]++;
    } else if (stdDeviations[i] > 6. && stdDeviations[i] <= 6.5) {
      devHist[12]++;
    } else if (stdDeviations[i] > 6.5 && stdDeviations[i] <= 7.) {
      devHist[13]++;
    } else if (stdDeviations[i] > 7. && stdDeviations[i] <= 7.5) {
      devHist[14]++;
    } else if (stdDeviations[i] > 7.5 && stdDeviations[i] <= 8.) {
      devHist[15]++;
    } else if (stdDeviations[i] > 8. && stdDeviations[i] <= 8.5) {
      devHist[16]++;
    } else if (stdDeviations[i] > 8.5 && stdDeviations[i] <= 9.) {
      devHist[17]++;
    } else if (stdDeviations[i] > 9. && stdDeviations[i] <= 9.5) {
      devHist[18]++;
    } else if (stdDeviations[i] > 9.5 && stdDeviations[i] <= 10.) {
      devHist[19]++;
    } else {
      std::cout << "above\n";
    }
  }
  std::ofstream output_file(fname);
  if (output_file.is_open()) {
    output_file << gname << "\n";
    for (size i = 0; i < devHist.size(); ++i) {
      output_file << i << ", " << devHist[i] << "\n";
    }
  }
}

template <typename T>
void printHistOut(const std::vector<T> &stdDeviations, std::string fname,
                  std::string gname) { // std::to_string(something you wnat
                                       // ot write at run time)
  std::ofstream output_file(fname);
  if (output_file.is_open()) {
    // output_file << gname << "\n";
    for (size i = 0; i < stdDeviations.size(); ++i) {
      output_file << stdDeviations[i] << "\n";
    }
  }
}


int main() {
  // Define the parameters for the network
  int nNode = 500; // Number of vertices in the graph
  int MinLink = 2;
  int MaxLink = 5;
  int nPart = nNode * 0;
  int nPart_fin = nNode * 20;
  int load = 250;
  // int fill_load = 100;
  // int empty_load = 25;
  int capacity = 20;
  int transport_Treshold = 5; 
  // caso normale
  int nIter_Load = (2 * nNode / load) * (nPart_fin - nPart) / nNode + 1; 
  // caso in cui la rete si scarica più lentamente
  // int nIter_Load = (1 + fill_load / empty_load) * (nNode / fill_load) * (nPart_fin - nPart) / nNode + 1;
  // int start_unloading = (nNode / fill_load) * (nPart_fin - nPart) / nNode;
  std::cout << "nIter_load: " << nIter_Load << "\n";
  // std::cout << "start_unloading: " << start_unloading << "\n";

  // Graph construction
  Graph graph(nNode);
  graph.randGraph(MinLink, MaxLink);
  graph.reconnectComponents(MaxLink);
  std::cout << "Created the graph.\n";
  std::vector<std::vector<double>> transitionMatrix =
      graph.createTransitionMatrix();
  std::cout << "Created the transition matrix.\n";

  // Synchronous initialization
  int nIter_synchronous = 100; // Time of evolution
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

  std::vector<int> synchronous_congested(nIter_Load + 1, 0);

  std::vector<double> synchronous_deviations;
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

  std::vector<int> asynchronous_congested(nIter_Load + 1, 0);

  // std::vector<std::vector<int>> asynchronous_totalFluxHistory;
  std::vector<double> asynchronous_deviations;
  std::vector<int> asynchronous_totalFluxHistory;
  asynchronous_totalFluxHistory.reserve(nIter_Load + 1);
  std::cout << "Finished asynchronous initilization\n\n";

  for (int i = 0; i < nIter_Load / 2; ++i) {
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
    std::vector<double> sync_stdDeviations(nNode, 0.);

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
    // int congestedCount = synchronous_dynamics.countCongested();
    // synchronous_congested.push_back(congestedCount);
    int average_flux = synchronous_dynamics.getAverageFlux(nIter_synchronous);
    // double average_syncDeviation = synchronous_dynamics.getStdDeviation(nIter_synchronous);
    double average_syncDeviation = synchronous_dynamics.getNodeDeviation(sync_stdDeviations);
    std::cout << "sync std deviation: " << average_syncDeviation << "\n";
    std::string string = "SyncHist_.csv";
    string.insert(9, std::to_string(i));
    // printDevHistOut(sync_stdDeviations, string, "Deviation, Occurence");
    printHistOut(sync_stdDeviations, string, "Deviation");

    synchronous_deviations.push_back(average_syncDeviation);
    synchronous_totalFluxHistory.push_back(average_flux);
    synchronous_dynamics.clearFluxHistory(nIter_synchronous);
    std::cout << "Synchornous Dynamics finished. \n";

    /////////////////////////////////////////////////////////
    // Asynchronous dynamics
    /////////////////////////////////////////////////////////
    // asynchronous_dynamics.clearPopulation();
    // // Normal case
    // asynchronous_dynamics.fillPopulation(i, nNode, nPart, nIter_Load, load);
    // std::cout << "filled\n";
    // // Slow unloading
    // // asynchronous_dynamics.fillPopulation(i, nNode, nPart, start_unloading, fill_load, empty_load);
    // // nIter_asynchronous = nIter_synchronous * nNode;

    // for (int ind = 0; ind < nIter_asynchronous; ++ind) {

    //   asynchronous_dynamics.createNewState();
    //   std::random_device rd;
    //   std::mt19937 gen_asynchronous(rd());
    //   std::uniform_int_distribution<int> int_distr(0, nNode - 1);
    //   std::uniform_real_distribution<double> real_distr(0.0, 1.0);
    //   std::vector<int> non_empty_nodes;
    //   non_empty_nodes.reserve(nNode);
    //   for (int node = 0; node < nNode; ++node) {
    //     int pop = asynchronous_dynamics.getNodePopulation(node);
    //     if (pop != 0) {
    //       non_empty_nodes.push_back(node);
    //     }
    //   }
    //   // Apply Markov dynamics if there are particles
    //   if (!non_empty_nodes.empty()) {
    //     std::uniform_int_distribution<int> non_empty_distr(
    //         0, non_empty_nodes.size() - 1);
    //     int myNode = non_empty_nodes[non_empty_distr(gen_asynchronous)];

    //     std::vector<int> neigh_labels =
    //         graph.getNeighbours(myNode); // extract vector of neighboring rows
    //                                      // of myNode
    //     double rand_real = real_distr(gen_asynchronous);
    //     std::vector<double> trans_probs = transitionMatrix[myNode];

    //     // move particle from myNode
    //     int myNodePop = asynchronous_dynamics.getNodePopulation(myNode);
    //     if (myNodePop > transport_Treshold) {
    //       myNodePop = transport_Treshold;
    //       for (int p = 0; p < myNodePop; ++p) {
    //         double rand_real = real_distr(gen_asynchronous);
    //         asynchronous_dynamics.updatePopulation(trans_probs, rand_real, myNode,
    //                                                neigh_labels);
    //       }
    //     } else {
    //       for (int p = 0; p < myNodePop; ++p) {
    //         double rand_real = real_distr(gen_asynchronous);
    //         asynchronous_dynamics.updatePopulation(trans_probs, rand_real, myNode,
    //                                                neigh_labels);
    //       }
    //     }
    //   }
    // }
    // // congestedCount = asynchronous_dynamics.countCongested();
    // // asynchronous_congested.push_back(congestedCount);
    // int average_flux2 = asynchronous_dynamics.getAverageFlux(nIter_synchronous);
    // double average_asyncDeviation = asynchronous_dynamics.getStdDeviation(nIter_synchronous, nNode);
    // std::cout << "async std deviation: " << average_asyncDeviation << "\n";
    // asynchronous_deviations.push_back(average_asyncDeviation);
    // // std::this_thread::sleep_for(std::chrono::milliseconds(100));
    // asynchronous_totalFluxHistory.push_back(average_flux2);
    // asynchronous_dynamics.clearFluxHistory(nIter_asynchronous);
    // std::cout << "Asynchornous Dynamics finished. \n\n";
  }
  // Synchronous measurement
  std::string strAdd = std::to_string(capacity) + "_" + std::to_string(nIter_Load);
  std::string strFlux1 = "SyncFlux_.csv";
  strFlux1.insert(9, strAdd);
  printFluxOut(synchronous_totalFluxHistory, strFlux1, "Load, Flux");
  std::string devFlux1 = "SyncDev_.csv";
  devFlux1.insert(8, strAdd);
  printFluxOut(synchronous_deviations, devFlux1, "Load, Deviation");

  // std::string strCong1 = "SyncCong_.csv";
  // strCong1.insert(9, strAdd);
  // printFluxOut(synchronous_congested, strCong1, "Load, Congested");
  // slowPrintFluxOut(synchronous_totalFluxHistory, start_unloading , strFlux1, "Load, Flux");
  std::cout << "synchronous dev: \n";
  printVector(synchronous_deviations);


  // Asynchronous measurement
  // std::string strFlux2 = "AsyncFlux_.csv";
  // strFlux2.insert(10, strAdd);
  // printFluxOut(asynchronous_totalFluxHistory, strFlux2, "Load, Flux");
  // std::string devFlux2 = "AsyncDev_.csv";
  // devFlux2.insert(8, strAdd);
  // printFluxOut(asynchronous_deviations, devFlux2, "Load, Deviation");
  // // std::string strCong2 = "AsyncCong_.csv";
  // // strCong2.insert(9, strAdd);
  // // printFluxOut(synchronous_congested, strCong2, "Load, Congested");
  // // slowPrintFluxOut(asynchronous_totalFluxHistory, start_unloading, strFlux2, "Load, Flux");
  // std::cout << "asynchronous dev: \n";
  // printVector(asynchronous_deviations);

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
