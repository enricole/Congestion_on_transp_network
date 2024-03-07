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

int main() {
  // Define the parameters for the network
  int nNode = 500; // Number of vertices in the graph
  int MinLink = 2;
  int MaxLink = 5;
  int capacity = 20;
  int transport_Treshold =
      1; 
  int l = 15; // load level, equivalent to l*nNode particles in the system
  // Graph construction
  Graph graph(nNode);
  graph.randGraph(MinLink, MaxLink);
  graph.reconnectComponents(MaxLink);
  std::cout << "Created the graph.\n";
  std::vector<std::vector<double>> transitionMatrix =
      graph.createTransitionMatrix();
  std::cout << "Created the transition matrix.\n";

  // Synchronous initialization
  int nIter_synchronous = 4000; // Time of evolution
  std::vector<std::vector<int>> population_synchronous;
  population_synchronous.reserve(nIter_synchronous + 1);
  std::vector<int> ini_state(nNode, 0);
  population_synchronous.push_back(ini_state);

  std::vector<int> synchronous_fluxHistory(nIter_synchronous, 0);
  synchronous_fluxHistory.reserve(nIter_synchronous + 1);
  std::vector<long double> synchronous_RhoHistory(capacity + 5, 0);
  GraphDynamics synchronous_dynamics(capacity, population_synchronous,
                                     synchronous_fluxHistory,
                                     synchronous_RhoHistory);




  // Asynchronous initialization
  int nIter_asynchronous = nIter_synchronous * nNode;
  std::vector<std::vector<int>> population_asynchronous;
  population_asynchronous.reserve(nIter_asynchronous + 1);
  population_asynchronous.push_back(ini_state);

  std::vector<int> asynchronous_fluxHistory(nIter_asynchronous, 0);
  asynchronous_fluxHistory.reserve(nIter_asynchronous + 1);
  std::vector<long double> asynchronous_RhoHistory(capacity + 1, 0);
  GraphDynamics asynchronous_dynamics(capacity, population_asynchronous,
                                      asynchronous_fluxHistory,
                                      asynchronous_RhoHistory);


    /////////////////////////////////////////////////////////
    // Synchronous dynamics
    /////////////////////////////////////////////////////////
    std::cout << "\nSynchornous Dynamics: \n\n";
    synchronous_dynamics.fillPopulation(l, nNode);

    for (int ind = 0; ind < nIter_synchronous; ++ind) {
      if (ind % (nIter_synchronous / 2) == 0) {
        std::cout << "Point in the program: " << ind << "\n";
      } 
      int min = nIter_synchronous * 5 / 20;
      if (ind > min && ind % 2 == 0) {
        synchronous_dynamics.updateRhoHistory();
      }
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
          double rand_real = real_distr(gen_synchronous);
          
          if (pop > transport_Treshold) {
            pop = transport_Treshold;
            for (int p = 0; p < pop; ++p) {
              synchronous_dynamics.updatePopulation(trans_probs, rand_real,
                                                    myNode, neigh_labels);
            }
          } else {
            for (int p = 0; p < pop; ++p) {
              synchronous_dynamics.updatePopulation(trans_probs, rand_real,
                                                    myNode, neigh_labels);
            }
          }
        }
      }
    }

    /////////////////////////////////////////////////////////
    // Asynchronous dynamics
    /////////////////////////////////////////////////////////
    std::cout << "\nAsynchornous Dynamics: \n\n";
    asynchronous_dynamics.fillPopulation(l, nNode);

    for (int ind = 0; ind < nIter_asynchronous; ++ind) {
      if (ind % (nIter_asynchronous / 2) == 0) {
        std::cout << "Point in the program: " << ind << "\n";
      }    
      int min = nIter_asynchronous * 5 / 20;
      if (ind > min && ind % 100 == 0) {
        asynchronous_dynamics.updateRhoHistory();
      }
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
            asynchronous_dynamics.updatePopulation(trans_probs, rand_real, myNode,
                                                   neigh_labels);
          }
        } else {
          for (int p = 0; p < myNodePop; ++p) {           
            asynchronous_dynamics.updatePopulation(trans_probs, rand_real, myNode,
                                                   neigh_labels);
          }
        }
      }

    }
    std::cout << "end\n\n";

  
  
  // Synchronous measurement
  std::vector<long double> syn_rhoHistory = synchronous_dynamics.getRhoHistory();
  Normalize(syn_rhoHistory);
  std::string string1 = "SyncRho_.csv";
  std::string lstring = std::to_string(l);
  string1.insert(8, lstring);
  printVectorOut(syn_rhoHistory, string1, "State, Probability");
  std::cout << " synchronous rho: \n";
  printVector(syn_rhoHistory);

  // std::cout << "Sync deviations: \n";
  // std::vector<long double> sync_deviations = synchronous_dynamics.getDeviations(1);
  // std::string devString1 = "SyncDeviations_.csv";
  // devString1.insert(15, lstring);
  // printVectorOut(sync_deviations, devString1, "iterations, deviation");
  // printVector(sync_deviations);



  // Asynchronous measurement
  std::vector<long double> asyn_rhoHistory = asynchronous_dynamics.getRhoHistory();
  Normalize(asyn_rhoHistory);
  std::string string2 = "AsyncRho_.csv";
  string2.insert(9, lstring);
  printVectorOut(asyn_rhoHistory, string2, "State, Probability");
  std::cout << "asynchrnous rho: \n";
  printVector(asyn_rhoHistory);

  // std::cout << "Async deviations: \n";
  // std::vector<long double> async_deviations = asynchronous_dynamics.getDeviations(nNode);
  // std::string devString2 = "AsyncDeviations_.csv";
  // devString2.insert(16, lstring);
  // printVectorOut(async_deviations, devString2, "iterations, deviation");
  // printVector(async_deviations);

  return 0;
}

