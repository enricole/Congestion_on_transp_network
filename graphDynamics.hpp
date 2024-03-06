#ifndef GRAPHDYNAMICS_HPP
#define GRAPHDYNAMICS_HPP

#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <stdexcept>
#include <thread>
#include <vector>
#include <cmath>

class GraphDynamics {
  int capacity;
  std::vector<std::vector<int>> population;
  std::vector<int> fluxHistory;
  std::vector<long double> rhoHistory;

public:
  GraphDynamics(int cap, std::vector<std::vector<int>> pop,
                std::vector<int> fluxHis, std::vector<long double> rhoHis)
      : capacity(cap), population(pop), fluxHistory(fluxHis),
        rhoHistory(rhoHis) {}

  std::vector<std::vector<int>> getPopulation() { return population; }
  // Crea un nuovo stato della rete, un nuovo elemento del vettore population
  void createNewState();
  // Recupera la popolazione attuale di un nodo
  int getNodePopulation(int node);
  int sumPopulation();
  // Elimina tutti gli elementi del vettore population tranne l'ultimo
  void clearPopulation();
  // Fornisce un vettore contenente i nomi di tutti i nodi non-vuoti
  std::vector<int> findEmptyableNodes();
  // Carica o scarica la rete in base a quale momento della simulazione ci
  // troviamo
  void fillPopulation(const int i, int nNode, int nPart, int nIter_Load, int load);
  // caso in cui scarica più lentamente
  void fillPopulation(const int i, int nNode, int nPart, int start_unloading, int fill_load, int empty_load);
  // caso in cui vogliamo misurare la distr di prob della densità
  void fillPopulation(int l, int nNode);
  bool nodeIsNegative();
  // Crea il vettore di prob da usare in updatePopulation copiando la matrie di
  // transizione
  std::vector<double>
  createProbVector(const std::vector<std::vector<double>> &transMatrix,
                   int myNode);
  // Effettua lo scambio probabilistico fra due nodi
  void updatePopulation(const std::vector<double> &probs, double rand,
                       int myNode, const std::vector<int> &neigh_labels);

  // FLUX MEASURES
  // Svuota il vettore di flusso
  void clearFluxHistory(int nIter);
  // Fa una media fra tutti i flussi misurati nella corrente iterazione
  int getAverageFlux(int nIter);
  int getAverageFlux(int nIter, int start, int end);
  // synchronous version
  double getStdDeviation(int nIter);
  // asynchronous version
  double getStdDeviation(int nIter, int nNode);
  std::vector<int> getFluxHistory() { return fluxHistory; }

  // DENSITY DISTRIBUTION MEAUSRES
  // Svuota il vettore contenente le misure di RHO(n)
  void clearRhoHistory();
  // std::vector<int> nodesOverCapacity();
  // Effettua una misura sullo stato corrente di Rho(n)
  void updateRhoHistory();
  std::vector<long double> getRhoHistory() { return rhoHistory; }
  // Misura lo scarto
  std::vector<long double> getMeanDistribution(int start, int end);
  // To get the deviation between a set of iteration and the next one
  std::vector<long double> getDeviations(int coefficient);
  // Standard deviation calculated from the single nodes
  double getNodeDeviation(std::vector<double>& stdDeviations);
  int countCongested();

};
#endif
