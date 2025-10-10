#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <string>

// Parâmetros do Algoritmo Genético
struct GA_Params {
    int popSize = 40;
    int nGen = 200;
    double pCrossover = 0.9;
    double pMutation = 0.15; // Probabilidade da mutação simples
    int tournamentK = 3;
    int stagnation_threshold = 20;
    int crossover_max_tries = 10; // Tentativas para gerar filho factível
};

// Parâmetros do ACO
struct ACO_Params {
    int nAnts = 20;
    int nIter = 100;
    double alpha = 1.0;
    double beta = 2.0;
    double rho = 0.1;
    double Q = 1000.0;
    double pLocalSearch = 1.0;
};

void load_parameters_from_file(const std::string& filename, GA_Params& ga_params, ACO_Params& aco_params);

#endif // PARAMETERS_HPP