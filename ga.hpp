#ifndef GA_HPP
#define GA_HPP

#include "irp.hpp"
#include "aco.hpp"
#include "parameters.hpp"
#include "individual.hpp"
#include <vector>
#include <string>

Individual run_genetic_algorithm(const IRP& irp, const GA_Params& ga_params, const ACO_Params& aco_params, bool verbose);

// --- NOVAS FUNÇÕES DO CICLO DO GA ---
// Nova (e única) função de geração de população
Individual make_simple_random_individual(const IRP& irp);

// Novos operadores de crossover
Individual one_point_crossover_customer(const Individual& a, const Individual& b, const IRP& irp);
Individual two_point_crossover_customer(const Individual& a, const Individual& b, const IRP& irp);
void advance_portion_mutation(Individual& ind, const IRP& irp);
// Nova (e única) função de mutação
void simple_random_mutation(Individual& ind, const IRP& irp);

Individual tournamentSelect(const std::vector<Individual>& pop, int k);

// --- FUNÇÕES DE RELATÓRIO/AUXILIARES ---
void printDeliveriesMatrix(const Individual& ind, const IRP& irp);
void exportAndPlotRoutes(
    const IRP& irp,
    const Individual& best,
    const ACO_Params& acoParams,
    const std::string& dataFilename = "routes_data.txt",
    const std::string& pyScript = "plot_routes.py"
);

#endif // GA_HPP