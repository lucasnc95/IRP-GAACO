#ifndef GA_HPP
#define GA_HPP

#include "irp.hpp"
#include "aco.hpp"
#include "parameters.hpp"
#include "individual.hpp"
#include <vector>
#include <string>

Individual run_genetic_algorithm(const IRP& irp, const GA_Params& ga_params, const ACO_Params& aco_params, bool verbose);

// Funções do ciclo do GA
Individual make_stock_up_individual(const IRP& irp);
Individual make_just_in_time_individual(const IRP& irp);
Individual crossover(const Individual& a, const Individual& b, const IRP& irp);
void advance_portion_mutation(Individual& ind, const IRP& irp); // <-- MUDANÇA: Nova e única mutação
Individual tournamentSelect(const std::vector<Individual>& pop, int k);


// Funções de relatório/auxiliares
void printDeliveriesMatrix(const Individual& ind, const IRP& irp);
void exportAndPlotRoutes(
    const IRP& irp,
    const Individual& best,
    const ACO_Params& acoParams,
    const std::string& dataFilename = "routes_data.txt",
    const std::string& pyScript = "plot_routes.py"
);

#endif // GA_HPP