#ifndef GA_HPP
#define GA_HPP

#include "irp.hpp"
#include "aco.hpp"
#include "parameters.hpp"
#include "individual.hpp"
#include <vector>
#include <string>

Individual run_genetic_algorithm(const IRP& irp, const GA_Params& ga_params, const ACO_Params& aco_params, bool verbose);
bool check_aco_consistency(const std::vector<int>& deliveries, const std::vector<Route>& routes);


std::pair<Individual, Individual> one_point_crossover_customer(const Individual& a, const Individual& b, const IRP& irp);
std::pair<Individual, Individual> two_point_crossover_customer(const Individual& a, const Individual& b, const IRP& irp);
void advance_portion_mutation(Individual& ind, const IRP& irp, double pMutation);
Individual tournamentSelect(const std::vector<Individual>& pop, int k);
void generate_routes_for_individual(Individual& ind, const IRP& irp, const ACO_Params& aco_params);

void printDeliveriesMatrix(const Individual& ind, const IRP& irp);
void exportAndPlotRoutes(
    const IRP& irp,
    const Individual& best,
    const ACO_Params& acoParams,
    const std::string& dataFilename = "routes_data.txt",
    const std::string& pyScript = "plot_routes.py"
);

#endif // GA_HPP