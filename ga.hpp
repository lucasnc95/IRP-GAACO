#ifndef GA_HPP
#define GA_HPP

#include "irp.hpp"
#include "aco.hpp"
#include "parameters.hpp"
#include "individual.hpp"
#include <vector>
#include <string>

// Declaração da função principal do GA
Individual run_genetic_algorithm(const IRP& irp, const GA_Params& ga_params, const ACO_Params& aco_params, bool verbose);

// Função de criação de indivíduo (Genótipo)
Individual make_new_heuristic_individual(const IRP& irp);

// Operadores genéticos
std::pair<Individual, Individual> one_point_crossover_customer(const Individual& a, const Individual& b, const IRP& irp);
std::pair<Individual, Individual> two_point_crossover_customer(const Individual& a, const Individual& b, const IRP& irp);
void advance_portion_mutation(Individual& ind, const IRP& irp, double pMutation);
Individual tournamentSelect(const std::vector<Individual>& pop, int k);
// Funções utilitárias
void printDeliveriesMatrix(const Individual& ind, const IRP& irp);
void exportAndPlotRoutes(
    const IRP& irp,
    const Individual& best,
    const ACO_Params& acoParams,
    const std::string& dataFilename = "routes_data.txt",
    const std::string& pyScript = "plot_routes.py"
);

// --- NOVA FUNÇÃO (RENOMEADA) ---
/**
 * @brief (Função APENAS para construir rotas)
 * Preenche o 'routes_per_period' de um indivíduo usando o ACO,
 * com base no seu 'deliveries'.
 */
void build_routes_for_individual(Individual& ind, const IRP& irp, const ACO_Params& aco_params);
// --- FIM DA NOVA FUNÇÃO ---

#endif // GA_HPP