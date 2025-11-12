/*
 * ARQUIVO MODIFICADO: ds_operator.hpp
 *
 * Adicionada a declaração de 'run_simple_reinsertion_search'.
 */

#ifndef DS_OPERATOR_HPP
#define DS_OPERATOR_HPP

#include <vector>
#include <gurobi_c++.h> 
#include "irp.hpp"
#include "individual.hpp"
#include "parameters.hpp"

// ... (Structs GurobiPWLCurve, GurobiInventoryResult, ReinsertionData permanecem iguais) ...
struct GurobiPWLCurve {
    std::vector<double> q_points;    
    std::vector<double> cost_points; 
};

struct GurobiInventoryResult {
    std::vector<double> q;  
    std::vector<double> S;  
    double totalCost;       
};

struct ReinsertionData {
    int customer_id = -1; 
    double cost_before_removal;
    std::vector<std::vector<Route>> routes_before_removal;
    double cost_after_removal; 
    std::vector<std::vector<Route>> routes_after_removal;
    Individual solution_without_customer; 
    std::vector<GurobiPWLCurve> insertion_cost_curves; 
    std::vector<long> max_q_inventory; 

    ReinsertionData(int nPeriods = 0) {
        cost_before_removal = 0.0;
        cost_after_removal = 0.0;
        routes_before_removal.resize(nPeriods);
        routes_after_removal.resize(nPeriods);
        insertion_cost_curves.resize(nPeriods);
        max_q_inventory.resize(nPeriods, 0);
    }
};

// --- Funções Existentes ---
GurobiInventoryResult computePerfectInventory_Gurobi(
    const Customer& cust,
    const Depot& depot,
    const std::vector<double>& reinsertion_cost_fixed,
    const std::vector<GurobiPWLCurve>& reinsertion_curves,
    int T, double qmax
);

ReinsertionData calculate_reinsertion_data(
    const Individual& original_ind, 
    const IRP& irp, 
    const ACO_Params& aco_params,
    int c_id_internal
);

void run_vnd_ds_operator(
    Individual& original_sol, 
    const IRP& irp, 
    const ACO_Params& aco_params,
    bool verbose
);

// --- NOVA FUNÇÃO ---
/**
 * @brief Busca Local Rápida (Simple Reinsertion).
 * Remove o cliente, otimiza o inventário com Gurobi baseando-se nas
 * rotas existentes e reinsere usando Cheapest Insertion (sem reroteirizar).
 */
void run_simple_reinsertion_search(
    Individual& sol, 
    const IRP& irp, 
    const ACO_Params& aco_params,
    bool verbose
);

void print_routes_for_period(const std::vector<Route>& routes, const IRP& irp, int t);

#endif // DS_OPERATOR_HPP