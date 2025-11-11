#ifndef DS_OPERATOR_HPP
#define DS_OPERATOR_HPP

#include <vector>
#include <gurobi_c++.h> // Incluído para o Gurobi
#include "route_builder.hpp"
#include "irp.hpp"
#include "individual.hpp"
#include "parameters.hpp"


struct GurobiPWLCurve {
    std::vector<double> q_points;    // pontos x (quantidade)
    std::vector<double> cost_points; // pontos y (custo de frete)
};


struct GurobiInventoryResult {
    std::vector<double> q;  // entregas ótimas por período
    std::vector<double> S;  // inventários ótimos por período
    double totalCost;       // Custo total 
};


struct ReinsertionData {
    int customer_id = -1; 
    
    // Dados "ANTES"
    double cost_before_removal;
    std::vector<std::vector<Route>> routes_before_removal;
    
    // Dados "DEPOIS" (sem o cliente)
    double cost_after_removal;
    std::vector<std::vector<Route>> routes_after_removal;
    Individual solution_without_customer; // Solução temporária (genótipo + fenótipo)
    
    // Dados para o Gurobi
    std::vector<GurobiPWLCurve> insertion_cost_curves; // Vetor de F_t(q_t)
    std::vector<long> max_q_inventory; // Vetor de U_i - I_{i,t-1}

    ReinsertionData(int nPeriods = 0) {
        cost_before_removal = 0.0;
        cost_after_removal = 0.0;
        routes_before_removal.resize(nPeriods);
        routes_after_removal.resize(nPeriods);
        insertion_cost_curves.resize(nPeriods);
        max_q_inventory.resize(nPeriods, 0);
    }
};



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
    int c_id_internal // ID do cliente a ser removido (0-based)
);


void busca_local(Individual& original_sol, const IRP& irp, const ACO_Params& aco_params, bool verbose);


void print_routes_for_period(const std::vector<Route>& routes, const IRP& irp, int t);

#endif 