/*
 * ARQUIVO MODIFICADO: ds_operator.hpp
 * A struct ReinsertionData foi expandida para armazenar custos e rotas
 * antes e depois da remoção do cliente, permitindo a impressão
 * detalhada dos resultados da operação.
 */

#ifndef DS_OPERATOR_HPP
#define DS_OPERATOR_HPP

#include "irp.hpp"
#include "individual.hpp"
#include "parameters.hpp"
#include <vector>

/**
 * @struct PiecewiseLinearCost
 * @brief Armazena a função de custo de inserção F_t(q_t) como breakpoints.
 * (Struct permanece como na resposta anterior)
 */
struct PiecewiseLinearCost {
    std::vector<int> q_breakpoints; 
    std::vector<double> cost_values; 
};

/**
 * @struct ReinsertionData
 * @brief Contém todos os dados necessários para o DSI de um cliente.
 *
 * Armazena o custo de inserção e a capacidade de estoque
 * para cada período do horizonte.
 */
struct ReinsertionData {
    int customer_id = -1;
    
    // --- NOVOS CAMPOS PARA LOGGING ---
    double cost_before_removal; // Custo (fitness) da solução original
    std::vector<std::vector<Route>> routes_before_removal; // Rotas da solução original
    double cost_after_removal; // Custo da solução (sem o cliente)
    std::vector<std::vector<Route>> routes_after_removal; // Rotas reroteirizadas
    // --- FIM DOS NOVOS CAMPOS ---
    
    // Dados para a DP (como antes)
    std::vector<PiecewiseLinearCost> insertion_cost_functions; 
    std::vector<long> max_q_inventory; 

    ReinsertionData(int nPeriods = 0) {
        insertion_cost_functions.resize(nPeriods);
        max_q_inventory.resize(nPeriods, 0);
        cost_before_removal = 0.0;
        cost_after_removal = 0.0;
        routes_before_removal.resize(nPeriods);
        routes_after_removal.resize(nPeriods);
    }
};


/**
 * @brief Sorteia um cliente, o remove e calcula seus dados de reinserção.
 */
ReinsertionData calculate_reinsertion_data(
    const Individual& original_ind, 
    const IRP& irp, 
    const ACO_Params& aco_params
);

void print_routes_for_period(const std::vector<Route>& routes, const IRP& irp, int t);

#endif // DS_OPERATOR_HPP