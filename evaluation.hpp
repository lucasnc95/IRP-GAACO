#ifndef EVALUATION_HPP
#define EVALUATION_HPP

#include "irp.hpp"
#include "individual.hpp"
#include "parameters.hpp"

// Estrutura para conter os resultados detalhados de uma avaliação
struct EvaluationResult {
    double routing_cost = 0.0;
    double customer_holding_cost = 0.0;
    double depot_holding_cost = 0.0;
    double final_inventory_penalty = 0.0; // <-- O campo para a penalidade
    bool is_feasible = true; // Para o relatório final

    // Retorna o custo operacional (sem a penalidade final)
    double operational_cost() const {
        return routing_cost + customer_holding_cost + depot_holding_cost;
    }

    // Retorna o fitness completo (custo operacional + penalidade)
    double total_fitness() const {
        return operational_cost() + final_inventory_penalty;
    }
};

// Função central que simula um indivíduo e retorna seus custos
EvaluationResult simulate_and_evaluate(
    const Individual& ind, 
    const IRP& irp, 
    const ACO_Params& aco_params
);

#endif // EVALUATION_HPP