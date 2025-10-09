#ifndef EVALUATION_HPP
#define EVALUATION_HPP

#include "irp.hpp"
#include "individual.hpp"
#include "parameters.hpp"

// Estrutura para conter os resultados de custo de uma solução
struct EvaluationResult {
    double routing_cost = 0.0;
    double customer_holding_cost = 0.0;
    double depot_holding_cost = 0.0;
    double final_inventory_penalty = 0.0;
    bool is_feasible = true;

    // Retorna o custo operacional (custos reais incorridos durante os períodos)
    double operational_cost() const {
        return routing_cost + customer_holding_cost + depot_holding_cost;
    }

    // Retorna o fitness completo (custo operacional + penalidade para guiar o GA)
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

// Função dedicada para checar a factibilidade
bool check_feasibility(const Individual& ind, const IRP& irp);

#endif // EVALUATION_HPP