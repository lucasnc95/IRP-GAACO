#include "evaluation.hpp"
#include "aco.hpp"
#include <vector>
#include <numeric>
#include <iostream>

using std::vector;



EvaluationResult simulate_and_evaluate(const Individual& ind, const IRP& irp, const ACO_Params& aco_params) {
    EvaluationResult result;
    const double LARGE_COST = 1e18;
    
    vector<long> customer_inv(irp.nCustomers);
    for (int i = 0; i < irp.nCustomers; ++i) customer_inv[i] = irp.customers[i].initialInv;
    
    long depot_inv = irp.depots[0].initialInv;

    for (int t = 0; t < irp.nPeriods; ++t) {
        const auto& del = ind.deliveries[t];
        long period_delivery_sum = std::accumulate(del.begin(), del.end(), 0L);
        
        depot_inv += irp.depots[0].production[t];

        if (period_delivery_sum > depot_inv) {
            result.routing_cost = LARGE_COST;
            result.is_feasible = false;
            return result;
        }
        
        if (period_delivery_sum > 0) {
            ACO_Result ares = runACO_for_period(irp, del, aco_params, false);
            if (ares.bestCost >= LARGE_COST / 10.0) {
                result.routing_cost = LARGE_COST;
                result.is_feasible = false;
                return result;
            }
            result.routing_cost += ares.bestCost;
        }

        depot_inv -= period_delivery_sum;
        if (depot_inv > 0) {
            result.depot_holding_cost += (double)depot_inv * irp.depots[0].invCost;
        }

        for (int c = 0; c < irp.nCustomers; ++c) {
            customer_inv[c] += (long)del[c];
            if(customer_inv[c] > irp.customers[c].maxLevelInv) result.is_feasible = false;
            
            customer_inv[c] -= (long)irp.customers[c].demand[t];
            
            if (customer_inv[c] < irp.customers[c].minLevelInv) {
                 result.is_feasible = false;
            }
            
            if (customer_inv[c] > 0) {
                result.customer_holding_cost += (double)customer_inv[c] * irp.customers[c].invCost;
            }
        }
    }
    
    // Cálculo da penalidade de estoque final
    double operational_cost = result.routing_cost + result.customer_holding_cost + result.depot_holding_cost;
    long total_leftover_inventory = 0;
    for (int c = 0; c < irp.nCustomers; ++c) {
        if (customer_inv[c] > 0) {
            total_leftover_inventory += customer_inv[c];
        }
    }
 //   if (total_leftover_inventory > 0) {
   //     result.final_inventory_penalty = operational_cost * 0.10 * total_leftover_inventory;
   // }
    
    return result;
}



void evaluate_and_fill(Individual& ind, const IRP& irp, const ACO_Params& aco_params) {
    const double LARGE_COST = 1e18;
    
    // Reseta os custos antes de recalcular
    ind.routing_cost = 0.0;
    ind.customer_holding_cost = 0.0;
    ind.depot_holding_cost = 0.0;
    ind.final_inventory_penalty = 0.0;

    // --- CÁLCULO DE ROTEAMENTO (Executa o ACO se necessário) ---
    for (int t = 0; t < irp.nPeriods; ++t) {
        long period_delivery_sum = std::accumulate(ind.deliveries[t].begin(), ind.deliveries[t].end(), 0L);
        if (period_delivery_sum > 0) {
            // Só executa o ACO se as rotas para este período não foram geradas ainda
            if (ind.routes_per_period[t].empty()) {
                ACO_Result ares = runACO_for_period(irp, ind.deliveries[t], aco_params, false);
                if (ares.bestCost >= LARGE_COST / 10.0) {
                    ind.fitness = LARGE_COST; // Falha no roteamento
                    return;
                }
                ind.routes_per_period[t] = ares.bestRoutes;
                ind.routing_cost += ares.bestCost;
            } else {
                 // Se as rotas já existem, apenas recalcula o custo (muito rápido)
                 for(const auto& route : ind.routes_per_period[t]) {
                     for(size_t i = 0; i < route.size() - 1; ++i) {
                         ind.routing_cost += irp.costMatrix[route[i]][route[i+1]];
                     }
                 }
            }
        }
    }
    
    // --- CÁLCULO DOS CUSTOS DE ESTOQUE (SIMULAÇÃO) ---
    vector<long> customer_inv(irp.nCustomers);
    for (int i = 0; i < irp.nCustomers; ++i) customer_inv[i] = irp.customers[i].initialInv;
    long depot_inv = irp.depots[0].initialInv;

    for (int t = 0; t < irp.nPeriods; ++t) {
        long period_delivery_sum = std::accumulate(ind.deliveries[t].begin(), ind.deliveries[t].end(), 0L);
        
        depot_inv += irp.depots[0].production[t];
        if (depot_inv > 0) ind.depot_holding_cost += (double)depot_inv * irp.depots[0].invCost;
        depot_inv -= period_delivery_sum;

        for (int c = 0; c < irp.nCustomers; ++c) {
            customer_inv[c] += (long)ind.deliveries[t][c];
            if (customer_inv[c] > 0) ind.customer_holding_cost += (double)customer_inv[c] * irp.customers[c].invCost;
            customer_inv[c] -= (long)irp.customers[c].demand[t];
        }
    }
    
    // Penalidade de estoque final
    long total_leftover = 0;
    for (int c = 0; c < irp.nCustomers; ++c) if (customer_inv[c] > 0) total_leftover += customer_inv[c];
    if (total_leftover > 0) {
        double op_cost = ind.routing_cost + ind.customer_holding_cost + ind.depot_holding_cost;
        ind.final_inventory_penalty = op_cost * 0.10 * total_leftover;
    }

    // O fitness final é a soma de todos os componentes
    ind.fitness = ind.routing_cost + ind.customer_holding_cost + ind.depot_holding_cost + ind.final_inventory_penalty;
    ind.is_feasible = check_feasibility(ind, irp); // Verifica a factibilidade
}



bool check_feasibility(const Individual& ind, const IRP& irp) {
    
    long total_fleet_capacity = (long)irp.nVehicles * irp.Capacity;
    for (int t = 0; t < irp.nPeriods; ++t) {
        long period_load = 0;
        for (int c = 0; c < irp.nCustomers; ++c) {
            if (ind.deliveries[t][c] > irp.Capacity) {
                return false; 
            }
            period_load += ind.deliveries[t][c];
        }
        // Verifica se a soma das entregas no período excede a capacidade da frota
        if (period_load > total_fleet_capacity) {
            return false;
        }
    }

    // --- Checagem 2: Simulação de Inventário (Depósito e Clientes) ---
    vector<long> customer_inv(irp.nCustomers);
    for (int i = 0; i < irp.nCustomers; ++i) customer_inv[i] = irp.customers[i].initialInv;
    
    long depot_inv = irp.depots[0].initialInv;

    for (int t = 0; t < irp.nPeriods; ++t) {
        depot_inv += irp.depots[0].production[t];
        long period_delivery_sum = std::accumulate(ind.deliveries[t].begin(), ind.deliveries[t].end(), 0L);

        // Verifica se há estoque suficiente no depósito
        if (period_delivery_sum > depot_inv) {
            return false;
        }
        depot_inv -= period_delivery_sum;

        for (int c = 0; c < irp.nCustomers; ++c) {
            customer_inv[c] += (long)ind.deliveries[t][c];
            // Verifica se o estoque máximo do cliente foi excedido
            if(customer_inv[c] > irp.customers[c].maxLevelInv) {
                return false;
            }
            
            customer_inv[c] -= (long)irp.customers[c].demand[t];
            
            // Verifica se houve ruptura de estoque (estoque negativo)
            if (customer_inv[c] < irp.customers[c].minLevelInv) {
                 return false;
            }
        }
    }

    // Se passou por todas as checagens, a solução é factível
    return true;
}