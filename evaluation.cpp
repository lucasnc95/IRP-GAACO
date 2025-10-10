#include "evaluation.hpp"
#include "aco.hpp"
#include <vector>
#include <numeric>

using std::vector;

void evaluate_and_fill(Individual& ind, const IRP& irp, const ACO_Params& aco_params) {
    const double LARGE_COST = 1e18;
    
    ind.routing_cost = 0.0;
    ind.customer_holding_cost = 0.0;
    ind.depot_holding_cost = 0.0;
    ind.final_inventory_penalty = 0.0; // Zera a penalidade
    ind.routes_per_period.assign(irp.nPeriods, {});

    for (int t = 0; t < irp.nPeriods; ++t) {
        long period_delivery_sum = std::accumulate(ind.deliveries[t].begin(), ind.deliveries[t].end(), 0L);
        if (period_delivery_sum > 0) {
            ACO_Result ares = runACO_for_period(irp, ind.deliveries[t], aco_params, false);
            if (ares.bestCost >= LARGE_COST / 10.0) {
                ind.fitness = LARGE_COST;
                ind.is_feasible = false;
                return;
            }
            ind.routes_per_period[t] = ares.bestRoutes;
            ind.routing_cost += ares.bestCost;
        }
    }
    
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
    
    // O fitness é simplesmente o custo operacional
    ind.fitness = ind.routing_cost + ind.customer_holding_cost + ind.depot_holding_cost;
    ind.is_feasible = check_feasibility(ind, irp);
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