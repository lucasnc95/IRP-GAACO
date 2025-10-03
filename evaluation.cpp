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
            double period_depot_cost = (double)depot_inv * irp.depots[0].invCost;
            
            // --- LINHA DE DEPURAÇÃO ATIVA ---
         /*   std::cout << "DEBUG -> Período " << t 
                      << ": DepotInv Fim=" << depot_inv 
                      << ", InvCost=" << irp.depots[0].invCost 
                      << ", Custo Adicionado=" << period_depot_cost << std::endl;*/
            
            result.depot_holding_cost += period_depot_cost;
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
    if (total_leftover_inventory > 0) {
        result.final_inventory_penalty = operational_cost * 0.10 * total_leftover_inventory;
    }
    
    // A linha "result.final_fitness = result.total_fitness();" foi removida daqui.
    
    return result;
}