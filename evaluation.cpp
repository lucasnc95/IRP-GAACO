#include "evaluation.hpp"
#include "aco.hpp"
#include <vector>
#include <numeric>
#include <set>

using std::vector;







void evaluate_and_fill(Individual& ind, const IRP& irp, const ACO_Params& aco_params) {
    const double LARGE_COST = 1e18;
    
    // Reseta todos os campos de resultado antes de recalcular
    ind.routing_cost = 0.0;
    ind.customer_holding_cost = 0.0;
    ind.depot_holding_cost = 0.0;
    ind.final_inventory_penalty = 0.0;
    ind.routes_per_period.assign(irp.nPeriods, {});
    ind.is_feasible = true; // Assume factibilidade até que uma violação seja encontrada

    // --- PASSO 1: GERAR ROTAS E VALIDAR ---
    for (int t = 0; t < irp.nPeriods; ++t) {
        long period_delivery_sum = std::accumulate(ind.deliveries[t].begin(), ind.deliveries[t].end(), 0L);
        if (period_delivery_sum > 0) {
            ACO_Result ares = runACO_for_period(irp, ind.deliveries[t], aco_params, false);

            // Se o ACO não encontrou solução, a rota é infactível
            if (ares.bestCost >= LARGE_COST / 10.0) {
                ind.fitness = LARGE_COST;
                ind.is_feasible = false;
                return;
            }

            // --- VERIFICAÇÃO DE CLIENTES NÃO ATENDIDOS ---
            std::set<int> customers_requiring_delivery;
            for (int c = 0; c < irp.nCustomers; ++c) {
                if (ind.deliveries[t][c] > 0) {
                    customers_requiring_delivery.insert(c + 1); // IDs no ACO são 1-based
                }
            }

            std::set<int> customers_in_routes;
            for (const auto& route : ares.bestRoutes) {
                for (int node : route) {
                    if (node > 0) { // Ignora o depósito (nó 0)
                        customers_in_routes.insert(node);
                    }
                }
            }

            // Se a lista de clientes que precisam de entrega não for igual à lista de clientes nas rotas, o ACO falhou.
            if (customers_requiring_delivery != customers_in_routes) {
                ind.fitness = LARGE_COST;
                ind.is_feasible = false;
                return; // Falha crítica, a solução é infactível
            }
            // --- FIM DA VERIFICAÇÃO ---

            ind.routes_per_period[t] = ares.bestRoutes;
        }
    }
    
    // --- PASSO 2: CALCULAR CUSTOS (SE AS ROTAS FORAM GERADAS COM SUCESSO) ---
    
    // Calcula o custo de roteamento a partir das rotas já geradas e validadas
    for (const auto& period_routes : ind.routes_per_period) {
        for (const auto& route : period_routes) {
            for (size_t i = 0; i < route.size() - 1; ++i) {
                ind.routing_cost += irp.costMatrix[route[i]][route[i+1]];
            }
        }
    }

    // Simula para calcular custos de estoque
    vector<long> customer_inv(irp.nCustomers);
    for (int i = 0; i < irp.nCustomers; ++i) customer_inv[i] = irp.customers[i].initialInv;
    long depot_inv = irp.depots[0].initialInv;

    for (int t = 0; t < irp.nPeriods; ++t) {
        long period_delivery_sum = std::accumulate(ind.deliveries[t].begin(), ind.deliveries[t].end(), 0L);
        depot_inv += irp.depots[0].production[t];
        depot_inv -= period_delivery_sum;
        if (depot_inv > 0) ind.depot_holding_cost += (double)depot_inv * irp.depots[0].invCost;

        for (int c = 0; c < irp.nCustomers; ++c) {
            customer_inv[c] += (long)ind.deliveries[t][c];
            customer_inv[c] -= (long)irp.customers[c].demand[t];
            if (customer_inv[c] > 0) ind.customer_holding_cost += (double)customer_inv[c] * irp.customers[c].invCost;
        }
    }
    
    // Se a simulação passou, mas a checagem geral falha, ainda é infactível
    if (!check_feasibility(ind, irp)) {
        ind.is_feasible = false;
        ind.fitness = LARGE_COST;
        return;
    }
    
    // Se tudo estiver OK, calcula o fitness final
    ind.fitness = ind.routing_cost + ind.customer_holding_cost + ind.depot_holding_cost;
}






// void evaluate_and_fill(Individual& ind, const IRP& irp, const ACO_Params& aco_params) {
//     const double LARGE_COST = 1e18;
    
//     ind.routing_cost = 0.0;
//     ind.customer_holding_cost = 0.0;
//     ind.depot_holding_cost = 0.0;
//     ind.final_inventory_penalty = 0.0;
//     ind.routes_per_period.assign(irp.nPeriods, {});

//     for (int t = 0; t < irp.nPeriods; ++t) {
//         long period_delivery_sum = std::accumulate(ind.deliveries[t].begin(), ind.deliveries[t].end(), 0L);
//         if (period_delivery_sum > 0) {
//             ACO_Result ares = runACO_for_period(irp, ind.deliveries[t], aco_params, false);
//             if (ares.bestCost >= LARGE_COST / 10.0) {
//                 ind.fitness = LARGE_COST;
//                 ind.is_feasible = false;
//                 return;
//             }
//             ind.routes_per_period[t] = ares.bestRoutes;
//         }
//     }
    
//     // Calcula o custo de roteamento a partir das rotas geradas
//     for (const auto& period_routes : ind.routes_per_period) {
//         for (const auto& route : period_routes) {
//             for (size_t i = 0; i < route.size() - 1; ++i) {
//                 ind.routing_cost += irp.costMatrix[route[i]][route[i+1]];
//             }
//         }
//     }

//     // Simula para calcular custos de estoque
//     vector<long> customer_inv(irp.nCustomers);
//     for (int i = 0; i < irp.nCustomers; ++i) customer_inv[i] = irp.customers[i].initialInv;
//     long depot_inv = irp.depots[0].initialInv;

//     for (int t = 0; t < irp.nPeriods; ++t) {
//         long period_delivery_sum = std::accumulate(ind.deliveries[t].begin(), ind.deliveries[t].end(), 0L);
//         depot_inv += irp.depots[0].production[t];
//         depot_inv -= period_delivery_sum;
//         if (depot_inv > 0) ind.depot_holding_cost += (double)depot_inv * irp.depots[0].invCost;

//         for (int c = 0; c < irp.nCustomers; ++c) {
//             customer_inv[c] += (long)ind.deliveries[t][c];
//             customer_inv[c] -= (long)irp.customers[c].demand[t];
//             if (customer_inv[c] > 0) ind.customer_holding_cost += (double)customer_inv[c] * irp.customers[c].invCost;
//         }
//     }
    
//     // Penalidade de estoque final
//     double operational_cost = ind.routing_cost + ind.customer_holding_cost + ind.depot_holding_cost;
//     long total_leftover_inventory = 0;
//     for (int c = 0; c < irp.nCustomers; ++c) {
//         if (customer_inv[c] > 0) {
//             total_leftover_inventory += customer_inv[c];
//         }
//     }
//     if (total_leftover_inventory > 0) {
//         ind.final_inventory_penalty = operational_cost * 0.10 * total_leftover_inventory;
//     }
    
//     ind.fitness = operational_cost;
//     ind.is_feasible = check_feasibility(ind, irp);
// }

// bool check_feasibility(const Individual& ind, const IRP& irp) {
    
//     long total_fleet_capacity = (long)irp.nVehicles * irp.Capacity;
//     for (int t = 0; t < irp.nPeriods; ++t) {
//         long period_load = 0;
//         for (int c = 0; c < irp.nCustomers; ++c) {
//             if (ind.deliveries[t][c] > irp.Capacity) {
//                 return false; 
//             }
//             period_load += ind.deliveries[t][c];
//         }
//         // Verifica se a soma das entregas no período excede a capacidade da frota
//         if (period_load > total_fleet_capacity) {
//             return false;
//         }
//     }

//     // --- Checagem 2: Simulação de Inventário (Depósito e Clientes) ---
//     vector<long> customer_inv(irp.nCustomers);
//     for (int i = 0; i < irp.nCustomers; ++i) customer_inv[i] = irp.customers[i].initialInv;
    
//     long depot_inv = irp.depots[0].initialInv;

//     for (int t = 0; t < irp.nPeriods; ++t) {
//         depot_inv += irp.depots[0].production[t];
//         long period_delivery_sum = std::accumulate(ind.deliveries[t].begin(), ind.deliveries[t].end(), 0L);

//         // Verifica se há estoque suficiente no depósito
//         if (period_delivery_sum > depot_inv) {
//             return false;
//         }
//         depot_inv -= period_delivery_sum;

//         for (int c = 0; c < irp.nCustomers; ++c) {
//             customer_inv[c] += (long)ind.deliveries[t][c];
//             // Verifica se o estoque máximo do cliente foi excedido
//             if(customer_inv[c] > irp.customers[c].maxLevelInv) {
//                 return false;
//             }
            
//             customer_inv[c] -= (long)irp.customers[c].demand[t];
            
//             // Verifica se houve ruptura de estoque (estoque negativo)
//             if (customer_inv[c] < irp.customers[c].minLevelInv) {
//                  return false;
//             }
//         }
//     }

//     // Se passou por todas as checagens, a solução é factível
//     return true;
// }


// bool check_feasibility(const Individual& ind, const IRP& irp) {
    
//     // --- Checagem 1: Capacidade de Veículos e Frota ---
//     long total_fleet_capacity = (long)irp.nVehicles * irp.Capacity;
//     for (int t = 0; t < irp.nPeriods; ++t) {
//         long period_load = 0;
//         for (int c = 0; c < irp.nCustomers; ++c) {
//             // Verifica se uma única entrega excede a capacidade de um veículo
//             if (ind.deliveries[t][c] > irp.Capacity) {
//                 return false; 
//             }
//             period_load += ind.deliveries[t][c];
//         }
//         // Verifica se a soma das entregas no período excede a capacidade da frota
//         if (period_load > total_fleet_capacity) {
//             return false;
//         }
//     }

//     // --- Checagem 2: Simulação de Inventário (Depósito e Clientes) ---
//     vector<long> customer_inv(irp.nCustomers);
//     for (int i = 0; i < irp.nCustomers; ++i) customer_inv[i] = irp.customers[i].initialInv;
    
//     long depot_inv = irp.depots[0].initialInv;

//     for (int t = 0; t < irp.nPeriods; ++t) {
//         depot_inv += irp.depots[0].production[t];
//         long period_delivery_sum = std::accumulate(ind.deliveries[t].begin(), ind.deliveries[t].end(), 0L);

//         // Verifica se há estoque suficiente no depósito para as entregas do período
//         if (period_delivery_sum > depot_inv) {
//             return false;
//         }
//         depot_inv -= period_delivery_sum;

//         for (int c = 0; c < irp.nCustomers; ++c) {
//             // O estoque no início do período 't' é o que sobrou do período 't-1'
//             // 1. A entrega chega
//             customer_inv[c] += (long)ind.deliveries[t][c];
            
//             // 2. A demanda do período é consumida
//             customer_inv[c] -= (long)irp.customers[c].demand[t];
            
//             // <-- MUDANÇA: As checagens de inventário são feitas AGORA, sobre o estoque final do período -->
//             // Verifica se houve ruptura de estoque (abaixo do mínimo)
//             if (customer_inv[c] < irp.customers[c].minLevelInv) {
//                  return false;
//             }
            
//             // Verifica se o estoque máximo foi excedido
//             if(customer_inv[c] > irp.customers[c].maxLevelInv) {
//                 return false;
//             }
//         }
//     }

//     // Se passou por todas as checagens, a solução é factível
//     return true;
// }

bool check_feasibility(const Individual& ind, const IRP& irp) {
    
    // --- Checagem 1: Capacidade de Veículos e Frota ---
    long total_fleet_capacity = (long)irp.nVehicles * irp.Capacity;
    for (int t = 0; t < irp.nPeriods; ++t) {
        long period_load = 0;
        // Mapeia alocação em veículos discretos
        vector<long> vehicle_loads(irp.nVehicles, 0);
        
        for (int c = 0; c < irp.nCustomers; ++c) {
            int delivery_q = ind.deliveries[t][c];
            if (delivery_q == 0) continue;

            if (delivery_q > irp.Capacity) {
                return false; // Entrega individual excede capacidade de um veículo
            }
            
            // Tenta alocar a entrega a um veículo
            bool allocated = false;
            for(int v=0; v < irp.nVehicles; ++v) {
                if(vehicle_loads[v] + delivery_q <= irp.Capacity) {
                    vehicle_loads[v] += delivery_q;
                    allocated = true;
                    break;
                }
            }
            
            if(!allocated) {
                return false; // Não há veículos suficientes para alocar as entregas individuais
            }
            
            period_load += delivery_q;
        }
        
        if (period_load > total_fleet_capacity) {
            return false; // Carga total do período excede capacidade da frota
        }
    }

    // --- Checagem 2: Simulação de Inventário (Depósito e Clientes) ---
    vector<long> customer_inv(irp.nCustomers);
    for (int i = 0; i < irp.nCustomers; ++i) customer_inv[i] = irp.customers[i].initialInv;
    
    long depot_inv = irp.depots[0].initialInv;

    for (int t = 0; t < irp.nPeriods; ++t) {
        depot_inv += irp.depots[0].production[t];
        long period_delivery_sum = std::accumulate(ind.deliveries[t].begin(), ind.deliveries[t].end(), 0L);

        if (period_delivery_sum > depot_inv) {
            return false;
        }
        depot_inv -= period_delivery_sum;

        for (int c = 0; c < irp.nCustomers; ++c) {
            customer_inv[c] += (long)ind.deliveries[t][c];
            
            // VERIFICAÇÃO DE PICO DE ESTOQUE (antes da demanda)
            if(customer_inv[c] > irp.customers[c].maxLevelInv) {
                return false;
            }
            
            customer_inv[c] -= (long)irp.customers[c].demand[t];
            
            // VERIFICAÇÃO DE ESTOQUE MÍNIMO (após a demanda)
            if (customer_inv[c] < irp.customers[c].minLevelInv) {
                 return false;
            }
        }
    }

    return true;
}