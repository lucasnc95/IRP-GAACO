/*
 * ARQUIVO MODIFICADO: evaluation.cpp
 *
 * Implementa as novas funções 'check_feasibility' e 'calculate_total_cost'.
 */

#include "evaluation.hpp"
#include "route.hpp"
#include "irp.hpp"
#include <vector>
#include <numeric>     // Para std::accumulate
#include <set>         // Para std::set
#include <map>
#include <iostream>

using std::vector;


void check_feasibility(Individual& ind, const IRP& irp) {
    // --- Checagem 1: Validação de Rota e Capacidade (Por Período) ---
    long total_fleet_capacity = (long)irp.nVehicles * irp.Capacity;

    for (int t = 0; t < irp.nPeriods; ++t) {
        
        // 1a: Checa se o número de rotas excede a frota
        if (ind.routes_per_period[t].size() > irp.nVehicles) {
            ind.is_feasible = false;
            return; // Violação de Tamanho de Frota
        }

        std::set<int> customers_requiring_delivery;
        long total_period_delivery_plan = 0;
        
        for (int c = 0; c < irp.nCustomers; ++c) {
            if (ind.deliveries[t][c] > 0) {
                customers_requiring_delivery.insert(c + 1); // 1-based
                total_period_delivery_plan += ind.deliveries[t][c];
                
                // 1b: Checa se entrega individual excede a capacidade do veículo
                if (ind.deliveries[t][c] > irp.Capacity) {
                    ind.is_feasible = false;
                    return; 
                }
            }
        }
        
        // 1c: Checa se o plano total do dia excede a capacidade da frota
        if (total_period_delivery_plan > total_fleet_capacity) {
            ind.is_feasible = false;
            return; // Violação: Capacidade da frota excedida no plano
        }

        std::set<int> customers_in_routes;
        
        // 1d: Valida cada rota individualmente
        for (const Route& route : ind.routes_per_period[t]) {
            long actual_route_load = 0;
            
            if (!route.visits.empty()) {
                // 1e: Checa se a rota está bem formada (começa/termina no depósito)
                if (route.visits.front() != 0 || route.visits.back() != 0) {
                    ind.is_feasible = false;
                    return; // Rota malformada
                }

                for (size_t i = 1; i < route.visits.size() - 1; ++i) { // Pula depósito
                    int node = route.visits[i];
                    if (node <= 0 || node > irp.nCustomers) {
                         ind.is_feasible = false; 
                         return; // ID de cliente inválido
                    }
                    
                    customers_in_routes.insert(node);
                    actual_route_load += (long)ind.deliveries[t][node - 1]; // node é 1-based
                }

                // 1f: Checa se a carga real da rota excede a capacidade do veículo
                if (actual_route_load > irp.Capacity) {
                    ind.is_feasible = false;
                    return; // Violação de Capacidade do Veículo
                }
            }
        } // Fim da validação de rotas individuais

        // 1g: CHECAGEM DE CONSISTÊNCIA (A correção do seu bug)
        // Os clientes no plano 'deliveries' são os mesmos nas 'routes'?
        if (customers_requiring_delivery != customers_in_routes) {
            ind.is_feasible = false;
            return; // Inconsistência entre plano e rotas
        }
    } // Fim da checagem de rotas (loop 't')

    
    // --- Checagem 2: Simulação de Inventário (Clientes e Depósito) ---
    vector<long> customer_inv(irp.nCustomers);
    for (int i = 0; i < irp.nCustomers; ++i) customer_inv[i] = irp.customers[i].initialInv;
    
    long depot_inv = irp.depots[0].initialInv;

    for (int t = 0; t < irp.nPeriods; ++t) {
        depot_inv += irp.depots[0].production[t];
        long period_delivery_sum = std::accumulate(ind.deliveries[t].begin(), ind.deliveries[t].end(), 0L);

        // 2a: Checa estoque do depósito ANTES da saída
        if (period_delivery_sum > depot_inv) {
            ind.is_feasible = false;
            return; // Violação: Stock-out no Depósito
        }
        depot_inv -= period_delivery_sum;

        for (int c = 0; c < irp.nCustomers; ++c) {
            customer_inv[c] += (long)ind.deliveries[t][c];
            
            // 2b: Checa capacidade máxima do cliente (pico de estoque)
            if(customer_inv[c] > irp.customers[c].maxLevelInv) {
                ind.is_feasible = false;
                return; // Violação: Nível Máximo do Cliente
            }
            
            customer_inv[c] -= (long)irp.customers[c].demand[t];
            
            // 2c: Checa estoque mínimo (fim do dia)
            if (customer_inv[c] < irp.customers[c].minLevelInv) {
                 ind.is_feasible = false;
                 return; // Violação: Stock-out no Cliente
            }
        }
    }

    // --- SUCESSO ---
    // Se passou por todas as checagens, a solução é factível.
    ind.is_feasible = true;
    return;
}


/**
 * @brief (Função APENAS de Cálculo de Custo)
 * Calcula os custos de um indivíduo, assumindo que ele é factível.
 * Preenche 'routing_cost', 'customer_holding_cost', 'depot_holding_cost', e 'fitness'.
 * NÃO checa factibilidade e NÃO mexe no campo 'is_feasible'.
 */
void calculate_total_cost(Individual& ind, const IRP& irp) {
    
    // 1. Reseta custos
    ind.routing_cost = 0.0;
    ind.customer_holding_cost = 0.0;
    ind.depot_holding_cost = 0.0;
    ind.final_inventory_penalty = 0.0; // (Não aplicável se factível)
    
    // 2. Calcula Custo de Roteamento (lendo o custo pré-calculado da struct)
    for (const auto& period_routes : ind.routes_per_period) {
        for (const Route& route : period_routes) {
            ind.routing_cost += route.cost;
        }
    }
    
    // 3. Simula para calcular Custos de Estoque
    // (Esta simulação é segura, pois 'check_feasibility' já validou)
    vector<long> customer_inv(irp.nCustomers);
    for (int i = 0; i < irp.nCustomers; ++i) customer_inv[i] = irp.customers[i].initialInv;
    long depot_inv = irp.depots[0].initialInv;

    for (int t = 0; t < irp.nPeriods; ++t) {
        depot_inv += irp.depots[0].production[t];
        long period_delivery_sum = std::accumulate(ind.deliveries[t].begin(), ind.deliveries[t].end(), 0L);
        depot_inv -= period_delivery_sum;

        if (depot_inv > irp.depots[0].minLevelInv) {
            ind.depot_holding_cost += (double)depot_inv * irp.depots[0].invCost;
        }

        for (int c = 0; c < irp.nCustomers; ++c) {
            customer_inv[c] += (long)ind.deliveries[t][c];
            customer_inv[c] -= (long)irp.customers[c].demand[t];
            
            if (customer_inv[c] > irp.customers[c].minLevelInv) { 
                long inventory_to_cost = customer_inv[c] - irp.customers[c].minLevelInv;
                ind.customer_holding_cost += (double)(inventory_to_cost) * irp.customers[c].invCost;
            }
        }
    }
    
    // 4. Define o Fitness Final
    ind.fitness = ind.routing_cost + ind.customer_holding_cost + ind.depot_holding_cost;
}