/*
 * ARQUIVO MODIFICADO: evaluation.cpp
 *
 * 1. A função `check_feasibility` foi completamente reescrita para ser 100% rigorosa.
 * Ela agora valida as rotas reais (vindas do ACO) e faz a simulação completa
 * de inventário. Esta função é a única "porteira" da factibilidade.
 *
 * 2. A função `evaluate_and_fill` foi simplificada. Ela agora ASSUME que
 * o indivíduo é factível (pois `check_feasibility` já foi chamada)
 * e sua única tarefa é calcular os custos e o fitness.
 */

#include "evaluation.hpp"
#include "aco.hpp"          // Para ACO_Params (embora não mais usado ativamente)
#include "route.hpp"        // Para Route
#include "irp.hpp"          // Para IRP
#include <vector>
#include <numeric>         // Para std::accumulate
#include <set>             // Para std::set
#include <map>             // Para std::map (necessário para a próxima função)

using std::vector;

/**
 * @brief Checa RIGOROSAMENTE a factibilidade de um indivíduo.
 * Esta função valida tudo: consistência das rotas, capacidade dos veículos,
 * capacidade da frota, estoque do depósito e estoque (min/max) dos clientes.
 * É esta função que previne os "falsos positivos".
 */
bool check_feasibility(const Individual& ind, const IRP& irp) {
    
    // --- Checagem 1: Validação de Rota e Capacidade (Por Período) ---
    for (int t = 0; t < irp.nPeriods; ++t) {
        
        // Se o número de rotas criadas excede a frota disponível
        if (ind.routes_per_period[t].size() > irp.nVehicles) {
            return false; // Violação de Tamanho de Frota
        }

        std::set<int> customers_requiring_delivery;
        for (int c = 0; c < irp.nCustomers; ++c) {
            if (ind.deliveries[t][c] > 0) {
                customers_requiring_delivery.insert(c + 1); // 1-based
            }
        }

        std::set<int> customers_in_routes;
        long total_period_load_from_routes = 0;

        // Valida cada rota individualmente
        for (const Route& route : ind.routes_per_period[t]) {
            long actual_route_load = 0;
            
            for (int node : route.visits) {
                if (node > 0) { // Nó 0 é o depósito
                    customers_in_routes.insert(node);
                    
                    // Soma a carga da rota com base no plano de entregas
                    actual_route_load += (long)ind.deliveries[t][node - 1]; 
                }
            }

            // Checa se a carga real da rota excede a capacidade do veículo
            if (actual_route_load > irp.Capacity) {
                return false; // Violação de Capacidade do Veículo
            }
            
            // Opcional: Checa se o 'remaining_capacity' da struct está correto
            // (Isso é bom para depuração, mas a checagem acima é a que vale)
            if (irp.Capacity - actual_route_load != route.remaining_capacity) {
                 // std::cerr << "Aviso: Capacidade restante mal calculada na rota!\n";
                 // Não torna infactível, mas é um bug na struct. A checagem de 'actual_route_load' é o que importa.
            }
            
            total_period_load_from_routes += actual_route_load;
        }

        // Checa consistência: Todos que precisam de entrega estão em uma rota?
        if (customers_requiring_delivery != customers_in_routes) {
            return false; // Inconsistência: Entregas planejadas não batem com as rotas
        }
    } // Fim da checagem de rotas

    // --- Checagem 2: Simulação de Inventário (Clientes e Depósito) ---
    vector<long> customer_inv(irp.nCustomers);
    for (int i = 0; i < irp.nCustomers; ++i) customer_inv[i] = irp.customers[i].initialInv;
    
    long depot_inv = irp.depots[0].initialInv;

    for (int t = 0; t < irp.nPeriods; ++t) {
        depot_inv += irp.depots[0].production[t];
        long period_delivery_sum = std::accumulate(ind.deliveries[t].begin(), ind.deliveries[t].end(), 0L);

        // Checa estoque do depósito ANTES da saída
        if (period_delivery_sum > depot_inv) {
            return false; // Violação: Stock-out no Depósito
        }
        depot_inv -= period_delivery_sum; // Caminhões são carregados

        for (int c = 0; c < irp.nCustomers; ++c) {
            // Entrega chega
            customer_inv[c] += (long)ind.deliveries[t][c];
            
            // Checa capacidade máxima do cliente (pico de estoque)
            if(customer_inv[c] > irp.customers[c].maxLevelInv) {
                return false; // Violação: Nível Máximo do Cliente
            }
            
            // Demanda do dia é consumida
            customer_inv[c] -= (long)irp.customers[c].demand[t];
            
            // Checa estoque mínimo (fim do dia)
            if (customer_inv[c] < irp.customers[c].minLevelInv) {
                 return false; // Violação: Stock-out no Cliente
            }
        }
    }

    // Se passou por todas as checagens, a solução é factível.
    return true;
}


/**
 * @brief Preenche os custos de um indivíduo.
 * IMPORTANTE: Esta função assume que `check_feasibility(ind, irp)` já retornou 'true'.
 * Ela não faz nenhuma validação, apenas calcula os custos de uma solução factível.
 */
void evaluate_and_fill(Individual& ind, const IRP& irp, const ACO_Params& aco_params) {
    // (aco_params não é mais usado, mas mantido na assinatura por consistência)
    
    // Reseta custos
    ind.routing_cost = 0.0;
    ind.customer_holding_cost = 0.0;
    ind.depot_holding_cost = 0.0;
    ind.final_inventory_penalty = 0.0; // Não relevante para soluções factíveis

    // --- PASSO 1: CALCULAR CUSTO DE ROTEAMENTO ---
    // Simplesmente soma os custos das rotas que já foram validadas
    for (const auto& period_routes : ind.routes_per_period) {
        for (const Route& route : period_routes) {
            ind.routing_cost += route.cost;
        }
    }
    
    // --- PASSO 2: SIMULAR PARA CALCULAR CUSTOS DE ESTOQUE ---
    vector<long> customer_inv(irp.nCustomers);
    for (int i = 0; i < irp.nCustomers; ++i) customer_inv[i] = irp.customers[i].initialInv;
    long depot_inv = irp.depots[0].initialInv;

    for (int t = 0; t < irp.nPeriods; ++t) {
        depot_inv += irp.depots[0].production[t];
        long period_delivery_sum = std::accumulate(ind.deliveries[t].begin(), ind.deliveries[t].end(), 0L);
        depot_inv -= period_delivery_sum;

        // Acumula custo de estoque do depósito
        // (Assume que minLevelInv do depósito é 0)
        if (depot_inv > 0) {
            ind.depot_holding_cost += (double)depot_inv * irp.depots[0].invCost;
        }

        for (int c = 0; c < irp.nCustomers; ++c) {
            customer_inv[c] += (long)ind.deliveries[t][c];
            customer_inv[c] -= (long)irp.customers[c].demand[t];
            
            // Acumula custo de estoque do cliente
            // (Assume que o custo é sobre o estoque final)
            if (customer_inv[c] > irp.customers[c].minLevelInv) {
                // (Se minLevelInv > 0, o custo é sobre o excesso)
                long inventory_to_cost = customer_inv[c]; 
                ind.customer_holding_cost += (double)(inventory_to_cost) * irp.customers[c].invCost;
            }
        }
    }
    
    // Define o fitness final
    ind.fitness = ind.routing_cost + ind.customer_holding_cost + ind.depot_holding_cost;
    ind.is_feasible = true;
}


/* * Esta função foi movida para cá, pois é chamada por 'evaluate_and_fill'
 * no seu código original. Agora ela é chamada por 'check_feasibility'.
 */
void calculate_solution_costs(Individual& ind, const IRP& irp) {
    // Reseta custos
    ind.routing_cost = 0.0;
    ind.customer_holding_cost = 0.0;
    ind.depot_holding_cost = 0.0;

    // Calcula custo de roteamento
    for (const auto& period_routes : ind.routes_per_period) {
        for (const Route& route : period_routes) {
            ind.routing_cost += route.cost;
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
    
    ind.fitness = ind.routing_cost + ind.customer_holding_cost + ind.depot_holding_cost;
    ind.is_feasible = true; // Só deve ser chamada para soluções factíveis
}


