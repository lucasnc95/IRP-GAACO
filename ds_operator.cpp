/*
 * ARQUIVO MODIFICADO: ds_operator.cpp
 * * 1. Adicionada a função 'print_routes_for_period' para impressão formatada.
 * 2. 'calculate_reinsertion_data' agora salva os custos e rotas
 * originais E reroteirizados dentro da struct ReinsertionData.
 * 3. 'calculate_reinsertion_data' agora calcula o custo da solução temporária
 * (sem o cliente removido) para fins de comparação.
 */

#include "ds_operator.hpp"
#include "aco.hpp"      // Para runACO_for_period
#include "utils.hpp"    // Para randint
#include "route.hpp"
#include "evaluation.hpp" // Para calculate_solution_costs
#include <algorithm>      // Para std::sort, std::find
#include <map>
#include <set>
#include <iostream>       // Para std::cout
#include <iomanip>        // Para std::fixed, std::setprecision

using std::vector;

// --- ESTRUTURAS AUXILIARES (como na resposta anterior) ---
struct InsertionOption {
    double cost;  // Custo de inserção (gamma)
    int capacity; // Capacidade disponível (kappa)
    
    bool operator<(const InsertionOption& other) const {
        if (cost != other.cost) return cost < other.cost;
        return capacity > other.capacity;
    }
};

std::pair<int, double> find_best_insertion(const Route& route, int customer_id, const IRP& irp) {
    if (route.visits.empty() || route.visits.size() < 2) return {-1, 1e18};
    double min_delta_cost = 1e18;
    int best_pos = -1;
    for (size_t i = 0; i < route.visits.size() - 1; ++i) {
        int node_j = route.visits[i];
        int node_k = route.visits[i+1];
        double delta_cost = irp.costMatrix[node_j][customer_id] + 
                            irp.costMatrix[customer_id][node_k] - 
                            irp.costMatrix[node_j][node_k];
        if (delta_cost < min_delta_cost) {
            min_delta_cost = delta_cost;
            best_pos = i + 1;
        }
    }
    return {best_pos, min_delta_cost};
}
// --- FIM DAS ESTRUTURAS AUXILIARES ---


/**
 * @brief Implementação principal da função.
 */
ReinsertionData calculate_reinsertion_data(
    const Individual& original_ind, 
    const IRP& irp, 
    const ACO_Params& aco_params) 
{
    const double LARGE_COST = 1e18;

    // --- Ação 1: Sortear Cliente e Criar Solução Temporária ---
    int c_id_internal = randint(0, irp.nCustomers - 1); // 0-based
    int c_id_global = c_id_internal + 1; // 1-based
    
    Individual temp_sol = original_ind;
    ReinsertionData data(irp.nPeriods);
    data.customer_id = c_id_internal;

    // --- Ação 2: Salvar Dados Originais (ANTES da remoção) ---
    data.cost_before_removal = original_ind.fitness; 
    data.routes_before_removal = original_ind.routes_per_period;

    // --- Ação 3: Remover Cliente e Reroteirizar com ACO ---
    for (int t = 0; t < irp.nPeriods; ++t) {
        temp_sol.deliveries[t][c_id_internal] = 0; // Zera a entrega
        
        long period_load = std::accumulate(temp_sol.deliveries[t].begin(), temp_sol.deliveries[t].end(), 0L);
        
        if (period_load > 0) {
            ACO_Result ares = runACO_for_period(irp, temp_sol.deliveries[t], aco_params, false);
            if (ares.bestCost < LARGE_COST / 10.0) {
                temp_sol.routes_per_period[t] = ares.bestRoutes;
            } else {
                temp_sol.routes_per_period[t].clear(); 
            }
        } else {
            temp_sol.routes_per_period[t].clear(); 
        }
    }

    // --- Ação 4: Calcular Custo DEPOIS da Remoção (simulação parcial) ---
    double temp_routing_cost = 0.0;
    double temp_customer_holding_cost = 0.0;
    double temp_depot_holding_cost = 0.0;
    
    vector<long> customer_inv(irp.nCustomers);
    for (int i = 0; i < irp.nCustomers; ++i) customer_inv[i] = irp.customers[i].initialInv;
    long depot_inv = irp.depots[0].initialInv;

    for (int t = 0; t < irp.nPeriods; ++t) {
        // Custo de Rota (baseado nas rotas reroteirizadas em temp_sol)
        for(const auto& route : temp_sol.routes_per_period[t]) {
            temp_routing_cost += route.cost;
        }

        long period_delivery_sum = std::accumulate(temp_sol.deliveries[t].begin(), temp_sol.deliveries[t].end(), 0L);
        depot_inv += irp.depots[0].production[t];
        depot_inv -= period_delivery_sum;
        if (depot_inv > 0) temp_depot_holding_cost += (double)depot_inv * irp.depots[0].invCost;

        for (int c = 0; c < irp.nCustomers; ++c) {
            // Ignora o cliente removido na simulação de custo
            if (c == c_id_internal) continue; 
            
            customer_inv[c] += (long)temp_sol.deliveries[t][c];
            customer_inv[c] -= (long)irp.customers[c].demand[t];
            if (customer_inv[c] > 0) temp_customer_holding_cost += (double)customer_inv[c] * irp.customers[c].invCost;
        }
    }
    
    // Salva os dados "depois"
    data.cost_after_removal = temp_routing_cost + temp_customer_holding_cost + temp_depot_holding_cost;
    data.routes_after_removal = temp_sol.routes_per_period; 

    
    // --- Ação 5: Calcular Max Inventário e Custos de Inserção F_t(q_t) ---
    long current_inv_c = irp.customers[c_id_internal].initialInv;
    for (int t = 0; t < irp.nPeriods; ++t) {
        
        // 1. Calcula o máximo que pode ser entregue (Restrição de Inventário)
        data.max_q_inventory[t] = irp.customers[c_id_internal].maxLevelInv - current_inv_c;
        if (data.max_q_inventory[t] < 0) data.max_q_inventory[t] = 0; 

        // 2. Calcula a função de custo de inserção F_t(q_t) (Restrição de Rota/Veículo)
        vector<InsertionOption> options;

        // Opção 2a: Inserir em rotas existentes (agora usa as rotas de temp_sol)
        for (const Route& route : temp_sol.routes_per_period[t]) {
            std::pair<int, double> insertion = find_best_insertion(route, c_id_global, irp);
            if (insertion.first != -1) { 
                options.push_back({insertion.second, route.remaining_capacity});
            }
        }
        
        // Opção 2b: Criar nova rota direta
        double new_route_cost = irp.costMatrix[0][c_id_global] + irp.costMatrix[c_id_global][0];
        options.push_back({new_route_cost, (int)irp.Capacity});

        // 3. Constrói a Função de Custo Linear por Partes
        std::sort(options.begin(), options.end()); 

        PiecewiseLinearCost& ft = data.insertion_cost_functions[t];
        int max_kappa_so_far = -1;
        
        for (const auto& opt : options) {
            if (opt.capacity > max_kappa_so_far) {
                ft.q_breakpoints.push_back(opt.capacity);
                ft.cost_values.push_back(opt.cost);
                max_kappa_so_far = opt.capacity;
            }
        }
        
        if (max_kappa_so_far < (int)irp.Capacity) {
             ft.q_breakpoints.push_back((int)irp.Capacity);
             ft.cost_values.push_back(new_route_cost);
        }

        // 4. Simula o inventário de 'c' (assumindo NENHUMA entrega)
        current_inv_c -= irp.customers[c_id_internal].demand[t];
        if (current_inv_c < irp.customers[c_id_internal].minLevelInv) {
             current_inv_c = irp.customers[c_id_internal].minLevelInv;
        }
    } // Fim do loop 't'

    return data;
}


// --- NOVA FUNÇÃO HELPER DE IMPRESSÃO ---
void print_routes_for_period(const std::vector<Route>& routes, const IRP& irp, int t) {
    std::cout << "  Periodo " << t << ":\n";
    if (routes.empty()) {
        std::cout << "    sem entregas.\n";
        return;
    }
    for (size_t r = 0; r < routes.size(); ++r) {
        const Route& route = routes[r];
        std::cout << "    Rota " << r << ": ";
        for (size_t i = 0; i < route.visits.size(); ++i) {
            std::cout << route.visits[i] << (i + 1 < route.visits.size() ? " " : "");
        }
        int load = irp.Capacity - route.remaining_capacity;
        std::cout << " (Carga: " << load << "/" << irp.Capacity
                  << ", Custo: " << route.cost << ")\n";
    }
}