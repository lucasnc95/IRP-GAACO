/*
 * ARQUIVO MODIFICADO: local_search.cpp
 * Funções de busca local refatoradas para usar a struct Route.
 * A passagem do vetor `route_loads` foi removida, pois
 * a capacidade agora é rastreada dentro de cada struct Route.
 */
#include "local_search.hpp" // <-- MUDANÇA: Inclui o novo cabeçalho
#include "irp.hpp"
#include "route.hpp"
#include <vector>
#include <numeric>
#include <iostream>
#include <algorithm>

using std::vector;

// --- Nova Função Auxiliar ---
/**
 * @brief Recalcula o custo e a capacidade restante de uma única rota.
 */
void recalculate_route_metrics(Route& route, const IRP& irp, 
                               const vector<int>& deliveries, 
                               const vector<vector<double>>& dist) 
{
    route.cost = 0.0;
    int current_load = 0;
    for (size_t i = 0; i < route.visits.size() - 1; ++i) {
        int u = route.visits[i];
        int v = route.visits[i + 1];
        route.cost += dist[u][v];
        if (i < route.visits.size() - 1 && v != 0) { // Não soma carga no retorno ao depósito
             current_load += deliveries[v - 1]; // -1 pois 'deliveries' é 0-based
        }
    }
    route.remaining_capacity = irp.Capacity - current_load;
}


// <-- MUDANÇA: Assinatura e lógica atualizadas -->
// Agora simplesmente soma os custos pré-calculados de cada rota.
static double calculate_total_cost(const vector<Route>& routes) {
    double total_cost = 0.0;
    for (const auto& route : routes) {
        total_cost += route.cost;
    }
    return total_cost;
}

// <-- MUDANÇA: Assinatura e lógica atualizadas -->
// Agora opera em uma 'Route&' e atualiza seu custo interno.
bool apply_2opt_on_route(Route& route, const vector<vector<double>>& dist, const IRP& irp, const vector<int>& deliveries) {
    if (route.visits.size() <= 4) return false;

    bool improved = false;
    bool made_improvement_in_pass = true;
    while (made_improvement_in_pass) {
        made_improvement_in_pass = false;
        for (size_t i = 1; i < route.visits.size() - 2; ++i) {
            for (size_t j = i + 1; j < route.visits.size() - 1; ++j) {
                double current_cost = dist[route.visits[i - 1]][route.visits[i]] + dist[route.visits[j]][route.visits[j + 1]];
                double new_cost = dist[route.visits[i - 1]][route.visits[j]] + dist[route.visits[i]][route.visits[j + 1]];

                if (new_cost < current_cost - 1e-9) {
                    std::reverse(route.visits.begin() + i, route.visits.begin() + j + 1);
                    made_improvement_in_pass = true;
                    improved = true;
                    goto next_pass;
                }
            }
        }
        next_pass:;
    }

    if (improved) {
        // Recalcula métricas se a rota mudou
        recalculate_route_metrics(route, irp, deliveries, dist);
    }
    return improved;
}


// <-- MUDANÇA: Assinatura e lógica atualizadas -->
// Remove 'route_loads' e usa 'route.remaining_capacity'
bool apply_inter_route_swap(vector<Route>& routes, const IRP& irp, const vector<int>& deliveries, const vector<vector<double>>& dist) {
    double best_saving = 1e-9;
    int best_r1 = -1, best_i = -1;
    int best_r2 = -1, best_j = -1;

    for (size_t r1_idx = 0; r1_idx < routes.size(); ++r1_idx) {
        for (size_t r2_idx = r1_idx + 1; r2_idx < routes.size(); ++r2_idx) {
            for (size_t i = 1; i < routes[r1_idx].visits.size() - 1; ++i) {
                for (size_t j = 1; j < routes[r2_idx].visits.size() - 1; ++j) {
                    int u = routes[r1_idx].visits[i];
                    int v = routes[r2_idx].visits[j];
                    int delivery_u = deliveries[u - 1];
                    int delivery_v = deliveries[v - 1];

                    // Verifica capacidade usando a struct
                    if ((irp.Capacity - routes[r1_idx].remaining_capacity - delivery_u + delivery_v > irp.Capacity) ||
                        (irp.Capacity - routes[r2_idx].remaining_capacity - delivery_v + delivery_u > irp.Capacity)) {
                        continue;
                    }
                    
                    int prev_u = routes[r1_idx].visits[i-1], next_u = routes[r1_idx].visits[i+1];
                    int prev_v = routes[r2_idx].visits[j-1], next_v = routes[r2_idx].visits[j+1];
                    
                    double current_cost = dist[prev_u][u] + dist[u][next_u] + dist[prev_v][v] + dist[v][next_v];
                    double new_cost = dist[prev_u][v] + dist[v][next_u] + dist[prev_v][u] + dist[u][next_v];
                    double saving = current_cost - new_cost;

                    if (saving > best_saving) {
                        best_saving = saving;
                        best_r1 = r1_idx; best_i = i;
                        best_r2 = r2_idx; best_j = j;
                    }
                }
            }
        }
    }

    if (best_r1 != -1) {
        int u = routes[best_r1].visits[best_i];
        int v = routes[best_r2].visits[best_j];
        
        // Aplica a troca
        std::swap(routes[best_r1].visits[best_i], routes[best_r2].visits[best_j]);
        
        // Recalcula métricas para ambas as rotas (mais seguro)
        recalculate_route_metrics(routes[best_r1], irp, deliveries, dist);
        recalculate_route_metrics(routes[best_r2], irp, deliveries, dist);
        
        return true;
    }
    return false;
}

// <-- MUDANÇA: Assinatura e lógica atualizadas -->
// Remove 'route_loads' e usa 'route.remaining_capacity'
bool apply_inter_route_relocate(vector<Route>& routes, const IRP& irp, const vector<int>& deliveries, const vector<vector<double>>& dist) {
    double best_saving = 1e-9;
    int best_r_from = -1, best_i = -1;
    int best_r_to = -1, best_j = -1;

    for (size_t r_from_idx = 0; r_from_idx < routes.size(); ++r_from_idx) {
        if (routes[r_from_idx].visits.size() <= 2) continue;

        for (size_t i = 1; i < routes[r_from_idx].visits.size() - 1; ++i) {
            int u = routes[r_from_idx].visits[i];
            int delivery_u = deliveries[u - 1];

            for (size_t r_to_idx = 0; r_to_idx < routes.size(); ++r_to_idx) {
                if (r_from_idx == r_to_idx) continue;

                // Verifica capacidade usando a struct
                if (routes[r_to_idx].remaining_capacity < delivery_u) continue;

                for (size_t j = 0; j < routes[r_to_idx].visits.size() - 1; ++j) {
                    int prev_u = routes[r_from_idx].visits[i-1];
                    int next_u = routes[r_from_idx].visits[i+1];
                    double cost_removed = dist[prev_u][u] + dist[u][next_u] - dist[prev_u][next_u];
                    
                    int v = routes[r_to_idx].visits[j];
                    int next_v = routes[r_to_idx].visits[j+1];
                    double cost_added = dist[v][u] + dist[u][next_v] - dist[v][next_v];
                    double saving = cost_removed - cost_added;

                    if (saving > best_saving) {
                        best_saving = saving;
                        best_r_from = r_from_idx; best_i = i;
                        best_r_to = r_to_idx; best_j = j;
                    }
                }
            }
        }
    }

    if (best_r_from != -1) {
        // Aplica a realocação
        int u = routes[best_r_from].visits[best_i];
        routes[best_r_from].visits.erase(routes[best_r_from].visits.begin() + best_i);
        routes[best_r_to].visits.insert(routes[best_r_to].visits.begin() + best_j + 1, u);

        // Recalcula métricas para ambas as rotas
        recalculate_route_metrics(routes[best_r_from], irp, deliveries, dist);
        recalculate_route_metrics(routes[best_r_to], irp, deliveries, dist);

        // Remove a rota se ficou vazia
        if (routes[best_r_from].visits.size() <= 2) {
            routes.erase(routes.begin() + best_r_from);
        }
        return true;
    }
    return false;
}

// <-- MUDANÇA: Assinatura e lógica atualizadas -->
// Agora recebe 'vector<Route>&' e não precisa de 'route_loads'
double improve_routes(vector<Route>& routes, const IRP& irp, const vector<int>& deliveries_for_period, const vector<vector<double>>& dist, bool verbose) {
    
    // Calcula o custo inicial somando os custos das rotas recebidas
    double initial_cost = calculate_total_cost(routes);

    bool improvement_found = true;
    while (improvement_found) {
        improvement_found = false;
        
        if (apply_inter_route_swap(routes, irp, deliveries_for_period, dist)) {
            improvement_found = true;
            continue;
        }
        if (apply_inter_route_relocate(routes, irp, deliveries_for_period, dist)) {
            improvement_found = true;
            continue;
        }
        
        for (Route& route : routes) {
            if (apply_2opt_on_route(route, dist, irp, deliveries_for_period)) {
                improvement_found = true;
            }
        }
    }

    double final_cost = calculate_total_cost(routes);
    
    if (verbose && final_cost < initial_cost - 1e-9) {
        std::cout << "      [Busca Local] Custo da rota melhorado: " << initial_cost << " -> " << final_cost << std::endl;
    }

    return initial_cost - final_cost;
}