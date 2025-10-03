#include "local_search.hpp"
#include <vector>
#include <numeric>
#include <iostream>
#include <algorithm>

using std::vector;

// Função auxiliar para calcular o custo total de um conjunto de rotas
static double calculate_total_cost(const vector<vector<int>>& routes, const vector<vector<double>>& dist) {
    double total_cost = 0.0;
    for (const auto& route : routes) {
        for (size_t i = 0; i + 1 < route.size(); ++i) {
            total_cost += dist[route[i]][route[i+1]];
        }
    }
    return total_cost;
}

// 2-Opt: já é exaustivo para uma rota, mantido como estava
bool apply_2opt_on_route(vector<int>& route, const vector<vector<double>>& dist) {
    if (route.size() <= 4) return false;

    bool improved = false;
    bool made_improvement_in_pass = true;
    while (made_improvement_in_pass) {
        made_improvement_in_pass = false;
        for (size_t i = 1; i < route.size() - 2; ++i) {
            for (size_t j = i + 1; j < route.size() - 1; ++j) {
                double current_cost = dist[route[i - 1]][route[i]] + dist[route[j]][route[j + 1]];
                double new_cost = dist[route[i - 1]][route[j]] + dist[route[i]][route[j + 1]];

                if (new_cost < current_cost - 1e-9) {
                    std::reverse(route.begin() + i, route.begin() + j + 1);
                    made_improvement_in_pass = true;
                    improved = true;
                    goto next_pass;
                }
            }
        }
        next_pass:;
    }
    return improved;
}

// <-- MUDANÇA: Swap agora implementa "Best Improvement" -->
// A função agora procura a melhor troca possível e a aplica no final.
bool apply_inter_route_swap(vector<vector<int>>& routes, vector<int>& route_loads, const IRP& irp, const vector<int>& deliveries, const vector<vector<double>>& dist) {
    double best_saving = 1e-9;
    int best_r1 = -1, best_i = -1;
    int best_r2 = -1, best_j = -1;

    for (size_t r1_idx = 0; r1_idx < routes.size(); ++r1_idx) {
        for (size_t r2_idx = r1_idx + 1; r2_idx < routes.size(); ++r2_idx) {
            for (size_t i = 1; i < routes[r1_idx].size() - 1; ++i) {
                for (size_t j = 1; j < routes[r2_idx].size() - 1; ++j) {
                    int u = routes[r1_idx][i];
                    int v = routes[r2_idx][j];
                    int delivery_u = deliveries[u - 1];
                    int delivery_v = deliveries[v - 1];

                    if ((route_loads[r1_idx] - delivery_u + delivery_v > irp.Capacity) ||
                        (route_loads[r2_idx] - delivery_v + delivery_u > irp.Capacity)) {
                        continue;
                    }
                    
                    int prev_u = routes[r1_idx][i-1], next_u = routes[r1_idx][i+1];
                    int prev_v = routes[r2_idx][j-1], next_v = routes[r2_idx][j+1];
                    
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
        int u = routes[best_r1][best_i];
        int v = routes[best_r2][best_j];
        int delivery_u = deliveries[u - 1];
        int delivery_v = deliveries[v - 1];

        std::swap(routes[best_r1][best_i], routes[best_r2][best_j]);
        route_loads[best_r1] = route_loads[best_r1] - delivery_u + delivery_v;
        route_loads[best_r2] = route_loads[best_r2] - delivery_v + delivery_u;
        return true;
    }
    return false;
}

// <-- MUDANÇA: Relocate agora implementa "Best Improvement" -->
// A função agora procura a melhor realocação possível e a aplica no final.
bool apply_inter_route_relocate(vector<vector<int>>& routes, vector<int>& route_loads, const IRP& irp, const vector<int>& deliveries, const vector<vector<double>>& dist) {
    double best_saving = 1e-9;
    int best_r_from = -1, best_i = -1;
    int best_r_to = -1, best_j = -1;

    for (size_t r_from_idx = 0; r_from_idx < routes.size(); ++r_from_idx) {
        if (routes[r_from_idx].size() <= 2) continue;

        for (size_t i = 1; i < routes[r_from_idx].size() - 1; ++i) {
            int u = routes[r_from_idx][i];
            int delivery_u = deliveries[u - 1];

            for (size_t r_to_idx = 0; r_to_idx < routes.size(); ++r_to_idx) {
                if (r_from_idx == r_to_idx) continue;

                if (route_loads[r_to_idx] + delivery_u > irp.Capacity) continue;

                for (size_t j = 0; j < routes[r_to_idx].size() - 1; ++j) {
                    int prev_u = routes[r_from_idx][i-1];
                    int next_u = routes[r_from_idx][i+1];
                    double cost_removed = dist[prev_u][u] + dist[u][next_u] - dist[prev_u][next_u];
                    
                    int v = routes[r_to_idx][j];
                    int next_v = routes[r_to_idx][j+1];
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
        int u = routes[best_r_from][best_i];
        int delivery_u = deliveries[u-1];

        routes[best_r_from].erase(routes[best_r_from].begin() + best_i);
        routes[best_r_to].insert(routes[best_r_to].begin() + best_j + 1, u);

        route_loads[best_r_from] -= delivery_u;
        route_loads[best_r_to] += delivery_u;

        if (routes[best_r_from].size() <= 2) {
            routes.erase(routes.begin() + best_r_from);
            route_loads.erase(route_loads.begin() + best_r_from);
        }
        return true;
    }
    return false;
}

// Orquestrador principal, que continua a chamar os operadores até que não haja mais melhorias.
double improve_routes(vector<vector<int>>& routes, const IRP& irp, const vector<int>& deliveries_for_period, const vector<vector<double>>& dist, bool verbose) {
    double initial_cost = calculate_total_cost(routes, dist);
    
    // Calcula as cargas iniciais de cada rota
    vector<int> route_loads(routes.size(), 0);
    for(size_t i = 0; i < routes.size(); ++i) {
        for(size_t j = 1; j < routes[i].size() - 1; ++j) {
            int customer_node = routes[i][j];
            route_loads[i] += deliveries_for_period[customer_node - 1];
        }
    }

    bool improvement_found = true;
    while (improvement_found) {
        improvement_found = false;
        
        if (apply_inter_route_swap(routes, route_loads, irp, deliveries_for_period, dist)) {
            improvement_found = true;
            continue;
        }
        if (apply_inter_route_relocate(routes, route_loads, irp, deliveries_for_period, dist)) {
            improvement_found = true;
            continue;
        }
        // 2-Opt é aplicado por último, pois é intra-rota e não afeta outras rotas
        for (auto& route : routes) {
            if (apply_2opt_on_route(route, dist)) {
                improvement_found = true;
            }
        }
    }

    double final_cost = calculate_total_cost(routes, dist);
    
    if (verbose && final_cost < initial_cost - 1e-9) {
        std::cout << "      [Busca Local] Custo da rota melhorado: " << initial_cost << " -> " << final_cost << std::endl;
    }

    return initial_cost - final_cost;
}