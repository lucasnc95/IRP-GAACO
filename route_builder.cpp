/*
 * ARQUIVO MODIFICADO: route_builder.cpp
 *
 * 1. Corrigido: Removida a linha errada "std::vector<int> customers_not_routed = ..."
 * 2. Corrigido: Alterado 'cust_angle.customer_id' para
 * 'cust_angle.customer_id_1based' para corresponder à struct.
 */

#include "route_builder.hpp"
#include "local_search.hpp" // Para improve_routes e recalculate_route_metrics
#include "utils.hpp"
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <map>

using std::vector;

// Estrutura auxiliar para ordenação
struct CustomerAngle {
    int customer_id_1based; // ID do cliente (1-based, como no ACO)
    double angle;

    bool operator<(const CustomerAngle& other) const {
        return angle < other.angle;
    }
};

/**
 * @brief Constrói rotas usando o algoritmo Sweep.
 */
std::vector<Route> build_routes_with_sweep(
    const IRP& irp, 
    const std::vector<int>& deliveries_for_period,
    const ACO_Params& aco_params)
{
    vector<Route> generated_routes;
    const int nCust = irp.nCustomers;
    const int depot_node = 0; // ID do depósito
    const double depot_x = irp.depots[0].x;
    const double depot_y = irp.depots[0].y;

    // --- PASSO 1: Calcular Ângulos e Lista de Clientes ---
    std::vector<CustomerAngle> customers_to_serve;
    for (int c = 0; c < nCust; ++c) {
        if (deliveries_for_period[c] > 0) {
            const auto& cust = irp.customers[c];
            // Calcula o ângulo polar em radianos
            double angle = atan2(cust.y - depot_y, cust.x - depot_x);
            customers_to_serve.push_back({c + 1, angle}); // Usa ID 1-based
        }
    }

    if (customers_to_serve.empty()) {
        return generated_routes;
    }

    std::sort(customers_to_serve.begin(), customers_to_serve.end());

    // --- PASSO 2: Construir Rotas (A "Varredura") ---
    Route current_route;
    current_route.visits.push_back(depot_node);
    current_route.remaining_capacity = irp.Capacity;
    
    // Itera sobre os clientes ordenados por ângulo
    for (const auto& cust_angle : customers_to_serve) {
        
        // CORREÇÃO DO BUG 2
        int cust_id = cust_angle.customer_id_1based; 
        
        int delivery_q = deliveries_for_period[cust_id - 1]; // 0-based

        if (delivery_q <= current_route.remaining_capacity) {
            current_route.visits.push_back(cust_id);
            current_route.remaining_capacity -= delivery_q;
        } 
        else {
            current_route.visits.push_back(depot_node);
            generated_routes.push_back(current_route);
            
            current_route = Route();
            current_route.visits.push_back(depot_node);
            current_route.visits.push_back(cust_id);
            current_route.remaining_capacity = irp.Capacity - delivery_q;
        }
    }
    
    current_route.visits.push_back(depot_node);
    generated_routes.push_back(current_route);

    // --- PASSO 3: Checagem da Frota ---
    if (generated_routes.size() > irp.nVehicles) {
        return {}; // Retorna vetor vazio (infactível)
    }

    // --- PASSO 4: Refinamento (Busca Local) ---
    
    // 4a. Construir a matriz de distâncias
    int N = 1 + nCust;
    vector<vector<double>> dist(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            int ii = (i == 0 ? 0 : irp.nDepots + (i - 1));
            int jj = (j == 0 ? 0 : irp.nDepots + (j - 1));
            dist[i][j] = irp.costMatrix[ii][jj];
        }
    }
    
    // 4b. Calcular métricas iniciais
    for(Route& route : generated_routes) {
        recalculate_route_metrics(route, irp, deliveries_for_period, dist);
    }
    
    // 4c. Chamar a Busca Local (se habilitada)
    if (aco_params.pLocalSearch > 0.0) {
        // (Assume que pLocalSearch > 0 significa 'true')
        improve_routes(generated_routes, irp, deliveries_for_period, dist, false);
    }
    
    return generated_routes;
}