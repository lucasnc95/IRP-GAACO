/*
 * ARQUIVO MODIFICADO: aco.cpp
 * Reescrito para construir e avaliar a nova struct Route.
 */
#include "aco.hpp"
#include "utils.hpp"
#include "local_search.hpp"
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>

using std::vector;

// --- Função auxiliar interna para calcular o custo de uma rota ---
// (Necessária para a construção inicial da rota pela formiga)
static double calculate_route_cost_helper(const vector<int>& visits, const vector<vector<double>>& dist) {
    double cost = 0.0;
    for (size_t i = 0; i < visits.size() - 1; ++i) {
        cost += dist[visits[i]][visits[i+1]];
    }
    return cost;
}


ACO_Result runACO_for_period(const IRP& irp,
                             const vector<int>& deliveries_for_period,
                             const ACO_Params& params,
                             bool verbose) {

    int nCust = irp.nCustomers;
    const double LARGE = 1e14;

    // Checagens iniciais de factibilidade (sem mudanças)
    long totalCapacity = (long)irp.nVehicles * irp.Capacity;
    long sumDeliveries = 0;
    for (int c = 0; c < nCust; ++c) sumDeliveries += deliveries_for_period[c];
    if (sumDeliveries > totalCapacity) {
        ACO_Result bad; bad.bestCost = LARGE; bad.bestRoutes.clear();
        return bad;
    }
    for (int c = 0; c < nCust; ++c) {
        if (deliveries_for_period[c] > irp.Capacity) {
            ACO_Result bad; bad.bestCost = LARGE; bad.bestRoutes.clear();
            return bad;
        }
    }

    // Matriz de distância local (sem mudanças)
    int depotIndex = 0;
    int baseIdx = irp.nDepots;
    int N = 1 + nCust;
    vector<vector<double>> dist(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) {
        int ii = (i==0 ? depotIndex : baseIdx + (i-1));
        int jj = (j==0 ? depotIndex : baseIdx + (j-1));
        dist[i][j] = irp.costMatrix[ii][jj];
    }

    // Matrizes de feromônio e heurística (sem mudanças)
    vector<vector<double>> tau(N, vector<double>(N, 1.0));
    vector<vector<double>> eta(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j)
        if (i != j) eta[i][j] = 1.0 / std::max(1e-9, dist[i][j]);

    double best_global_cost = LARGE;
    
    // <-- MUDANÇA: Agora é um vetor de structs Route -->
    vector<Route> best_global_routes;
    
    bool anyFeasibleFound = false;

    // Loop principal do ACO
    for (int iter = 0; iter < params.nIter; ++iter) {
        vector<vector<double>> delta(N, vector<double>(N, 0.0));

        // Loop das formigas
        for (int ant = 0; ant < params.nAnts; ++ant) {
            vector<int> unserved;
            for (int i = 0; i < nCust; ++i) if (deliveries_for_period[i] > 0) unserved.push_back(i+1);
            
            // <-- MUDANÇA: Agora é um vetor de structs Route -->
            vector<Route> antRoutes;
            double antCost = 0.0;
            bool feasibleAnt = true;

            // Constrói as rotas para a formiga atual
            for (int veh = 0; veh < irp.nVehicles && !unserved.empty(); ++veh) {
                
                // <-- MUDANÇA: Cria e inicializa a struct Route -->
                Route current_route;
                current_route.remaining_capacity = irp.Capacity;
                current_route.visits.push_back(0); // Começa no depósito
                
                int cur = 0;
                int safety = 0;
                const int SAFETY_LIMIT = std::max(100, nCust * 10);
                
                while (true) {
                    if (++safety > SAFETY_LIMIT) { feasibleAnt = false; break; }

                    vector<int> candidates;
                    for (int node : unserved) {
                        // <-- MUDANÇA: Checa a capacidade restante da struct -->
                        if (deliveries_for_period[node - 1] <= current_route.remaining_capacity) {
                            candidates.push_back(node);
                        }
                    }
                    if (candidates.empty()) break; // Nenhum cliente cabe, fecha a rota

                    // Lógica de seleção de candidatos (sem mudanças)
                    double sum = 0.0;
                    vector<double> weights(candidates.size());
                    for (size_t k = 0; k < candidates.size(); ++k) {
                        int node = candidates[k];
                        double val = pow(tau[cur][node], params.alpha) * pow(eta[cur][node], params.beta);
                        weights[k] = val; sum += val;
                    }
                    
                    int chosenNode;
                    if (sum <= 1e-9) {
                        double bestD = 1e18; int bestNode = candidates[0];
                        for (int node : candidates) if (dist[cur][node] < bestD) { bestD = dist[cur][node]; bestNode = node; }
                        chosenNode = bestNode;
                    } else {
                        double r = randreal() * sum;
                        double acc = 0.0;
                        chosenNode = candidates.back();
                        for (size_t k = 0; k < candidates.size(); ++k) {
                            acc += weights[k];
                            if (r <= acc) { chosenNode = candidates[k]; break; }
                        }
                    }
                    
                    auto it = find(unserved.begin(), unserved.end(), chosenNode);
                    if (it != unserved.end()) unserved.erase(it);
                    
                    // <-- MUDANÇA: Atualiza a struct Route -->
                    current_route.visits.push_back(chosenNode);
                    current_route.remaining_capacity -= deliveries_for_period[chosenNode - 1];
                    cur = chosenNode;
                }
                
                if (!feasibleAnt) break;

                // <-- MUDANÇA: Finaliza e armazena a struct Route -->
                if (current_route.visits.size() > 1) { // Se visitou alguém
                    current_route.visits.push_back(0); // Retorna ao depósito
                    // Calcula o custo inicial da rota
                    current_route.cost = calculate_route_cost_helper(current_route.visits, dist);
                    antRoutes.push_back(current_route);
                }
            } // Fim do loop de veículos

            if (!unserved.empty()) feasibleAnt = false;

            if (feasibleAnt && !antRoutes.empty()) {
                
                // <-- MUDANÇA: A busca local agora opera em 'vector<Route>' -->
                if (randreal() < params.pLocalSearch) {
                    improve_routes(antRoutes, irp, deliveries_for_period, dist, verbose);
                }
                
                // <-- MUDANÇA: Custo total é somado dos custos da struct -->
                antCost = 0.0;
                for (const auto& route : antRoutes) {
                    antCost += route.cost;
                }
                
                anyFeasibleFound = true;
                if (antCost < best_global_cost) {
                    best_global_cost = antCost;
                    best_global_routes = antRoutes;
                }
                
                // Depósito de feromônio (lógica interna sem mudanças)
                double deposit = params.Q / std::max(1.0, antCost);
                for (const auto& route : antRoutes) {
                    for (size_t p = 0; p + 1 < route.visits.size(); ++p) {
                        int u = route.visits[p], v = route.visits[p+1];
                        delta[u][v] += deposit;
                        delta[v][u] += deposit;
                    }
                }
            }
        } // Fim do loop das formigas

        // Evaporação (sem mudanças)
        for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) {
            tau[i][j] = (1.0 - params.rho) * tau[i][j] + delta[i][j];
            if (tau[i][j] < 1e-12) tau[i][j] = 1e-12;
        }
    } // Fim do loop de iterações do ACO
    
    ACO_Result res;
    if (!anyFeasibleFound) {
        res.bestCost = LARGE;
        res.bestRoutes.clear();
    } else {
        res.bestCost = best_global_cost;
        res.bestRoutes = best_global_routes; // Retorna o vetor<Route>
    }
    return res;
}