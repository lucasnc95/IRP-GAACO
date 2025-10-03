#include "aco.hpp"
#include "utils.hpp"
#include "local_search.hpp"
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>

using std::vector;

ACO_Result runACO_for_period(const IRP& irp,
                             const vector<int>& deliveries_for_period,
                             const ACO_Params& params,
                             bool verbose) {

    int nCust = irp.nCustomers;
    const double LARGE = 1e14;

    // Checagens iniciais de factibilidade
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

    // Monta a matriz de distância local para o ACO (0=depot, 1..N=customers)
    int depotIndex = 0;
    int baseIdx = irp.nDepots;
    int N = 1 + nCust;
    vector<vector<double>> dist(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) {
        int ii = (i==0 ? depotIndex : baseIdx + (i-1));
        int jj = (j==0 ? depotIndex : baseIdx + (j-1));
        dist[i][j] = irp.costMatrix[ii][jj];
    }

    // Inicializa matrizes de feromônio e heurística
    vector<vector<double>> tau(N, vector<double>(N, 1.0));
    vector<vector<double>> eta(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j)
        if (i != j) eta[i][j] = 1.0 / std::max(1e-9, dist[i][j]);

    double best_global_cost = LARGE;
    vector<vector<int>> best_global_routes;
    bool anyFeasibleFound = false;

    // Loop principal do ACO
    for (int iter = 0; iter < params.nIter; ++iter) {
        vector<vector<double>> delta(N, vector<double>(N, 0.0));

        // Loop das formigas
        for (int ant = 0; ant < params.nAnts; ++ant) {
            vector<int> unserved;
            for (int i = 0; i < nCust; ++i) if (deliveries_for_period[i] > 0) unserved.push_back(i+1);
            
            vector<vector<int>> antRoutes;
            double antCost = 0.0;
            bool feasibleAnt = true;

            // Constrói as rotas para a formiga atual
            for (int veh = 0; veh < irp.nVehicles && !unserved.empty(); ++veh) {
                int remainingCap = irp.Capacity;
                int cur = 0;
                vector<int> route;
                route.push_back(0);
                int safety = 0;
                const int SAFETY_LIMIT = std::max(100, nCust * 10);
                
                while (true) {
                    if (++safety > SAFETY_LIMIT) { feasibleAnt = false; break; }

                    vector<int> candidates;
                    for (int node : unserved) {
                        if (deliveries_for_period[node - 1] <= remainingCap) candidates.push_back(node);
                    }
                    if (candidates.empty()) break;

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
                    
                    route.push_back(chosenNode);
                    remainingCap -= deliveries_for_period[chosenNode - 1];
                    cur = chosenNode;
                }
                
                if (route.size() > 1) {
                    route.push_back(0);
                    antRoutes.push_back(route);
                }
            }

            if (!unserved.empty()) feasibleAnt = false;

            if (feasibleAnt && !antRoutes.empty()) {
                
                antCost = 0.0;
                for (const auto& route : antRoutes) {
                    for (size_t p = 0; p + 1 < route.size(); ++p) {
                        antCost += dist[route[p]][route[p+1]];
                    }
                }
                
                if (randreal() < params.pLocalSearch) {
                    double savings = improve_routes(antRoutes, irp, deliveries_for_period, dist, verbose);
                    antCost -= savings;
                }
                
                anyFeasibleFound = true;
                if (antCost < best_global_cost) {
                    best_global_cost = antCost;
                    best_global_routes = antRoutes;
                }
                
                double deposit = params.Q / std::max(1.0, antCost);
                for (const auto& route : antRoutes) {
                    for (size_t p = 0; p + 1 < route.size(); ++p) {
                        int u = route[p], v = route[p+1];
                        delta[u][v] += deposit;
                        delta[v][u] += deposit;
                    }
                }
            }
        } // Fim do loop das formigas

        // Evaporação e atualização do feromônio
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
        res.bestRoutes = best_global_routes;
    }
    return res;
}