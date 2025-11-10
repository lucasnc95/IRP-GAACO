#include "ds_operator.hpp"
#include "aco.hpp"      // Para runACO_for_period
#include "utils.hpp"    // Para randint
#include "route.hpp"
#include "evaluation.hpp" // Para calculate_total_cost, build_routes_for_individual
#include "ga.hpp"           // Para 'build_routes_for_individual'
#include <algorithm>      // Para std::sort, std::find
#include <map>
#include <set>
#include <iostream>       // Para std::cout
#include <iomanip>        // Para std::fixed, std::setprecision

using std::vector;


const double PENALTY_COST = 1e7;

struct InsertionOption {
    double cost;  
    int capacity; 
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





GurobiInventoryResult computePerfectInventory_Gurobi(
    const Customer& cust,
    const Depot& depot,
    const std::vector<double>& reinsertion_cost_fixed,
    const std::vector<GurobiPWLCurve>& reinsertion_curves,
    int T, double qmax) // qmax é a capacidade do VEÍCULO
{
    try {
        GRBEnv env(true);
        env.set("LogToConsole", "0"); 
        env.set(GRB_IntParam_Threads, 0); 
        env.set(GRB_IntParam_NumericFocus, 3); 
        env.start();
        GRBModel model(env);

        model.set(GRB_DoubleParam_MIPGap, 0.0);
        
        std::vector<GRBVar> q(T), S(T);
        for (int t = 0; t < T; ++t) {
            q[t] = model.addVar(0.0, qmax, 0.0, GRB_CONTINUOUS, "q_" + std::to_string(t+1));
            S[t] = model.addVar(cust.minLevelInv, cust.maxLevelInv, 0.0, GRB_CONTINUOUS, "S_" + std::to_string(t+1));
        }

        for (int t = 0; t < T; ++t) {
            GRBLinExpr bal_lhs, peak_lhs;
            if (t == 0) {
                bal_lhs = cust.initialInv + q[t] - cust.demand[t];
                peak_lhs = cust.initialInv + q[t];
            } else {
                bal_lhs = S[t-1] + q[t] - cust.demand[t];
                peak_lhs = S[t-1] + q[t];
            }
            model.addConstr(S[t] == bal_lhs, "bal_" + std::to_string(t+1));
            model.addConstr(peak_lhs <= cust.maxLevelInv, "peak_inv_constr_" + std::to_string(t+1));
        }

        // --- INÍCIO DA MODIFICAÇÃO (FUNÇÃO OBJETIVO) ---
        GRBLinExpr obj = 0.0;
        for (int t = 0; t < T; ++t) {
            // Custo de manter estoque no CLIENTE
            obj += cust.invCost * S[t];     
            
            // Custo de transportar para o cliente (Frete)
            if (!reinsertion_curves.empty() && !reinsertion_curves[t].q_points.empty()) {
                GRBVar c_q = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "c_q_" + std::to_string(t+1));
                const auto& curve = reinsertion_curves[t];
                model.addGenConstrPWL(q[t], c_q,
                    (int)curve.q_points.size(),
                    curve.q_points.data(), curve.cost_points.data(),
                    "pwl_reinsert_t" + std::to_string(t+1));
                obj += c_q; 
            }
            
            // (CORREÇÃO) Custo EVITADO no DEPÓSITO
            // Ao entregar 'q[t]', economizamos o custo de manter 'q[t]' no depósito.
            obj -= depot.invCost * q[t];    
        }
        // --- FIM DA MODIFICAÇÃO ---

        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize();

        GurobiInventoryResult result;
        result.q.resize(T);
        result.S.resize(T);
        
        int status = model.get(GRB_IntAttr_Status);
        if (status == GRB_OPTIMAL) {
            result.totalCost = model.get(GRB_DoubleAttr_ObjVal);
            if (result.totalCost >= PENALTY_COST * 0.9) {
                 result.totalCost = PENALTY_COST;
            } else {
                for (int t = 0; t < T; ++t) {
                    result.q[t] = q[t].get(GRB_DoubleAttr_X);
                    result.S[t] = S[t].get(GRB_DoubleAttr_X);
                }
            }
        } else {
            result.totalCost = PENALTY_COST; 
        }
        return result;
    }
    catch (GRBException& e) { std::cerr << "Gurobi Error: " << e.getMessage() << std::endl; exit(1); }
    catch (...) { std::cerr << "Unknown error in computePerfectInventory_Gurobi()" << std::endl; exit(1); }
}


ReinsertionData calculate_reinsertion_data(
    const Individual& original_ind, 
    const IRP& irp, 
    const ACO_Params& aco_params,
    int c_id_internal
) 
{
    int c_id_global = c_id_internal + 1;
    ReinsertionData data(irp.nPeriods);
    data.customer_id = c_id_internal;

    Individual original_for_eval = original_ind;
    calculate_total_cost(original_for_eval, irp); 
    data.cost_before_removal = original_for_eval.fitness; 
    data.routes_before_removal = original_ind.routes_per_period;

    Individual temp_sol = original_ind;
    double temp_routing_cost = 0.0;
    for (int t = 0; t < irp.nPeriods; ++t) {
        if (original_ind.deliveries[t][c_id_internal] > 0) {
            temp_sol.deliveries[t][c_id_internal] = 0;
            long period_load = std::accumulate(temp_sol.deliveries[t].begin(), temp_sol.deliveries[t].end(), 0L);
            if (period_load > 0) {
                ACO_Result ares = runACO_for_period(irp, temp_sol.deliveries[t], aco_params, false);
                temp_sol.routes_per_period[t] = (ares.bestCost < PENALTY_COST / 10.0) ? ares.bestRoutes : std::vector<Route>();
            } else {
                temp_sol.routes_per_period[t].clear();
            }
        } 
        for (const auto& route : temp_sol.routes_per_period[t]) {
            temp_routing_cost += route.cost;
        }
    }
    data.routes_after_removal = temp_sol.routes_per_period;
    data.solution_without_customer = temp_sol; 

    double temp_customer_holding_cost = 0.0, temp_depot_holding_cost = 0.0;
    vector<long> customer_inv(irp.nCustomers);
    for (int i = 0; i < irp.nCustomers; ++i) customer_inv[i] = irp.customers[i].initialInv;
    long depot_inv = irp.depots[0].initialInv;
    for (int t = 0; t < irp.nPeriods; ++t) {
        depot_inv += irp.depots[0].production[t];
        long period_delivery_sum = std::accumulate(temp_sol.deliveries[t].begin(), temp_sol.deliveries[t].end(), 0L);
        depot_inv -= period_delivery_sum;
        if (depot_inv > 0) temp_depot_holding_cost += (double)depot_inv * irp.depots[0].invCost;
        for (int c = 0; c < irp.nCustomers; ++c) {
            if (c == c_id_internal) continue; 
            customer_inv[c] += (long)temp_sol.deliveries[t][c];
            customer_inv[c] -= (long)irp.customers[c].demand[t];
            if (customer_inv[c] > irp.customers[c].minLevelInv) {
                temp_customer_holding_cost += (double)(customer_inv[c] - irp.customers[c].minLevelInv) * irp.customers[c].invCost;
            }
        }
    }
    data.cost_after_removal = temp_routing_cost + temp_customer_holding_cost + temp_depot_holding_cost;
    
    long current_inv_c = irp.customers[c_id_internal].initialInv;
    for (int t = 0; t < irp.nPeriods; ++t) {
        data.max_q_inventory[t] = irp.customers[c_id_internal].maxLevelInv - current_inv_c;
        if (data.max_q_inventory[t] < 0) data.max_q_inventory[t] = 0; 
        vector<InsertionOption> options;
        for (const Route& route : temp_sol.routes_per_period[t]) {
            std::pair<int, double> insertion = find_best_insertion(route, c_id_global, irp);
            if (insertion.first != -1) { 
                options.push_back({insertion.second, route.remaining_capacity});
            }
        }
        int routes_used = temp_sol.routes_per_period[t].size();
        bool new_route_available = (routes_used < irp.nVehicles);
        double new_route_cost = irp.costMatrix[0][c_id_global] + irp.costMatrix[c_id_global][0];
        if (new_route_available) {
            options.push_back({new_route_cost, (int)irp.Capacity});
        }
        std::sort(options.begin(), options.end()); 
        GurobiPWLCurve& ft = data.insertion_cost_curves[t];
        int max_kappa_so_far = -1;
        ft.q_points.push_back(0.0);
        ft.cost_points.push_back(0.0);
        for (const auto& opt : options) {
            if (opt.capacity > max_kappa_so_far) {
                ft.q_points.push_back(opt.capacity);
                ft.cost_points.push_back(opt.cost);
                max_kappa_so_far = opt.capacity;
            }
        }
        if (max_kappa_so_far < (int)irp.Capacity) {
             if (max_kappa_so_far >= 0) {
                ft.q_points.push_back(max_kappa_so_far + 1.0);
                ft.cost_points.push_back(PENALTY_COST); 
             }
             ft.q_points.push_back(irp.Capacity);
             ft.cost_points.push_back(PENALTY_COST);
        }
        current_inv_c -= irp.customers[c_id_internal].demand[t];
        if (current_inv_c < irp.customers[c_id_internal].minLevelInv) {
             current_inv_c = irp.customers[c_id_internal].minLevelInv;
        }
    } 
    return data;
}

void busca_local(Individual& original_sol, const IRP& irp, const ACO_Params& aco_params, bool verbose) {
    
 
    // 1. Cria e embaralha a lista de clientes
    std::vector<int> customer_indices(irp.nCustomers);
    std::iota(customer_indices.begin(), customer_indices.end(), 0); // Preenche 0, 1, 2...
    std::shuffle(customer_indices.begin(), customer_indices.end(), rng);

    bool improvement_found = true;
    int loop_count = 0;

    // 2. Loop principal do VND (First Improvement)
    while (improvement_found) {
        improvement_found = false;
        loop_count++;
        
        // 3. Itera sobre a lista embaralhada de clientes
        for (int c_id_internal : customer_indices) {
            
            double cost_before_move = original_sol.fitness;

            // 4. CHAMA A FUNÇÃO DE PRÉ-PROCESSAMENTO
            ReinsertionData data = calculate_reinsertion_data(original_sol, irp, aco_params, c_id_internal);

            // 5. CHAMA O GUROBI
            GurobiInventoryResult result = computePerfectInventory_Gurobi(
                irp.customers[c_id_internal],
                irp.depots[0],
                {}, 
                data.insertion_cost_curves,
                irp.nPeriods,
                (double)irp.Capacity 
            );
            
            if (result.totalCost >= 1e17) {
                // Gurobi falhou em encontrar uma solução
                if(verbose) std::cout << "  Cliente " << (c_id_internal + 1) << ": Gurobi falhou (sem solução). Pulando.\n";
                continue;
            }

            // 6. GERA A NOVA SOLUÇÃO CANDIDATA
            Individual new_sol = data.solution_without_customer;
            for(int t=0; t<irp.nPeriods; ++t) {
                int delivery_q = static_cast<int>(std::round(result.q[t]));
                new_sol.deliveries[t][c_id_internal] = (delivery_q > 0) ? delivery_q : 0;
            }
            
            build_routes_for_individual(new_sol, irp, aco_params);
            
            // 7. AVALIA A NOVA SOLUÇÃO
            check_feasibility(new_sol, irp); 
            
            if (new_sol.is_feasible) {
                calculate_total_cost(new_sol, irp); 
                
                // 8. COMPARA E ATUALIZA (First Improvement)
                if (new_sol.fitness < cost_before_move) {
                    original_sol = new_sol; // Aceita a melhoria
                    improvement_found = true; // Sinaliza para reiniciar o loop
                    
                    if (verbose) {
                        std::cout << "\n>>> MELHORIA ENCONTRADA (Cliente " << (c_id_internal + 1) << ")! "
                                  << "Custo " << cost_before_move << " -> " << new_sol.fitness << " <<<\n";
                    }
                    // 9. REINICIA O LOOP
 
                    break; 
                
                } else {
                    if (verbose) std::cout << "  Cliente " << (c_id_internal + 1) << ": Testado, sem melhoria (Custo: " 
                                           << new_sol.fitness << " vs " << cost_before_move << ").\n";
                }
            } else {
                if (verbose) std::cout << "  Cliente " << (c_id_internal + 1) << ": Gurobi+ACO gerou solução INFACTÍVEL. Descartado.\n";
            }
        } 
    } 

}

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
                  << ", Custo: " << std::fixed << std::setprecision(2) << route.cost << ")\n";
    }
}

