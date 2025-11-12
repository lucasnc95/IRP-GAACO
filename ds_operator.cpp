#include "ds_operator.hpp"
//#include "aco.hpp"      // Para runACO_for_period
#include "utils.hpp"    // Para randint
#include "route.hpp"
#include "route_builder.hpp"
#include "evaluation.hpp" // Para calculate_total_cost, build_routes_for_individual
#include "ga.hpp"           // Para 'build_routes_for_individual'
#include <algorithm>      // Para std::sort, std::find
#include <map>
#include <set>
#include <iostream>       // Para std::cout
#include <iomanip>        // Para std::fixed, std::setprecision

using std::vector;


const double PENALTY_COST = 1e7;

// --- ESTRUTURAS AUXILIARES ---
struct InsertionOption {
    double cost;  
    int capacity; 
    
    // Ordena por custo crescente. Se empate, menor capacidade primeiro.
    bool operator<(const InsertionOption& other) const {
        if (std::abs(cost - other.cost) > 1e-9) return cost < other.cost;
        return capacity < other.capacity;
    }
};

// --- FUNÇÃO HELPER: Encontrar melhor inserção em uma rota ---
std::pair<int, double> find_best_insertion(const Route& route, int customer_id, const IRP& irp) {
    if (route.visits.empty() || route.visits.size() < 2) return {-1, PENALTY_COST}; 
    
    double min_delta_cost = PENALTY_COST;
    int best_pos = -1;
    
    // Itera sobre as arestas da rota (0->A, A->B, ..., Z->0)
    for (size_t i = 0; i < route.visits.size() - 1; ++i) {
        int node_j = route.visits[i];
        int node_k = route.visits[i+1];
        
        // Custo de inserir 'customer_id' entre j e k
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
    int T, double qmax) 
{
    try {
        GRBEnv env(true);
        env.set("LogToConsole", "0"); 
        env.set(GRB_IntParam_Threads, 0); 
        env.set(GRB_IntParam_NumericFocus, 3); 
        env.start();
        GRBModel model(env);

        // Foco em zerar o gap (agora é um MIP real, então isso é importante)
        model.set(GRB_DoubleParam_MIPGap, 0.0);
        
        std::vector<GRBVar> q(T), S(T);
        for (int t = 0; t < T; ++t) {
            
            // <-- MUDANÇA PRINCIPAL: GRB_INTEGER -->
            // Define a quantidade entregue como Inteira.
            q[t] = model.addVar(0.0, qmax, 0.0, GRB_INTEGER, "q_" + std::to_string(t+1));
            
            // O estoque S[t] pode continuar contínuo, pois é consequência
            // linear de variáveis inteiras (estoque inicial + entregas - demandas).
            // Mas pode ser GRB_CONTINUOUS sem problemas.
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

        GRBLinExpr obj = 0.0;
        for (int t = 0; t < T; ++t) {
            obj += cust.invCost * S[t];     
            obj -= depot.invCost * q[t]; // Custo evitado

            if (!reinsertion_curves.empty() && !reinsertion_curves[t].q_points.empty()) {
                GRBVar c_q = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "c_q_" + std::to_string(t+1));
                const auto& curve = reinsertion_curves[t];
                
                model.addGenConstrPWL(q[t], c_q,
                    (int)curve.q_points.size(),
                    curve.q_points.data(), curve.cost_points.data(),
                    "pwl_reinsert_t" + std::to_string(t+1));
                
                obj += c_q; 
            } 
        }

        model.setObjective(obj, GRB_MINIMIZE);
        model.optimize();

        GurobiInventoryResult result;
        result.q.resize(T);
        result.S.resize(T);
        
        int status = model.get(GRB_IntAttr_Status);
        if (status == GRB_OPTIMAL) {
            result.totalCost = model.get(GRB_DoubleAttr_ObjVal);
            // Verifica factibilidade baseada no custo de penalidade
            if (result.totalCost >= PENALTY_COST * 0.9) {
                 result.totalCost = PENALTY_COST;
            } else {
                for (int t = 0; t < T; ++t) {
                    // Agora q[t] já é inteiro, mas get() retorna double
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

    // 1. Salvar Dados "Antes"
    // Precisamos de uma cópia para calcular o fitness sem alterar a const original,
    // caso ela não tenha sido avaliada recentemente.
    Individual original_for_eval = original_ind;
    calculate_total_cost(original_for_eval, irp); 
    data.cost_before_removal = original_for_eval.fitness; 
    data.routes_before_removal = original_ind.routes_per_period;

    // 2. Criar Solução Temporária e Reroteirizar (USANDO SWEEP)
    Individual temp_sol = original_ind;
    double temp_routing_cost = 0.0;

    for (int t = 0; t < irp.nPeriods; ++t) {
        // Se o cliente tinha entrega, remove e refaz as rotas
        if (original_ind.deliveries[t][c_id_internal] > 0) {
            temp_sol.deliveries[t][c_id_internal] = 0; // Zera a entrega
            
            // Chama o SWEEP para reconstruir as rotas dos clientes restantes
            temp_sol.routes_per_period[t] = build_routes_with_sweep(
                irp, 
                temp_sol.deliveries[t], 
                aco_params // Passa params para decidir se usa busca local interna
            );
        } 
        // Se o cliente não estava no dia, as rotas originais (copiadas) são mantidas
        
        // Acumula custo das rotas (novas ou mantidas)
        for (const auto& route : temp_sol.routes_per_period[t]) {
            temp_routing_cost += route.cost;
        }
    }
    
    data.routes_after_removal = temp_sol.routes_per_period;
    data.solution_without_customer = temp_sol; 

    // 3. Calcular Custo Parcial "Depois" (Simulação de Inventário)
    double temp_customer_holding_cost = 0.0;
    double temp_depot_holding_cost = 0.0;
    
    vector<long> customer_inv(irp.nCustomers);
    for (int i = 0; i < irp.nCustomers; ++i) customer_inv[i] = irp.customers[i].initialInv;
    long depot_inv = irp.depots[0].initialInv;

    for (int t = 0; t < irp.nPeriods; ++t) {
        // Depósito
        depot_inv += irp.depots[0].production[t];
        long period_delivery_sum = std::accumulate(temp_sol.deliveries[t].begin(), temp_sol.deliveries[t].end(), 0L);
        depot_inv -= period_delivery_sum;
        if (depot_inv > 0) temp_depot_holding_cost += (double)depot_inv * irp.depots[0].invCost;

        // Clientes (exceto o removido)
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
    
    // 4. Construção das Curvas de Custo (F_t) e Max Inventário
    long current_inv_c = irp.customers[c_id_internal].initialInv;
    
    for (int t = 0; t < irp.nPeriods; ++t) {
        
        // 4a. Max Qtd permitida por Inventário
        data.max_q_inventory[t] = irp.customers[c_id_internal].maxLevelInv - current_inv_c;
        if (data.max_q_inventory[t] < 0) data.max_q_inventory[t] = 0; 

        // 4b. Coleta opções de inserção nas rotas do SWEEP
        vector<InsertionOption> options;
        
        // Inserção em rotas existentes
        for (const Route& route : temp_sol.routes_per_period[t]) {
            if (route.remaining_capacity > 0) {
                std::pair<int, double> insertion = find_best_insertion(route, c_id_global, irp);
                if (insertion.first != -1) { 
                    options.push_back({insertion.second, route.remaining_capacity});
                }
            }
        }
        
        // Nova rota (se houver veículo)
        int routes_used = temp_sol.routes_per_period[t].size();
        if (routes_used < irp.nVehicles) {
            double new_route_cost = irp.costMatrix[0][c_id_global] + irp.costMatrix[c_id_global][0];
            options.push_back({new_route_cost, (int)irp.Capacity});
        }

        // Ordena por menor custo
        std::sort(options.begin(), options.end()); 

        // 4c. Constrói a Curva PWL como FUNÇÃO ESCADA
        GurobiPWLCurve& ft = data.insertion_cost_curves[t];
        
        // Ponto (0, 0)
        ft.q_points.push_back(0.0);
        ft.cost_points.push_back(0.0);
        
        double epsilon = 0.1; 
        int max_cap_covered = 0;
        
        for (const auto& opt : options) {
            // Só adiciona se a opção estender a capacidade coberta
            if (opt.capacity > max_cap_covered) {
                
                // Degrau sobe logo após a capacidade anterior
                double start_q = std::max(epsilon, (double)max_cap_covered + epsilon);
                
                if (start_q <= opt.capacity) {
                    // Início do degrau
                    ft.q_points.push_back(start_q);
                    ft.cost_points.push_back(opt.cost);
                    
                    // Fim do degrau (custo constante)
                    ft.q_points.push_back((double)opt.capacity);
                    ft.cost_points.push_back(opt.cost);
                }
                
                max_cap_covered = opt.capacity;
            }
        }
        
        // Fecha com penalidade
        if (max_cap_covered < (int)irp.Capacity) {
             ft.q_points.push_back((double)max_cap_covered + epsilon);
             ft.cost_points.push_back(PENALTY_COST);
             
             ft.q_points.push_back((double)irp.Capacity);
             ft.cost_points.push_back(PENALTY_COST);
        }

        // 4d. Simula inventário para o próximo dia
        current_inv_c -= irp.customers[c_id_internal].demand[t];
        if (current_inv_c < irp.customers[c_id_internal].minLevelInv) {
             current_inv_c = irp.customers[c_id_internal].minLevelInv;
        }
    } 

    return data;
}

void run_vnd_ds_operator(Individual& original_sol, const IRP& irp, const ACO_Params& aco_params, bool verbose) {
    
 
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



double incremental_remove_customer(std::vector<Route>& routes, int c_id_global, const IRP& irp, int qty_to_remove) {
    for (auto it = routes.begin(); it != routes.end(); ++it) {
        Route& route = *it;
        
        for (size_t i = 1; i < route.visits.size() - 1; ++i) {
            if (route.visits[i] == c_id_global) {
                // Encontrou o cliente na posição 'i'
                int prev = route.visits[i-1];
                int curr = route.visits[i]; // c_id_global
                int next = route.visits[i+1];

                // Delta = (Custo antigo das arestas removidas) - (Custo da nova aresta)
                double cost_removed = irp.costMatrix[prev][curr] + irp.costMatrix[curr][next];
                double cost_added = irp.costMatrix[prev][next];
                double cost_reduction = cost_removed - cost_added;

                // Atualiza a rota
                route.cost -= cost_reduction;
                route.remaining_capacity += qty_to_remove;
                route.visits.erase(route.visits.begin() + i);

                // Se a rota ficou vazia (0->0), remove ela do vetor
                if (route.visits.size() <= 2) {
                    double empty_route_cost = route.cost; // Deve ser ~0 ou dist[0][0]
                    routes.erase(it);
                    return cost_reduction + empty_route_cost; // A redução total inclui o custo residual da rota vazia
                }

                return cost_reduction;
            }
        }
    }
    return 0.0; // Cliente não estava nas rotas deste dia
}

/**
 * @brief Insere cliente na melhor rota e atualiza custo/capacidade incrementalmente.
 * @return O valor do AUMENTO de custo.
 */
double incremental_insert_customer(std::vector<Route>& routes, int c_id_global, const IRP& irp, int qty_to_insert) {
    int best_r_idx = -1;
    int best_pos = -1;
    double min_cost_increase = 1e18;

    // 1. Tenta inserir em rotas existentes
    for (size_t r = 0; r < routes.size(); ++r) {
        if (routes[r].remaining_capacity >= qty_to_insert) {
            // Usa a lógica de Cheapest Insertion Delta
            for (size_t i = 0; i < routes[r].visits.size() - 1; ++i) {
                int u = routes[r].visits[i];
                int v = routes[r].visits[i+1];
                double delta = irp.costMatrix[u][c_id_global] + 
                               irp.costMatrix[c_id_global][v] - 
                               irp.costMatrix[u][v];
                
                if (delta < min_cost_increase) {
                    min_cost_increase = delta;
                    best_r_idx = r;
                    best_pos = i + 1; // Insere depois de u
                }
            }
        }
    }

    // 2. Tenta criar nova rota (se houver veículo)
    if (routes.size() < irp.nVehicles) {
        double new_route_cost = irp.costMatrix[0][c_id_global] + irp.costMatrix[c_id_global][0];
        if (new_route_cost < min_cost_increase) {
            min_cost_increase = new_route_cost;
            best_r_idx = -99; // Código para nova rota
        }
    }

    // 3. Aplica a inserção
    if (best_r_idx == -1 && min_cost_increase >= 1e17) {
        return 1e18; // Falha na inserção (sem capacidade)
    }

    if (best_r_idx == -99) {
        // Cria nova rota
        Route new_route;
        new_route.visits = {0, c_id_global, 0};
        new_route.cost = min_cost_increase;
        new_route.remaining_capacity = irp.Capacity - qty_to_insert;
        routes.push_back(new_route);
    } else {
        // Insere na rota existente
        Route& r = routes[best_r_idx];
        r.visits.insert(r.visits.begin() + best_pos, c_id_global);
        r.cost += min_cost_increase;
        r.remaining_capacity -= qty_to_insert;
    }
    
    return min_cost_increase;
}

// --- FIM DAS FUNÇÕES INCREMENTAIS ---


// --- FUNÇÃO AUXILIAR: CÁLCULO RÁPIDO DE CUSTO DE ESTOQUE DE UM CLIENTE ---
double calculate_single_customer_inventory_cost(const IRP& irp, int c_idx, const std::vector<int>& deliveries) {
    double cost = 0.0;
    long inv = irp.customers[c_idx].initialInv;
    for(int t=0; t<irp.nPeriods; ++t) {
        inv += deliveries[t];
        inv -= irp.customers[c_idx].demand[t];
        if(inv > irp.customers[c_idx].minLevelInv) {
            cost += (double)(inv - irp.customers[c_idx].minLevelInv) * irp.customers[c_idx].invCost;
        }
    }
    return cost;
}

// --- FUNÇÃO AUXILIAR: CÁLCULO RÁPIDO DE CUSTO DE ESTOQUE DO DEPÓSITO ---
double calculate_depot_inventory_cost(const IRP& irp, const std::vector<std::vector<int>>& deliveries_matrix) {
    double cost = 0.0;
    long inv = irp.depots[0].initialInv;
    for(int t=0; t<irp.nPeriods; ++t) {
        inv += irp.depots[0].production[t];
        long period_del = 0;
        for(int c=0; c<irp.nCustomers; ++c) period_del += deliveries_matrix[t][c];
        inv -= period_del;
        if(inv > irp.depots[0].minLevelInv) {
            cost += (double)inv * irp.depots[0].invCost;
        }
    }
    return cost;
}


// ============================================================================
// === BUSCA LOCAL: SIMPLE REINSERTION (OTIMIZADA) ===
// ============================================================================

void run_simple_reinsertion_search(
    Individual& sol, 
    const IRP& irp, 
    const ACO_Params& aco_params,
    bool verbose
) {
    if (verbose) std::cout << "\n=== BUSCA LOCAL: SIMPLE REINSERTION (Incremental) ===\n";

    std::vector<int> customer_indices(irp.nCustomers);
    std::iota(customer_indices.begin(), customer_indices.end(), 0);
    std::shuffle(customer_indices.begin(), customer_indices.end(), rng);

    bool improvement = true;
    while (improvement) {
        improvement = false;
        
        for (int c_id_internal : customer_indices) {
            int c_id_global = c_id_internal + 1;
            double current_fitness = sol.fitness;

            // --- A: REMOÇÃO E CÁLCULO DE DELTA DE CUSTO INICIAL ---
            Individual temp_sol = sol;
            double routing_delta = 0.0;
            
            std::vector<int> old_deliveries_c(irp.nPeriods);
            for (int t = 0; t < irp.nPeriods; ++t) {
                int qty = temp_sol.deliveries[t][c_id_internal];
                if (qty > 0) {
                    old_deliveries_c[t] = qty;
                    routing_delta -= incremental_remove_customer(temp_sol.routes_per_period[t], c_id_global, irp, qty);
                    temp_sol.deliveries[t][c_id_internal] = 0;
                }
            }
            
            // <-- CORREÇÃO DO ERRO DE COMPILAÇÃO AQUI -->
            // Em vez de sol.deliveries_per_customer(), extraímos manualmente
            double old_cust_inv_cost = calculate_single_customer_inventory_cost(irp, c_id_internal, old_deliveries_c);

            // --- B: PREPARAR DADOS PARA GUROBI ---
            std::vector<GurobiPWLCurve> curves(irp.nPeriods);
            long current_inv = irp.customers[c_id_internal].initialInv;
            
            for (int t = 0; t < irp.nPeriods; ++t) {
                // 1. Custo de Inserção
                vector<InsertionOption> options;
                for (const Route& route : temp_sol.routes_per_period[t]) {
                    std::pair<int, double> ins = find_best_insertion(route, c_id_global, irp);
                    if (ins.first != -1) {
                        options.push_back({ins.second, route.remaining_capacity});
                    }
                }
                if (temp_sol.routes_per_period[t].size() < irp.nVehicles) {
                    double direct_cost = irp.costMatrix[0][c_id_global] + irp.costMatrix[c_id_global][0];
                    options.push_back({direct_cost, (int)irp.Capacity});
                }
                std::sort(options.begin(), options.end());

                // Monta curva PWL
                GurobiPWLCurve& ft = curves[t];
                ft.q_points.push_back(0.0); ft.cost_points.push_back(0.0);
                int max_kappa = -1;
                for(auto& opt : options) {
                    if(opt.capacity > max_kappa) {
                        ft.q_points.push_back(opt.capacity);
                        ft.cost_points.push_back(opt.cost);
                        max_kappa = opt.capacity;
                    }
                }
                if (max_kappa < (int)irp.Capacity) {
                    if (max_kappa >= 0) {
                        ft.q_points.push_back(max_kappa + 1.0); ft.cost_points.push_back(PENALTY_COST);
                    }
                    ft.q_points.push_back(irp.Capacity); ft.cost_points.push_back(PENALTY_COST);
                }

                // 2. Simula inventário para consumo
                current_inv -= irp.customers[c_id_internal].demand[t];
                if (current_inv < irp.customers[c_id_internal].minLevelInv) 
                    current_inv = irp.customers[c_id_internal].minLevelInv;
            }

            // --- C: GUROBI ---
            GurobiInventoryResult res = computePerfectInventory_Gurobi(
                irp.customers[c_id_internal], irp.depots[0], {}, curves, irp.nPeriods, (double)irp.Capacity
            );

            if (res.totalCost >= PENALTY_COST * 0.9) continue; 

            // --- D: REINSERÇÃO INCREMENTAL ---
            bool possible = true;
            std::vector<int> new_deliveries_c(irp.nPeriods, 0);
            
            for (int t = 0; t < irp.nPeriods; ++t) {
                int q = static_cast<int>(std::round(res.q[t]));
                if (q > 0) {
                    new_deliveries_c[t] = q;
                    double increase = incremental_insert_customer(temp_sol.routes_per_period[t], c_id_global, irp, q);
                    if (increase >= 1e17) {
                        possible = false;
                        break;
                    }
                    routing_delta += increase;
                }
            }

            if (possible) {
                // --- E: ATUALIZAÇÃO DOS CUSTOS DA SOLUÇÃO ---
                temp_sol.routing_cost = sol.routing_cost + routing_delta;
                
                // Extrai entregas novas para cálculo de custo
                // (new_deliveries_c já está preenchido)
                
                double new_cust_inv_cost = calculate_single_customer_inventory_cost(irp, c_id_internal, new_deliveries_c);
                double cust_inv_delta = new_cust_inv_cost - old_cust_inv_cost;
                temp_sol.customer_holding_cost = sol.customer_holding_cost + cust_inv_delta;

                // Atualiza plano de entregas completo para recalcular depósito
                for(int t=0; t<irp.nPeriods; ++t) temp_sol.deliveries[t][c_id_internal] = new_deliveries_c[t];
                
                // Depósito (recalculado)
                temp_sol.depot_holding_cost = calculate_depot_inventory_cost(irp, temp_sol.deliveries);

                temp_sol.fitness = temp_sol.routing_cost + temp_sol.customer_holding_cost + temp_sol.depot_holding_cost;
                temp_sol.is_feasible = true; 

                // --- F: VERIFICAÇÃO DE MELHORIA ---
                if (temp_sol.fitness < sol.fitness - 1e-9) {
                    if (verbose) {
                        std::cout << "  [SimpleReinsertion] Cliente " << c_id_global 
                                  << ": Melhoria " << std::fixed << std::setprecision(2) << sol.fitness 
                                  << " -> " << temp_sol.fitness << "\n";
                    }
                    sol = temp_sol;
                    improvement = true;
                    std::shuffle(customer_indices.begin(), customer_indices.end(), rng);
                    break; 
                }
            }
        } 
    } 
}