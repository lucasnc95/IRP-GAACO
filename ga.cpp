#include "ds_operator.hpp"
#include "route_builder.hpp"
#include "ga.hpp"
#include "utils.hpp"
#include "evaluation.hpp" 
//#include "aco.hpp"
#include "route.hpp"
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <utility>
#include <set>

using std::vector;

const double PENALTY_COST = 1e7;

static bool check_aco_consistency(const vector<int>& deliveries, const vector<Route>& routes) {
    std::set<int> customers_requiring_delivery;
    for (int c = 0; c < deliveries.size(); ++c) {
        if (deliveries[c] > 0) {
            customers_requiring_delivery.insert(c + 1); // 1-based
        }
    }
    std::set<int> customers_in_routes;
    for (const Route& route : routes) {
        for (int node : route.visits) {
            if (node > 0) customers_in_routes.insert(node);
        }
    }
    return customers_requiring_delivery == customers_in_routes;
}




std::pair<Individual, Individual> one_point_crossover_customer(const Individual& a, const Individual& b, const IRP& irp) {
    Individual child1(irp.nPeriods, irp.nCustomers);
    Individual child2(irp.nPeriods, irp.nCustomers);
    int crosspoint = randint(0, irp.nCustomers);
    for (int c = 0; c < irp.nCustomers; ++c) {
        if (c < crosspoint) {
            for (int t = 0; t < irp.nPeriods; ++t) {
                child1.deliveries[t][c] = a.deliveries[t][c];
                child2.deliveries[t][c] = b.deliveries[t][c];
            }
        } else {
            for (int t = 0; t < irp.nPeriods; ++t) {
                child1.deliveries[t][c] = b.deliveries[t][c];
                child2.deliveries[t][c] = a.deliveries[t][c];
            }
        }
    }
    return {child1, child2};
}
std::pair<Individual, Individual> two_point_crossover_customer(const Individual& a, const Individual& b, const IRP& irp) {
    Individual child1(irp.nPeriods, irp.nCustomers);
    Individual child2(irp.nPeriods, irp.nCustomers);
    int p1 = randint(0, irp.nCustomers);
    int p2 = randint(0, irp.nCustomers);
    if (p1 > p2) std::swap(p1, p2);
    for (int c = 0; c < irp.nCustomers; ++c) {
        if (c >= p1 && c < p2) {
            for (int t = 0; t < irp.nPeriods; ++t) {
                child1.deliveries[t][c] = b.deliveries[t][c];
                child2.deliveries[t][c] = a.deliveries[t][c];
            }
        } else {
            for (int t = 0; t < irp.nPeriods; ++t) {
                child1.deliveries[t][c] = a.deliveries[t][c];
                child2.deliveries[t][c] = b.deliveries[t][c];
            }
        }
    }
    return {child1, child2};
}

Individual tournamentSelect(const vector<Individual>& pop, int k) {
    int n = pop.size();
    int best_idx = randint(0, n - 1);
    for (int i = 1; i < k; ++i) {
        int cand_idx = randint(0, n - 1);
        if (pop[cand_idx].fitness < pop[best_idx].fitness) {
            best_idx = cand_idx;
        }
    }
    return pop[best_idx];
}


// Individual make_new_heuristic_individual(const IRP& irp) {
//     Individual ind(irp.nPeriods, irp.nCustomers);
//     vector<long> current_inv(irp.nCustomers);
//     for(int i = 0; i < irp.nCustomers; ++i) current_inv[i] = irp.customers[i].initialInv;

//     long fleet_capacity = (long)irp.nVehicles * irp.Capacity;

//     for (int t = 0; t < irp.nPeriods; ++t) {
//         vector<int> min_deliveries_today(irp.nCustomers, 0);
//         long period_load_min = 0;

//         // --- PASSO 1: CALCULAR ENTREGAS MÍNIMAS OBRIGATÓRIAS ---
//         for (int c = 0; c < irp.nCustomers; ++c) {
//             long inv_after_demand = current_inv[c] - irp.customers[c].demand[t];
            
//             if (inv_after_demand < irp.customers[c].minLevelInv) {
//                 long needed = irp.customers[c].minLevelInv - inv_after_demand;
//                 long space_available = irp.customers[c].maxLevelInv - current_inv[c];
//                 long delivery_amount = std::min({needed, space_available, (long)irp.Capacity});
                
//                 if (delivery_amount > 0) {
//                     min_deliveries_today[c] = delivery_amount;
//                     period_load_min += delivery_amount;
//                 }
//             }
//         }
        
//         ind.deliveries[t] = min_deliveries_today;
//         long current_period_load = period_load_min;
        

//         if (current_period_load > fleet_capacity) {
//             for (int c = 0; c < irp.nCustomers; ++c) {
//                 current_inv[c] += (long)ind.deliveries[t][c] - irp.customers[c].demand[t];
//             }
//             continue; // Pula para o próximo dia
//         }

//         // --- PASSO 2: ENTREGAS OPORTUNISTAS (30% de chance) ---
//         long remaining_fleet_cap = fleet_capacity - current_period_load;
        
//         vector<int> customer_order(irp.nCustomers);
//         std::iota(customer_order.begin(), customer_order.end(), 0);
//         std::shuffle(customer_order.begin(), customer_order.end(), rng);

//         for (int c : customer_order) {
//             if (remaining_fleet_cap <= 0) break;
            
//             if (randreal() < 0.30) { 
//                 int N = randint(1, irp.nPeriods-t);
//                 long q_extra = 0;
//                 for (int t_future = t + 1; t_future <= t + N && t_future < irp.nPeriods; ++t_future) {
//                     q_extra += irp.customers[c].demand[t_future];
//                 }
//                 if (q_extra == 0) continue;
                
//                 long current_delivery = ind.deliveries[t][c];
//                 long space = irp.customers[c].maxLevelInv - (current_inv[c] + current_delivery);
//                 long max_q_vehicle = (long)irp.Capacity - current_delivery;
                
//                 q_extra = std::min({q_extra, space, max_q_vehicle, remaining_fleet_cap});

//                 if (q_extra > 0) {
//                     ind.deliveries[t][c] += q_extra;
//                     remaining_fleet_cap -= q_extra;
//                 }
//             }
//         } // Fim do Passo 2
        
//         // --- PASSO 3: ATUALIZA O INVENTÁRIO (para o próximo período t+1) ---
//         for (int c = 0; c < irp.nCustomers; ++c) {
//             current_inv[c] += (long)ind.deliveries[t][c] - irp.customers[c].demand[t];
//             // Garante que o inventário não fique "negativo" para a simulação do dia seguinte
//             if(current_inv[c] < irp.customers[c].minLevelInv) {
//                 current_inv[c] = irp.customers[c].minLevelInv;
//             }
//         }
//     } 
    
//     return ind;
// }



Individual make_new_heuristic_individual(const IRP& irp) {
    Individual ind(irp.nPeriods, irp.nCustomers);
    vector<long> current_inv(irp.nCustomers);
    for(int i = 0; i < irp.nCustomers; ++i) current_inv[i] = irp.customers[i].initialInv;

    long fleet_capacity = (long)irp.nVehicles * irp.Capacity;

    for (int t = 0; t < irp.nPeriods; ++t) {
        long period_load = 0;
        vector<int> mandatory_custs, optional_custs;

        for (int c = 0; c < irp.nCustomers; ++c) {
            if (current_inv[c] < irp.customers[c].demand[t]) {
                mandatory_custs.push_back(c);
            } else {
                optional_custs.push_back(c);
            }
        }

        // Atende clientes obrigatórios
        for (int c : mandatory_custs) {
            long needed = irp.customers[c].demand[t] - current_inv[c];
            long space = irp.customers[c].maxLevelInv - (current_inv[c] + needed);
            int extra = (space > 0) ? randint(0, std::min((long)20, space)) : 0; // Adiciona um extra aleatório pequeno
            int q = std::min((long)irp.Capacity, needed + extra);
            
            if (period_load + q <= fleet_capacity) {
                ind.deliveries[t][c] = q;
                period_load += q;
            }
        }
        
        // Atende clientes opcionais aleatoriamente
        std::shuffle(optional_custs.begin(), optional_custs.end(), rng);
        for (int c : optional_custs) {
            if (period_load >= fleet_capacity) break;
            if (randreal() < 0.3) { // Chance de atender um cliente opcional
                long space = irp.customers[c].maxLevelInv - current_inv[c];
                long max_q = std::min({space, (long)irp.Capacity, fleet_capacity - period_load});
                if (max_q > 0) {
                    int q = randint(1, max_q);
                    ind.deliveries[t][c] = q;
                    period_load += q;
                }
            }
        }

        // Atualiza o inventário para o próximo período
        for (int c = 0; c < irp.nCustomers; ++c) {
            current_inv[c] += (long)ind.deliveries[t][c] - irp.customers[c].demand[t];
        }
    }
    return ind;
}


// void build_routes_for_individual(Individual& ind, const IRP& irp, const ACO_Params& aco_params) {
//     ind.routes_per_period.assign(irp.nPeriods, std::vector<Route>());
//     const double LARGE_COST = 1e18;

//     for (int t = 0; t < irp.nPeriods; ++t) {
//         long period_delivery_sum = std::accumulate(ind.deliveries[t].begin(), ind.deliveries[t].end(), 0L);
        
//         if (period_delivery_sum > 0) {
//             ACO_Result ares = runACO_for_period(irp, ind.deliveries[t], aco_params, false);
            
//             // Armazena as rotas. Se o ACO falhou (custo alto), 
//             // as rotas estarão vazias ou inválidas, e 'check_feasibility' vai pegar.
//             if(ares.bestCost < LARGE_COST / 10.0) {
//                 ind.routes_per_period[t] = ares.bestRoutes;
//             } else {
//                 ind.routes_per_period[t].clear(); // Garante que está vazio se o ACO falhou
//             }
//         }
//     }
// }




void build_routes_for_individual(Individual& ind, const IRP& irp, const ACO_Params& aco_params) {
    ind.routes_per_period.assign(irp.nPeriods, std::vector<Route>());

    for (int t = 0; t < irp.nPeriods; ++t) {
        long period_delivery_sum = std::accumulate(ind.deliveries[t].begin(), ind.deliveries[t].end(), 0L);
        
        if (period_delivery_sum > 0) {
            // <-- MUDANÇA: Chama o SWEEP em vez do ACO -->
            std::vector<Route> routes = build_routes_with_sweep(
                irp, 
                ind.deliveries[t], 
                aco_params // Passa os parâmetros (para a flag pLocalSearch)
            );
            
            // Armazena as rotas. Se o Sweep falhou (ex: > K veículos),
            // 'routes' estará vazio e 'check_feasibility' irá falhar.
            ind.routes_per_period[t] = routes;
        }
        // Se a soma for 0, o vetor de rotas já está vazio (correto)
    }
}













Individual make_new_heuristic_individual(const IRP& irp, const ACO_Params& aco_params) {
    
    Individual ind(irp.nPeriods, irp.nCustomers);
    
    // --- PASSO 1: Sortear 20% dos clientes para reservar ---
    int num_reserved = static_cast<int>(irp.nCustomers * 0.7);
    if (num_reserved == 0 && irp.nCustomers > 0) num_reserved = 1; // Garante pelo menos 1
    
    vector<int> customer_indices(irp.nCustomers);
    std::iota(customer_indices.begin(), customer_indices.end(), 0); // 0, 1, 2...
    std::shuffle(customer_indices.begin(), customer_indices.end(), rng);
    
    std::set<int> reserved_customers; // Clientes a serem inseridos depois (0-based)
    for(int i=0; i < num_reserved; ++i) {
        reserved_customers.insert(customer_indices[i]);
    }
    
    // --- PASSO 2: Construir solução para os 80% (Base) ---
    vector<long> current_inv(irp.nCustomers);
    for(int i = 0; i < irp.nCustomers; ++i) current_inv[i] = irp.customers[i].initialInv;
    long fleet_capacity = (long)irp.nVehicles * irp.Capacity;

    for (int t = 0; t < irp.nPeriods; ++t) {
        vector<int> min_deliveries_today(irp.nCustomers, 0);
        long period_load_min = 0;

        for (int c = 0; c < irp.nCustomers; ++c) {
            // PULA CLIENTES RESERVADOS
            if (reserved_customers.count(c)) continue; 

            long inv_after_demand = current_inv[c] - irp.customers[c].demand[t];
            if (inv_after_demand < irp.customers[c].minLevelInv) {
                long needed = irp.customers[c].minLevelInv - inv_after_demand;
                long space_available = irp.customers[c].maxLevelInv - current_inv[c];
                long delivery_amount = std::min({needed, space_available, (long)irp.Capacity});
                
                if (delivery_amount > 0) {
                    min_deliveries_today[c] = delivery_amount;
                    period_load_min += delivery_amount;
                }
            }
        }
        
        ind.deliveries[t] = min_deliveries_today;
        long current_period_load = period_load_min;
        
        if (current_period_load > fleet_capacity) {
            // (Simula para o próximo dia)
            for (int c = 0; c < irp.nCustomers; ++c) {
                if (reserved_customers.count(c)) continue; // Pula
                current_inv[c] += (long)ind.deliveries[t][c] - irp.customers[c].demand[t];
                if(current_inv[c] < irp.customers[c].minLevelInv) current_inv[c] = irp.customers[c].minLevelInv;
            }
            continue; 
        }

        // Entregas oportunistas (apenas para os 80%)
        long remaining_fleet_cap = fleet_capacity - current_period_load;
        vector<int> customer_order;
        for(int c=0; c < irp.nCustomers; ++c) if(!reserved_customers.count(c)) customer_order.push_back(c);
        std::shuffle(customer_order.begin(), customer_order.end(), rng);

        for (int c : customer_order) {
            if (remaining_fleet_cap <= 0) break;
            if (randreal() < 0.30) { 
                int N = randint(1, 3);
                long q_extra = 0;
                for (int t_future = t + 1; t_future <= t + N && t_future < irp.nPeriods; ++t_future) {
                    q_extra += irp.customers[c].demand[t_future];
                }
                if (q_extra == 0) continue;
                
                long current_delivery = ind.deliveries[t][c];
                long space = irp.customers[c].maxLevelInv - (current_inv[c] + current_delivery);
                long max_q_vehicle = (long)irp.Capacity - current_delivery;
                q_extra = std::min({q_extra, space, max_q_vehicle, remaining_fleet_cap});

                if (q_extra > 0) {
                    ind.deliveries[t][c] += q_extra;
                    remaining_fleet_cap -= q_extra;
                }
            }
        } 
        
        // Simula inventário para o próximo dia (apenas 80%)
        for (int c = 0; c < irp.nCustomers; ++c) {
            if (reserved_customers.count(c)) continue; // Pula
            current_inv[c] += (long)ind.deliveries[t][c] - irp.customers[c].demand[t];
            if(current_inv[c] < irp.customers[c].minLevelInv) {
                current_inv[c] = irp.customers[c].minLevelInv;
            }
        }
    } // Fim do loop 't' (plano 80% está feito)
    
    // --- PASSO 3: Criar as rotas base (sem os clientes reservados) ---
    build_routes_for_individual(ind, irp, aco_params);
    // (Neste ponto, 'ind' é uma solução parcial factível para 80% dos clientes)

    // --- PASSO 4: Chamar a "busca local" (lógica Gurobi) para INSERIR os clientes ---
    
 
              
    for (int c_id_internal : reserved_customers) {
        
        // 1. Calcula dados de reinserção (custos F_t(q)) para 'c'
        //    com base nas rotas ATUAIS (que ainda não o contêm)
        ReinsertionData data = calculate_reinsertion_data(ind, irp, aco_params, c_id_internal);
        
        // 2. Chama Gurobi para encontrar o plano de inventário ótimo para 'c'
        GurobiInventoryResult result = computePerfectInventory_Gurobi(
            irp.customers[c_id_internal],
            irp.depots[0],
            {}, 
            data.insertion_cost_curves,
            irp.nPeriods,
            (double)irp.Capacity 
        );

        if (result.totalCost >= PENALTY_COST * 0.9) {
            // Gurobi falhou em encontrar um plano de inserção.
            // O indivíduo provavelmente será infactível.
             std::cerr << "AVISO: Falha do Gurobi ao inserir cliente " << (c_id_internal+1) << "\n";
            continue; // Tenta o próximo cliente
        }

        // 3. ATUALIZA O GENÓTIPO: Adiciona o plano do Gurobi ao 'ind.deliveries'
        for(int t=0; t < irp.nPeriods; ++t) {
            int delivery_q = static_cast<int>(std::round(result.q[t]));
            if (delivery_q > 0) {
                ind.deliveries[t][c_id_internal] = delivery_q;
            }
        }
        
        // 4. ATUALIZA O FENÓTIPO: Reconstrói as rotas
        // (Agora incluindo o cliente 'c')
        build_routes_for_individual(ind, irp, aco_params);
        
        // (Verificação de factibilidade não é necessária aqui, 
        // pois a próxima inserção (passo 1) recalculará tudo)
    }

    // O indivíduo 'ind' agora está completo (100% dos clientes)
    return ind;
}








Individual run_genetic_algorithm(const IRP& irp, const GA_Params& ga_params, const ACO_Params& aco_params, bool verbose) {
    
    const double LARGE_COST = 1e18;
    vector<Individual> pop;
    pop.reserve(ga_params.popSize);
    
    std::cout << "Inicializando população..." << std::endl;
    for (int i = 0; i < ga_params.popSize; ++i) {
        Individual ind = make_new_heuristic_individual(irp, aco_params);
        
        build_routes_for_individual(ind, irp, aco_params);
       // busca_local(ind, irp, aco_params, false); 
        check_feasibility(ind, irp); 
        
        if (ind.is_feasible) {
            calculate_total_cost(ind, irp);
        } else {
            ind.fitness = LARGE_COST;
        }
        pop.push_back(ind);
    }
    std::cout << "População inicial educada e avaliada." << std::endl;
    
    std::sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b){
        return a.fitness < b.fitness;
    });
    Individual bestOverall = pop.front();
    
    std::cout << "\nIniciando loop do GA para " << ga_params.nGen << " gerações..." << std::endl;
    for (int gen = 0; gen < ga_params.nGen; ++gen) {
        std::cout<<std::endl;
        std::cout<<std::endl;
        std::cout<<"Iniciando geração "<<gen+1<<std::endl;
        vector<Individual> newPop;
        newPop.reserve(ga_params.popSize);

        // 1. Elitismo
        int num_elites = (int)(ga_params.popSize * 0.10);
        for(int i = 0; i < num_elites && i < pop.size(); ++i) {
            newPop.push_back(pop[i]);
        }
        std::cout<<"Elitismo completo com "<<num_elites<<" indivíduos."<<std::endl;
        // 2. Geração de Filhos
        int children_target_count = num_elites + (int)(ga_params.popSize * 0.70);
        int crossover_tries = 0; 

        while (newPop.size() < children_target_count) {
            std::cout<<"Gerando filhos: "<<newPop.size()+1<<" de "<<children_target_count<<std::endl;   
            Individual parent1 = tournamentSelect(pop, ga_params.tournamentK);
            Individual parent2 = tournamentSelect(pop, ga_params.tournamentK);
            
            std::pair<Individual, Individual> children;
            if (randreal() < ga_params.pCrossover) {
                children = (randreal() < 0.5) ? 
                           one_point_crossover_customer(parent1, parent2, irp) :
                           two_point_crossover_customer(parent1, parent2, irp);
            } else {
                children = {parent1, parent2};
            }
            
            // --- Processa o Filho 1 ---
            
            build_routes_for_individual(children.first, irp, aco_params);
            if(gen % 50 == 0)busca_local(children.first, irp, aco_params, false);
            check_feasibility(children.first, irp); 
            
            if (children.first.is_feasible) {
                calculate_total_cost(children.first, irp); 
                newPop.push_back(children.first);
            }
            else {
                children.first.fitness = LARGE_COST;
                newPop.push_back(children.first);
            }
            
            if (newPop.size() >= children_target_count) break;


            build_routes_for_individual(children.second, irp, aco_params);
            if(gen % 50 == 0)busca_local(children.second, irp, aco_params, false);
            check_feasibility(children.second, irp);
            
            if (children.second.is_feasible) {
                calculate_total_cost(children.second, irp); 
                newPop.push_back(children.second);
            }
            else {
                children.second.fitness = LARGE_COST;
                newPop.push_back(children.second);
        }

        // 3. Geração de novas soluções
        while (newPop.size() < ga_params.popSize) {
            Individual ind = make_new_heuristic_individual(irp, aco_params);
            build_routes_for_individual(ind, irp, aco_params);
            if(gen % 50 == 0)busca_local(ind, irp, aco_params, false);

            check_feasibility(ind, irp);
            
            if (ind.is_feasible) {
                calculate_total_cost(ind, irp);
            } else {
                ind.fitness = LARGE_COST;
            }
            newPop.push_back(ind);
        }
        
        pop.swap(newPop);

        std::sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b){
            return a.fitness < b.fitness;
        });
        
        if (pop.front().fitness < bestOverall.fitness) {
            bestOverall = pop.front();
        }
        
        // Log da Geração
        double min_fit = pop.front().fitness;
        double avg_fit = 0.0;
        int feasible_count = 0;
        for(const auto& ind : pop) {
            if(ind.is_feasible) {
                avg_fit += ind.fitness;
                feasible_count++;
            }
        }
        if (feasible_count > 0) avg_fit /= feasible_count;
        
        if (true) {
            std::cout << "Gen " << std::setw(4) << gen + 1 << "/" << ga_params.nGen
                      << " | Fact.: " << std::setw(3) << feasible_count << "/" << (int)pop.size()
                      << " | Melhor: " << std::fixed << std::setprecision(2) << min_fit
                      << " | Média(fact): " << avg_fit
                      << " | Global: " << bestOverall.fitness << std::endl;
                      
        }
    }
    
}  
    return bestOverall;
}




void printDeliveriesMatrix(const Individual& ind, const IRP& irp) {
    std::cout << "Matriz de entregas (linhas=periodos, colunas=clientes):\n";
    for (int t = 0; t < irp.nPeriods; ++t) {
        for (int c = 0; c < irp.nCustomers; ++c) {
            std::cout << ind.deliveries[t][c] << (c + 1 < irp.nCustomers ? ' ' : '\n');
        }
    }
}

void exportAndPlotRoutes(
    const IRP& irp,
    const Individual& best,
    const ACO_Params& acoParams,
    const std::string& dataFilename,
    const std::string& pyScript) 
{
    std::ofstream out(dataFilename);
    if (!out.is_open()) return;

    int T = irp.nPeriods, nCust = irp.nCustomers, nVeh = irp.nVehicles;
    out << T << " " << nCust << " " << nVeh << " " << irp.Capacity << "\n";
    out << "DEPOT " << irp.depots[0].x << " " << irp.depots[0].y << "\n";

    vector<int> inv(nCust);
    for (int c = 0; c < nCust; ++c) inv[c] = irp.customers[c].initialInv;

    for (int t = 0; t < T; ++t) {
        out << "PERIOD " << t << "\n";
        out << "CUSTOMERS\n";
        const vector<int>& del = best.deliveries[t];
        for (int c = 0; c < nCust; ++c) {
            long inv_after = (long)inv[c] + (long)del[c] - (long)irp.customers[c].demand[t];
            out << (c+1) << " " << irp.customers[c].x << " " << irp.customers[c].y << " "
                << inv[c] << " " << del[c] << " " << irp.customers[c].demand[t] << " " << inv_after << "\n";
        }

        out << "ROUTES\n";
        if (!best.routes_per_period[t].empty()) {
            for (const Route& route : best.routes_per_period[t]) {
                for (size_t i = 0; i < route.visits.size(); ++i) {
                    out << route.visits[i] << (i + 1 < route.visits.size() ? " " : "");
                }
                out << "\n";
            }
        }
        out << "END_ROUTES\n";

        for (int c = 0; c < nCust; ++c) {
            inv[c] += del[c] - irp.customers[c].demand[t];
        }
    }
    out.close();

    std::ostringstream cmd;
    cmd << "python3 " << pyScript << " " << dataFilename << " > /dev/null 2>&1";
    std::system(cmd.str().c_str());
}

