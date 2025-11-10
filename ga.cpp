#include "ds_operator.hpp"
#include "ga.hpp"
#include "utils.hpp"
#include "evaluation.hpp" 
#include "aco.hpp"
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


void build_routes_for_individual(Individual& ind, const IRP& irp, const ACO_Params& aco_params) {
    ind.routes_per_period.assign(irp.nPeriods, std::vector<Route>());
    const double LARGE_COST = 1e18;

    for (int t = 0; t < irp.nPeriods; ++t) {
        long period_delivery_sum = std::accumulate(ind.deliveries[t].begin(), ind.deliveries[t].end(), 0L);
        
        if (period_delivery_sum > 0) {
            ACO_Result ares = runACO_for_period(irp, ind.deliveries[t], aco_params, false);
            
            // Armazena as rotas. Se o ACO falhou (custo alto), 
            // as rotas estarão vazias ou inválidas, e 'check_feasibility' vai pegar.
            if(ares.bestCost < LARGE_COST / 10.0) {
                ind.routes_per_period[t] = ares.bestRoutes;
            } else {
                ind.routes_per_period[t].clear(); // Garante que está vazio se o ACO falhou
            }
        }
    }
}


Individual run_genetic_algorithm(const IRP& irp, const GA_Params& ga_params, const ACO_Params& aco_params, bool verbose) {
    
    const double LARGE_COST = 1e18;
    vector<Individual> pop;
    pop.reserve(ga_params.popSize);
    
    std::cout << "Inicializando população..." << std::endl;
    for (int i = 0; i < ga_params.popSize; ++i) {
        Individual ind = make_new_heuristic_individual(irp);
        
        build_routes_for_individual(ind, irp, aco_params);
        check_feasibility(ind, irp); 
        
        if (ind.is_feasible) {
            calculate_total_cost(ind, irp);
            if (verbose) std::cout << "Educando Indivíduo Inicial " << i << "...\n";
            busca_local(ind, irp, aco_params, false); 
        } else {
            ind.fitness = LARGE_COST;
            busca_local(ind, irp, aco_params, false); 
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
            check_feasibility(children.first, irp); 
            
            if (children.first.is_feasible) {
                calculate_total_cost(children.first, irp); 
                // <-- MUDANÇA: "Educa" o filho
                if (randreal() < 0.2)busca_local(children.first, irp, aco_params, false);
                newPop.push_back(children.first);
            }
            else {
                children.first.fitness = LARGE_COST;
                busca_local(children.first, irp, aco_params, false);
                newPop.push_back(children.first);
            }
            
            if (newPop.size() >= children_target_count) break;


            build_routes_for_individual(children.second, irp, aco_params);
            check_feasibility(children.second, irp);
            
            if (children.second.is_feasible) {
                calculate_total_cost(children.second, irp); 
                // <-- MUDANÇA: "Educa" o filho
                if (randreal() < 0.2)busca_local(children.second, irp, aco_params, false);
                newPop.push_back(children.second);
            }
            else {
                children.second.fitness = LARGE_COST;
                busca_local(children.second, irp, aco_params, false);
                newPop.push_back(children.second);
        }

        // 3. Geração de novas soluções
        while (newPop.size() < ga_params.popSize) {
            Individual ind = make_new_heuristic_individual(irp);
            build_routes_for_individual(ind, irp, aco_params);
            check_feasibility(ind, irp);
            
            if (ind.is_feasible) {
                calculate_total_cost(ind, irp);
                if (randreal() < 0.2)busca_local(ind, irp, aco_params, false);
            } else {
                ind.fitness = LARGE_COST;
                busca_local(ind, irp, aco_params, false);
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

