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


/*
 * Crossover Baseado em Tempo (Adaptado de Zhao et al. 2025)
 * Em vez de cortar a matriz verticalmente (por clientes),
 * este operador corta horizontalmente (por dias).
 * * Lógica:
 * 1. Embaralha os dias.
 * 2. Define pontos de corte.
 * 3. Filho herda o "plano de entregas completo do dia" do Pai A ou Pai B.
 */
std::pair<Individual, Individual> crossover_time_based(const Individual& a, const Individual& b, const IRP& irp) {
    Individual child1(irp.nPeriods, irp.nCustomers);
    Individual child2(irp.nPeriods, irp.nCustomers);

    // 1. Cria lista de dias e embaralha
    std::vector<int> days(irp.nPeriods);
    std::iota(days.begin(), days.end(), 0);
    std::shuffle(days.begin(), days.end(), rng);

    // 2. Define ponto de corte (vamos usar 1 ponto para simplificar, ou 2 como no artigo)
    // Zhao usa 2 pontos para criar 3 zonas. Vamos replicar a lógica de 3 zonas.
    int j1 = randint(0, irp.nPeriods);
    int j2 = randint(0, irp.nPeriods);
    if (j1 > j2) std::swap(j1, j2);

    for (int idx = 0; idx < irp.nPeriods; ++idx) {
        int t = days[idx]; // O dia real correspondente a este índice embaralhado

        // Lógica de herança baseada nas zonas
        // Zona 1 (idx < j1) ou Zona 3 (idx >= j2): Herança Direta
        // Zona 2 (j1 <= idx < j2): Troca
        
        const Individual* source1_for_c1;
        const Individual* source2_for_c1;
        
        if (idx < j1 || idx >= j2) {
            // Zonas externas: Filho 1 herda de A, Filho 2 herda de B
            source1_for_c1 = &a;
            source2_for_c1 = &b;
        } else {
            // Zona interna: Inverte (Crossover)
            source1_for_c1 = &b;
            source2_for_c1 = &a;
        }

        // Copia o plano de entregas DO DIA INTEIRO
        for (int c = 0; c < irp.nCustomers; ++c) {
            child1.deliveries[t][c] = source1_for_c1->deliveries[t][c];
            child2.deliveries[t][c] = source2_for_c1->deliveries[t][c];
        }
    }
    
    // Nota: A factibilidade de estoque (stock-out) pode ser quebrada aqui
    // porque estamos misturando dias de planos diferentes.
    // O pipeline do GA (check_feasibility -> discard) cuidará disso,
    // mas este crossover pode ter taxa de rejeição mais alta que os outros.
    
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



// Individual make_new_heuristic_individual(const IRP& irp) {
//     Individual ind(irp.nPeriods, irp.nCustomers);
//     vector<long> current_inv(irp.nCustomers);
//     for(int i = 0; i < irp.nCustomers; ++i) current_inv[i] = irp.customers[i].initialInv;

//     long fleet_capacity = (long)irp.nVehicles * irp.Capacity;

//     for (int t = 0; t < irp.nPeriods; ++t) {
//         long period_load = 0;
//         vector<int> mandatory_custs, optional_custs;

//         for (int c = 0; c < irp.nCustomers; ++c) {
//             if (current_inv[c] < irp.customers[c].demand[t]) {
//                 mandatory_custs.push_back(c);
//             } else {
//                 optional_custs.push_back(c);
//             }
//         }

//         // Atende clientes obrigatórios
//         for (int c : mandatory_custs) {
//             long needed = irp.customers[c].demand[t] - current_inv[c];
//             long space = irp.customers[c].maxLevelInv - (current_inv[c] + needed);
//             int extra = (space > 0) ? randint(0, std::min((long)20, space)) : 0; // Adiciona um extra aleatório pequeno
//             int q = std::min((long)irp.Capacity, needed + extra);
            
//             if (period_load + q <= fleet_capacity) {
//                 ind.deliveries[t][c] = q;
//                 period_load += q;
//             }
//         }
        
//         // Atende clientes opcionais aleatoriamente
//         std::shuffle(optional_custs.begin(), optional_custs.end(), rng);
//         for (int c : optional_custs) {
//             if (period_load >= fleet_capacity) break;
//             if (randreal() < 0.3) { // Chance de atender um cliente opcional
//                 long space = irp.customers[c].maxLevelInv - current_inv[c];
//                 long max_q = std::min({space, (long)irp.Capacity, fleet_capacity - period_load});
//                 if (max_q > 0) {
//                     int q = randint(1, max_q);
//                     ind.deliveries[t][c] = q;
//                     period_load += q;
//                 }
//             }
//         }

//         // Atualiza o inventário para o próximo período
//         for (int c = 0; c < irp.nCustomers; ++c) {
//             current_inv[c] += (long)ind.deliveries[t][c] - irp.customers[c].demand[t];
//         }
//     }
//     return ind;
// }


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



// Individual make_new_heuristic_individual(const IRP& irp, const ACO_Params& aco_params) {
    
//     Individual ind(irp.nPeriods, irp.nCustomers);
    
//     // --- PASSO 1: Sortear 20% dos clientes para reservar ---
//     // int num_reserved = static_cast<int>(irp.nCustomers * 1);
//     double removal_rate = 0.10 + (randreal() * 0.40); // 0.10 a 0.50
//     removal_rate = 0.0f;
//     int num_reserved = static_cast<int>(irp.nCustomers * removal_rate);
//     //if (num_reserved == 0 && irp.nCustomers > 0) num_reserved = 1; // Garante pelo menos 1
    
//     vector<int> customer_indices(irp.nCustomers);
//     std::iota(customer_indices.begin(), customer_indices.end(), 0); // 0, 1, 2...
//     std::shuffle(customer_indices.begin(), customer_indices.end(), rng);
    
//     std::set<int> reserved_customers; // Clientes a serem inseridos depois (0-based)
//     for(int i=0; i < num_reserved; ++i) {
//         reserved_customers.insert(customer_indices[i]);
//     }
    
//     // --- PASSO 2: Construir solução para os 80% (Base) ---
//     vector<long> current_inv(irp.nCustomers);
//     for(int i = 0; i < irp.nCustomers; ++i) current_inv[i] = irp.customers[i].initialInv;
//     long fleet_capacity = (long)irp.nVehicles * irp.Capacity;

//     for (int t = 0; t < irp.nPeriods; ++t) {
//         vector<int> min_deliveries_today(irp.nCustomers, 0);
//         long period_load_min = 0;

//         for (int c = 0; c < irp.nCustomers; ++c) {
//             // PULA CLIENTES RESERVADOS
//             if (reserved_customers.count(c)) continue; 

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
//             // (Simula para o próximo dia)
//             for (int c = 0; c < irp.nCustomers; ++c) {
//                 if (reserved_customers.count(c)) continue; // Pula
//                 current_inv[c] += (long)ind.deliveries[t][c] - irp.customers[c].demand[t];
//                 if(current_inv[c] < irp.customers[c].minLevelInv) current_inv[c] = irp.customers[c].minLevelInv;
//             }
//             continue; 
//         }

//         // Entregas oportunistas (apenas para os 80%)
//         long remaining_fleet_cap = fleet_capacity - current_period_load;
//         vector<int> customer_order;
//         for(int c=0; c < irp.nCustomers; ++c) if(!reserved_customers.count(c)) customer_order.push_back(c);
//         std::shuffle(customer_order.begin(), customer_order.end(), rng);

//         for (int c : customer_order) {
//             if (remaining_fleet_cap <= 0) break;
//             if (randreal() < 0.15) { 
//                 int N = randint(1, 3);
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
//         } 
        
//         // Simula inventário para o próximo dia (apenas 80%)
//         for (int c = 0; c < irp.nCustomers; ++c) {
//             if (reserved_customers.count(c)) continue; // Pula
//             current_inv[c] += (long)ind.deliveries[t][c] - irp.customers[c].demand[t];
//             if(current_inv[c] < irp.customers[c].minLevelInv) {
//                 current_inv[c] = irp.customers[c].minLevelInv;
//             }
//         }
//     } // Fim do loop 't' (plano 80% está feito)
    
//     // --- PASSO 3: Criar as rotas base (sem os clientes reservados) ---
//     build_routes_for_individual(ind, irp, aco_params);
//     // (Neste ponto, 'ind' é uma solução parcial factível para 80% dos clientes)

//     // --- PASSO 4: Chamar a "busca local" (lógica Gurobi) para INSERIR os clientes ---
    
 
              
//     for (int c_id_internal : reserved_customers) {
        
//         // 1. Calcula dados de reinserção (custos F_t(q)) para 'c'
//         //    com base nas rotas ATUAIS (que ainda não o contêm)
//         ReinsertionData data = calculate_reinsertion_data(ind, irp, aco_params, c_id_internal);
        
//         // 2. Chama Gurobi para encontrar o plano de inventário ótimo para 'c'
//         GurobiInventoryResult result = computePerfectInventory_Gurobi(
//             irp.customers[c_id_internal],
//             irp.depots[0],
//             {}, 
//             data.insertion_cost_curves,
//             irp.nPeriods,
//             (double)irp.Capacity 
//         );

//         if (result.totalCost >= PENALTY_COST * 0.9) {
//             // Gurobi falhou em encontrar um plano de inserção.
//             // O indivíduo provavelmente será infactível.
//              std::cerr << "AVISO: Falha do Gurobi ao inserir cliente " << (c_id_internal+1) << "\n";
//             continue; // Tenta o próximo cliente
//         }

//         // 3. ATUALIZA O GENÓTIPO: Adiciona o plano do Gurobi ao 'ind.deliveries'
//         for(int t=0; t < irp.nPeriods; ++t) {
//             int delivery_q = static_cast<int>(std::round(result.q[t]));
//             if (delivery_q > 0) {
//                 ind.deliveries[t][c_id_internal] = delivery_q;
//             }
//         }
        
//         // 4. ATUALIZA O FENÓTIPO: Reconstrói as rotas
//         // (Agora incluindo o cliente 'c')
//         build_routes_for_individual(ind, irp, aco_params);
        
//         // (Verificação de factibilidade não é necessária aqui, 
//         // pois a próxima inserção (passo 1) recalculará tudo)
//     }

//     return ind;
// }


Individual make_new_heuristic_individual(const IRP& irp, const ACO_Params& aco_params) {
    
    Individual ind(irp.nPeriods, irp.nCustomers);
    
    vector<long> current_inv(irp.nCustomers);
    for(int i = 0; i < irp.nCustomers; ++i) current_inv[i] = irp.customers[i].initialInv;

    double service_probability = 0.3 + (randreal() * 0.5);

    for (int t = 0; t < irp.nPeriods; ++t) {
        // 1. Identificar Urgentes e Não-Urgentes
        vector<int> urgent_customers;
        vector<int> non_urgent_customers;
        
        for (int c = 0; c < irp.nCustomers; ++c) {
            long inv_next = current_inv[c] - irp.customers[c].demand[t];
            if (inv_next < irp.customers[c].minLevelInv) {
                urgent_customers.push_back(c);
            } else {
                if (current_inv[c] < irp.customers[c].maxLevelInv) {
                    non_urgent_customers.push_back(c);
                }
            }
        }

        // 2. Inicializar "Veículos Virtuais"
        // (Agora o compilador sabe o que é 'VehicleBin')
        vector<VehicleBin> fleet(irp.nVehicles);
        for(int k=0; k<irp.nVehicles; ++k) {
            fleet[k] = {k, 0, irp.Capacity};
        }

        std::shuffle(urgent_customers.begin(), urgent_customers.end(), rng);
        std::shuffle(non_urgent_customers.begin(), non_urgent_customers.end(), rng);

        // 3. Atender Urgentes (Mandatório)
        for (int c : urgent_customers) {
            long demand_today = irp.customers[c].demand[t];
            long min_needed = irp.customers[c].minLevelInv - (current_inv[c] - demand_today);
            long max_can_receive = irp.customers[c].maxLevelInv - current_inv[c];
            
            int best_v = -1;
            int min_waste = irp.Capacity + 1;

            // Tenta alocar a quantidade MÍNIMA (Best Fit)
            for(int k=0; k<irp.nVehicles; ++k) {
                if (fleet[k].current_load + min_needed <= fleet[k].capacity) {
                    int waste = fleet[k].capacity - (fleet[k].current_load + min_needed);
                    if (waste < min_waste) {
                        min_waste = waste;
                        best_v = k;
                    }
                }
            }

            if (best_v != -1) {
                // Coube. Tenta entregar mais (Order-Up-To ou encher o veículo)
                long space_in_vehicle = fleet[best_v].capacity - fleet[best_v].current_load;
                long delivery = std::min(max_can_receive, space_in_vehicle);
                // Garante que pelo menos o mínimo seja entregue
                delivery = std::max(delivery, min_needed); 
                
                ind.deliveries[t][c] = delivery;
                fleet[best_v].current_load += delivery;
            } else {
                // FALHA CRÍTICA (Plano mínimo não cabe na frota)
                // Força a entrega no veículo mais vazio (cria infactibilidade)
                int v_empty = 0; 
                for(int k=1; k<irp.nVehicles; ++k) if(fleet[k].current_load < fleet[v_empty].current_load) v_empty = k;
                
                ind.deliveries[t][c] = min_needed; 
                fleet[v_empty].current_load += min_needed; 
            }
        }

        // 4. Atender Não-Urgentes (Oportunista)
        for (int c : non_urgent_customers) {
            if (randreal() > service_probability) continue;

            long max_can_receive = irp.customers[c].maxLevelInv - current_inv[c];
            if (max_can_receive <= 0) continue;

            // Encontra um veículo com espaço (First Fit)
            int best_v = -1;
            for(int k=0; k<irp.nVehicles; ++k) {
                if (fleet[k].current_load < fleet[k].capacity) {
                    best_v = k;
                    break;
                }
            }

            if (best_v != -1) {
                long space_in_vehicle = fleet[best_v].capacity - fleet[best_v].current_load;
                long max_possible = std::min(max_can_receive, space_in_vehicle);
                
                if (max_possible > 0) {
                    double ratio = 0.5 + (randreal() * 0.5); 
                    long delivery = std::max(1L, (long)(max_possible * ratio));
                    
                    ind.deliveries[t][c] = delivery;
                    fleet[best_v].current_load += delivery;
                }
            }
        }

        // 5. Atualiza Estoque para t+1
        for (int c = 0; c < irp.nCustomers; ++c) {
            current_inv[c] += (long)ind.deliveries[t][c] - irp.customers[c].demand[t];
            if (current_inv[c] < irp.customers[c].minLevelInv) 
                current_inv[c] = irp.customers[c].minLevelInv;
        }
    }

    return ind;
}

// Individual run_genetic_algorithm(const IRP& irp, const GA_Params& ga_params, const ACO_Params& aco_params, bool verbose) {
    
//     const double LARGE_COST = 1e18;
//     vector<Individual> pop;
//     pop.reserve(ga_params.popSize);
    
//     std::cout << "Inicializando população..." << std::endl;
//     for (int i = 0; i < ga_params.popSize; ++i) {
//         Individual ind = make_new_heuristic_individual(irp, aco_params);
        
//         build_routes_for_individual(ind, irp, aco_params);
//        // busca_local(ind, irp, aco_params, false); 
//         check_feasibility(ind, irp); 
        
//         if (ind.is_feasible) {
//             calculate_total_cost(ind, irp);
//         } else {
//             ind.fitness = LARGE_COST;
//         }
//         pop.push_back(ind);
//     }
//     std::cout << "População inicial educada e avaliada." << std::endl;
    
//     std::sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b){
//         return a.fitness < b.fitness;
//     });
//     Individual bestOverall = pop.front();
    
//     std::cout << "\nIniciando loop do GA para " << ga_params.nGen << " gerações..." << std::endl;
//     for (int gen = 0; gen < ga_params.nGen; ++gen) {
//         std::cout<<std::endl;
//         std::cout<<std::endl;
//         std::cout<<"Iniciando geração "<<gen+1<<std::endl;
//         vector<Individual> newPop;
//         newPop.reserve(ga_params.popSize);

//         // 1. Elitismo
//         int num_elites = (int)(ga_params.popSize * 0.10);
//         for(int i = 0; i < num_elites && i < pop.size(); ++i) {
//             newPop.push_back(pop[i]);
//         }
//         std::cout<<"Elitismo completo com "<<num_elites<<" indivíduos."<<std::endl;
//         // 2. Geração de Filhos
//         int children_target_count = num_elites + (int)(ga_params.popSize * 0.70);
//         int crossover_tries = 0; 

//         while (newPop.size() < children_target_count) {
//             std::cout<<"Gerando filhos: "<<newPop.size()+1<<" de "<<children_target_count<<std::endl;   
//             Individual parent1 = tournamentSelect(pop, ga_params.tournamentK);
//             Individual parent2 = tournamentSelect(pop, ga_params.tournamentK);
            
//             std::pair<Individual, Individual> children;
//             if (randreal() < ga_params.pCrossover) {
//                 children = (randreal() < 0.5) ? 
//                            one_point_crossover_customer(parent1, parent2, irp) :
//                            two_point_crossover_customer(parent1, parent2, irp);
//             } else {
//                 children = {parent1, parent2};
//             }
            
//             // --- Processa o Filho 1 ---
            
//             build_routes_for_individual(children.first, irp, aco_params);
//             if(gen % 50 == 0)busca_local(children.first, irp, aco_params, false);
//             check_feasibility(children.first, irp); 
            
//             if (children.first.is_feasible) {
//                 calculate_total_cost(children.first, irp); 
//                 newPop.push_back(children.first);
//             }
//             else {
//                 children.first.fitness = LARGE_COST;
//                 newPop.push_back(children.first);
//             }
            
//             if (newPop.size() >= children_target_count) break;


//             build_routes_for_individual(children.second, irp, aco_params);
//             if(gen % 50 == 0)busca_local(children.second, irp, aco_params, false);
//             check_feasibility(children.second, irp);
            
//             if (children.second.is_feasible) {
//                 calculate_total_cost(children.second, irp); 
//                 newPop.push_back(children.second);
//             }
//             else {
//                 children.second.fitness = LARGE_COST;
//                 newPop.push_back(children.second);
//         }

//         // 3. Geração de novas soluções
//         while (newPop.size() < ga_params.popSize) {
//             Individual ind = make_new_heuristic_individual(irp, aco_params);
//             build_routes_for_individual(ind, irp, aco_params);
//             if(gen % 50 == 0)busca_local(ind, irp, aco_params, false);

//             check_feasibility(ind, irp);
            
//             if (ind.is_feasible) {
//                 calculate_total_cost(ind, irp);
//             } else {
//                 ind.fitness = LARGE_COST;
//             }
//             newPop.push_back(ind);
//         }
        
//         pop.swap(newPop);

//         std::sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b){
//             return a.fitness < b.fitness;
//         });
        
//         if (pop.front().fitness < bestOverall.fitness) {
//             bestOverall = pop.front();
//         }
        
//         // Log da Geração
//         double min_fit = pop.front().fitness;
//         double avg_fit = 0.0;
//         int feasible_count = 0;
//         for(const auto& ind : pop) {
//             if(ind.is_feasible) {
//                 avg_fit += ind.fitness;
//                 feasible_count++;
//             }
//         }
//         if (feasible_count > 0) avg_fit /= feasible_count;
        
//         if (true) {
//             std::cout << "Gen " << std::setw(4) << gen + 1 << "/" << ga_params.nGen
//                       << " | Fact.: " << std::setw(3) << feasible_count << "/" << (int)pop.size()
//                       << " | Melhor: " << std::fixed << std::setprecision(2) << min_fit
//                       << " | Média(fact): " << avg_fit
//                       << " | Global: " << bestOverall.fitness << std::endl;
                      
//         }
//     }
    
// }  
//     return bestOverall;
// }



// Individual run_genetic_algorithm(const IRP& irp, const GA_Params& ga_params, const ACO_Params& aco_params, bool verbose) {
    
//     const double LARGE_COST = 1e18;
//     vector<Individual> pop;
//     pop.reserve(ga_params.popSize);
    
//     std::cout << "Inicializando população..." << std::endl;
//     for (int i = 0; i < ga_params.popSize; ++i) {
//         Individual ind = make_new_heuristic_individual(irp);
//         build_routes_for_individual(ind, irp, aco_params);
//         check_feasibility(ind, irp); 
        
//         if (ind.is_feasible) {
//             calculate_total_cost(ind, irp);
//             // --- EDUCAÇÃO (Busca Leve) ---
//             if (verbose) std::cout << "Educando (Leve) Indivíduo Inicial " << i << "...\n";
//             run_simple_reinsertion_search(ind, irp, aco_params, true); 
//         } else {
//             ind.fitness = LARGE_COST;
//         }
//         pop.push_back(ind);
//     }
//     std::cout << "População inicial educada e avaliada." << std::endl;
    
//     std::sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b){
//         return a.fitness < b.fitness;
//     });
//     Individual bestOverall = pop.front();
    
//     std::cout << "\nIniciando loop do GA para " << ga_params.nGen << " gerações..." << std::endl;
//     for (int gen = 0; gen < ga_params.nGen; ++gen) {
//         vector<Individual> newPop;
//         newPop.reserve(ga_params.popSize);

//         // 1. Elitismo
//         std::cout<<"Elitismo"<<std::endl;
//         int num_elites = (int)(ga_params.popSize * 0.10);
//         for(int i = 0; i < num_elites && i < pop.size(); ++i) {
//             newPop.push_back(pop[i]);
//         }
//         std::cout<<"Elitismo finalizado para geração "<<gen<<std::endl;

//         // 2. Geração de Filhos
//         int children_target_count = num_elites + (int)(ga_params.popSize * 0.70);
//         int crossover_tries = 0; 
//         std::cout<<"Crossover"<<std::endl;

//         while (newPop.size() < children_target_count ) {
//             Individual parent1 = tournamentSelect(pop, ga_params.tournamentK);
//             Individual parent2 = tournamentSelect(pop, ga_params.tournamentK);
            
//             std::pair<Individual, Individual> children;
//             if (randreal() < ga_params.pCrossover) {
//                             double r = randreal();
//                             if (r < 0.33) {
//                                 children = one_point_crossover_customer(parent1, parent2, irp);
//                             } else if (r < 0.66) {
//                                 children = two_point_crossover_customer(parent1, parent2, irp);
//                             } else {
//                                 children = crossover_time_based(parent1, parent2, irp);
//                             }
//                         } else {
//                             children = {parent1, parent2};
//                         }
    
//             build_routes_for_individual(children.first, irp, aco_params);
//             check_feasibility(children.first, irp); 
//             std::cout<<"Filho 1"<<std::endl;

//             if (children.first.is_feasible) {
//                 calculate_total_cost(children.first, irp); 
//                 // --- EDUCAÇÃO (Busca Leve) ---
//                 run_simple_reinsertion_search(children.first, irp, aco_params, true);
//                 newPop.push_back(children.first);
//             }
//             if (newPop.size() >= children_target_count) break;

//             // --- Processa o Filho 2 ---
//             // advance_portion_mutation(children.second, irp, ga_params.pMutation);
//             build_routes_for_individual(children.second, irp, aco_params);
//             check_feasibility(children.second, irp);
            
//             if (children.second.is_feasible) {
//                 calculate_total_cost(children.second, irp); 
//                 // --- EDUCAÇÃO (Busca Leve) ---
//                 run_simple_reinsertion_search(children.second, irp, aco_params, true);
//                 newPop.push_back(children.second);
//             }
//         }
//         std::cout<<"Filho 2"<<std::endl;
//         // 3. Geração de Imigrantes
//         while (newPop.size() < ga_params.popSize) {
//             Individual ind = make_new_heuristic_individual(irp);
//             build_routes_for_individual(ind, irp, aco_params);
//             check_feasibility(ind, irp);
//             if (ind.is_feasible) {
//                 calculate_total_cost(ind, irp);
//                 // --- EDUCAÇÃO (Busca Leve) ---
//                 run_simple_reinsertion_search(ind, irp, aco_params, true);
//             } else {
//                 ind.fitness = LARGE_COST;
//             }
//             newPop.push_back(ind);
//         }
        
//         pop.swap(newPop);

//         // --- MUDANÇA: BUSCA LOCAL "PESADA" (a cada 100 gerações) ---
//         if (gen > 0 && gen % 100 == 0) {
//             if (verbose) {
//                 std::cout << "\n--- EXECUTANDO BUSCA LOCAL PESADA (VND-DS) NA POPULAÇÃO ---\n";
//             }
//             for (Individual& ind : pop) {
//                 if (ind.is_feasible) {
//                     // Chama a busca "pesada" (antiga 'sua_futura_busca_local')
//                     run_vnd_ds_operator(ind, irp, aco_params, false);
//                 }
//             }
//             if (verbose) {
//                 std::cout << "--- BUSCA PESADA CONCLUÍDA ---\n";
//             }
//         }
        
//         // Ordena e atualiza o melhor
//         std::sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b){
//             return a.fitness < b.fitness;
//         });
//         if (pop.front().fitness < bestOverall.fitness) {
//             bestOverall = pop.front();
//         }
        
//         // Log da Geração
//         double min_fit = pop.front().fitness;
//         double avg_fit = 0.0;
//         int feasible_count = 0;
//         for(const auto& ind : pop) {
//             if(ind.is_feasible) { avg_fit += ind.fitness; feasible_count++; }
//         }
//         if (feasible_count > 0) avg_fit /= feasible_count;
        
//         if (verbose || (gen % 10 == 0) || (gen == ga_params.nGen - 1)) {
//             std::cout << "Gen " << std::setw(4) << gen + 1 << "/" << ga_params.nGen
//                       << " | Fact.: " << std::setw(3) << feasible_count << "/" << (int)pop.size()
//                       << " | Melhor: " << std::fixed << std::setprecision(2) << min_fit
//                       << " | Média(fact): " << avg_fit
//                       << " | Global: " << bestOverall.fitness << "\n";
//         }
//     }
    
//     return bestOverall;
// }




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

/**
 * @brief Calcula a distância de diversidade (Eq. 21 Zhao et al. 2025)
 * Delta(P1, P2) = (1/n) * sum( 1 se padroes de visita diferem )
 */
double calculate_diversity_distance(const Individual& a, const Individual& b, int nCustomers, int nPeriods) {
    double diff_count = 0.0;

    for (int c = 0; c < nCustomers; ++c) {
        bool pattern_differs = false;
        for (int t = 0; t < nPeriods; ++t) {
            // Verifica se o status de visita (tem entrega ou não) é diferente
            bool visit_a = (a.deliveries[t][c] > 0);
            bool visit_b = (b.deliveries[t][c] > 0);
            
            if (visit_a != visit_b) {
                pattern_differs = true;
                break; // Basta um dia diferente para o padrão ser diferente
            }
        }
        if (pattern_differs) {
            diff_count += 1.0;
        }
    }
    return diff_count / (double)nCustomers;
}

/**
 * @brief Atualiza os ranks e o Biased Fitness da população.
 */
void update_biased_fitness(std::vector<Individual>& pop, int nCustomers, int nPeriods, int num_elites) {
    int pop_size = pop.size();
    if (pop_size == 0) return;

    // 1. RANK DE CUSTO
    // Ordena por custo (menor é melhor)
    std::sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b) {
        return a.fitness < b.fitness;
    });
    // Atribui rank de custo (0..N-1)
    for (int i = 0; i < pop_size; ++i) {
        pop[i].rank_cost = i;
    }

    // 2. RANK DE DIVERSIDADE
    // Calcula a distância média de cada indivíduo para os outros
    // (No HGS original, usa-se os k vizinhos mais próximos, aqui usamos média geral para simplificar
    // ou podemos usar os 50% mais próximos para ser mais fiel ao Vidal 2012)
    int n_closest = std::max(1, pop_size / 2); 

    for (int i = 0; i < pop_size; ++i) {
        std::vector<double> distances;
        distances.reserve(pop_size);
        for (int j = 0; j < pop_size; ++j) {
            if (i == j) continue;
            distances.push_back(calculate_diversity_distance(pop[i], pop[j], nCustomers, nPeriods));
        }
        // Pega a média dos n_closest mais próximos (menor distância = menos diverso)
        std::sort(distances.begin(), distances.end()); // Crescente
        double sum_dist = 0.0;
        for(int k=0; k < std::min((int)distances.size(), n_closest); ++k) {
            sum_dist += distances[k];
        }
        pop[i].diversity_contribution = (distances.empty()) ? 0.0 : (sum_dist / n_closest);
    }

    // Ordena por diversidade (MAIOR distância é melhor/mais diverso -> rank menor)
    // Queremos rank 0 para o mais diverso.
    std::sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b) {
        return a.diversity_contribution > b.diversity_contribution; // Decrescente
    });
    for (int i = 0; i < pop_size; ++i) {
        pop[i].rank_diversity = i;
    }

    // 3. BIASED FITNESS
    // BF = RankCost + (1 - nbElites/popSize) * RankDiversity
    // (Fórmula padrão HGS para equilibrar os dois objetivos)
    double w_div = 1.0 - ((double)num_elites / (double)pop_size);
    
    for (int i = 0; i < pop_size; ++i) {
        pop[i].biased_fitness = (double)pop[i].rank_cost + w_div * (double)pop[i].rank_diversity;
    }
}

/**
 * @brief Seleciona sobreviventes removendo os piores até atingir target_size.
 * Remove clones primeiro.
 */
void select_survivors(std::vector<Individual>& pop, int target_size, int nCustomers, int nPeriods, int num_elites) {
    
    // 1. Remove Clones (Soluções com distância 0 e mesmo custo)
    // Isso é crucial para evitar estagnação.
    bool clone_found = true;
    while (pop.size() > target_size && clone_found) {
        clone_found = false;
        update_biased_fitness(pop, nCustomers, nPeriods, num_elites);
        
        // Ordena por Biased Fitness (piores no final) para facilitar remoção,
        // mas mantendo o melhor custo protegido (rank_cost = 0 sempre sobrevive na lógica de sort do HGS,
        // mas aqui vamos garantir explicitamente).
        std::sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b) {
            return a.biased_fitness < b.biased_fitness; 
        });

        // Procura clones. Se achar, remove o que tiver pior biased fitness (o que está mais pro fim)
        // Compara todos com todos é O(N^2), ok para N=100.
        for (int i = pop.size() - 1; i >= 1; --i) {
            for (int j = 0; j < i; ++j) {
                if (std::abs(pop[i].fitness - pop[j].fitness) < 1e-4) { // Mesmo custo
                    if (calculate_diversity_distance(pop[i], pop[j], nCustomers, nPeriods) < 1e-9) {
                        // É clone. Remove 'i' (que tem pior ou igual biased fitness)
                        pop.erase(pop.begin() + i);
                        clone_found = true;
                        goto end_clone_search;
                    }
                }
            }
        }
        end_clone_search:;
    }

    // 2. Remove por Biased Fitness (Piores Indivíduos)
    while (pop.size() > target_size) {
        update_biased_fitness(pop, nCustomers, nPeriods, num_elites);
        
        // Ordena: Menor Biased Fitness (Melhor) -> Maior Biased Fitness (Pior)
        std::sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b) {
            return a.biased_fitness < b.biased_fitness;
        });

        // Remove o último (o pior), DESDE QUE não seja o melhor custo absoluto
        // (O melhor custo terá rank_cost=0, então dificilmente terá BF alto, mas por segurança)
        if (pop.back().rank_cost == 0) {
            // Caso patológico extremo onde o melhor custo é o pior em diversidade e a população é pequena.
            // Remove o penúltimo.
            pop.erase(pop.end() - 2);
        } else {
            pop.erase(pop.end() - 1);
        }
    }
}


// Individual run_genetic_algorithm(const IRP& irp, const GA_Params& ga_params, const ACO_Params& aco_params, bool verbose) {
    
//     const double LARGE_COST = 1e18;
//     vector<Individual> pop;
//     pop.reserve(ga_params.popSize * 2); // Reserva espaço extra para crescimento
    
//     // --- 1. INICIALIZAÇÃO (COM DIVERSIDADE) ---
//     std::cout << "Inicializando população..." << std::endl;
//     int initial_pop_size = ga_params.popSize; 
    
//     for (int i = 0; i < initial_pop_size; ++i) {
//         // Cria indivíduo "cru" com diversidade estrutural (10-50% removido)
//         Individual ind = make_new_heuristic_individual(irp, aco_params);
        
//         // Constrói rotas base
//         build_routes_for_individual(ind, irp, aco_params);
//         check_feasibility(ind, irp); 
        
//         if (ind.is_feasible) {
//             calculate_total_cost(ind, irp);
//             // MUDANÇA: NÃO aplica busca local aqui para manter diversidade inicial alta
//             // (Ou aplica uma busca MUITO leve se necessário, mas 'cru' é melhor para diversidade)
//         } else {
//             ind.fitness = LARGE_COST;
//         }
        
//         // Garante que só entram factíveis (ou tenta de novo)
//         if (ind.is_feasible) {
//             pop.push_back(ind);
//         } else {
//             i--; 
//         }
//         if (i % 10 == 0 && i > 0) std::cout << "." << std::flush;
//     }
//     std::cout << "\nPopulação inicial gerada." << std::endl;

//     // Encontra o melhor inicial
//     Individual bestOverall = pop[0];
//     for(const auto& ind : pop) if(ind.fitness < bestOverall.fitness) bestOverall = ind;

//     int num_elites = (int)(ga_params.popSize * 0.10); // Para cálculo do Biased Fitness

//     // --- 2. LOOP EVOLUCIONÁRIO (HGS) ---
//     std::cout << "\nIniciando loop do GA (HGS) para " << ga_params.nGen << " iterações..." << std::endl;
    
//     for (int gen = 0; gen < ga_params.nGen; ++gen) {
        
//         // a) Atualiza métricas de população (Ranks e Biased Fitness)
//         update_biased_fitness(pop, irp.nCustomers, irp.nPeriods, num_elites);

//         // b) Seleção de Pais (Torneio Binário via Biased Fitness)
//         int idx1 = randint(0, pop.size() - 1);
//         int idx2 = randint(0, pop.size() - 1);
//         const Individual& p1 = (pop[idx1].biased_fitness < pop[idx2].biased_fitness) ? pop[idx1] : pop[idx2];
        
//         idx1 = randint(0, pop.size() - 1);
//         idx2 = randint(0, pop.size() - 1);
//         const Individual& p2 = (pop[idx1].biased_fitness < pop[idx2].biased_fitness) ? pop[idx1] : pop[idx2];

//         // c) Crossover (Misto)
//         std::pair<Individual, Individual> children;
//         if (randreal() < ga_params.pCrossover) {
//             double r = randreal();
//             if (r < 0.33) {
//                 children = one_point_crossover_customer(p1, p2, irp);
//             } else if (r < 0.66) {
//                 children = two_point_crossover_customer(p1, p2, irp);
//             } else {
//                 // Novo Crossover Temporal (Zhao et al.)
//                 children = crossover_time_based(p1, p2, irp);
//             }
//         } else {
//             children = {p1, p2};
//         }

//         // d) Processamento do Filho 1 (Educação e Inserção)
//         build_routes_for_individual(children.first, irp, aco_params);
//         check_feasibility(children.first, irp);
//         if (children.first.is_feasible) {
//             calculate_total_cost(children.first, irp);
            
//             // Educação Leve (Sempre aplica nos filhos para garantir qualidade local)
//             run_simple_reinsertion_search(children.first, irp, aco_params, false); 
            
//             if (children.first.fitness < bestOverall.fitness) {
//                 bestOverall = children.first;
//                 std::cout << ">>> Nova Melhor Solução (Filho 1): " << std::fixed << std::setprecision(2) << bestOverall.fitness << "\n";
//             }
//             pop.push_back(children.first); // População cresce temporariamente
//         }

//         // e) Processamento do Filho 2
//         build_routes_for_individual(children.second, irp, aco_params);
//         check_feasibility(children.second, irp);
//         if (children.second.is_feasible) {
//             calculate_total_cost(children.second, irp);
            
//             run_simple_reinsertion_search(children.second, irp, aco_params, false);
            
//             if (children.second.fitness < bestOverall.fitness) {
//                 bestOverall = children.second;
//                 std::cout << ">>> Nova Melhor Solução (Filho 2): " << std::fixed << std::setprecision(2) << bestOverall.fitness << "\n";
//             }
//             pop.push_back(children.second);
//         }

//         // f) Seleção de Sobreviventes (Corte da População)
//         // Mantém tamanho constante removendo clones e piores (Biased Fitness)
//         if (pop.size() > ga_params.popSize) {
//             select_survivors(pop, ga_params.popSize, irp.nCustomers, irp.nPeriods, num_elites);
//         }

//         // g) Busca Local Pesada (VND-DS) / Shaking
//         // A cada 100 gerações, aplica a busca pesada (com Sweep) na elite
//         if (gen > 0 && gen % 100 == 0) {
//             if (verbose) std::cout << "\n--- Gen " << gen << ": Executando Busca Pesada (VND-DS) na Elite ---\n";
            
//             // Reordena por custo puro para pegar a elite verdadeira
//             std::sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b){
//                 return a.fitness < b.fitness;
//             });

//             // Aplica nos Top 5
//             for(int i=0; i<std::min((int)pop.size(), 5); ++i) {
//                 run_vnd_ds_operator(pop[i], irp, aco_params, false);
//                 if (pop[i].fitness < bestOverall.fitness) {
//                     bestOverall = pop[i];
//                     std::cout << ">>> Nova Melhor Solução (Busca Pesada): " << bestOverall.fitness << "\n";
//                 }
//             }
//         }

//         // Log periódico
//         if (verbose || (gen % 50 == 0) || gen == ga_params.nGen - 1) {
//             double avg_fit = 0.0;
//             for(const auto& ind : pop) avg_fit += ind.fitness;
//             avg_fit /= pop.size();
            
//             std::cout << "Gen " << std::setw(4) << gen 
//                       << " | Pop: " << pop.size()
//                       << " | Melhor: " << std::fixed << std::setprecision(2) << bestOverall.fitness
//                       << " | Média: " << avg_fit << "\n";
//         }
//     }
    
//     // Análise final na melhor solução (com verbose)
//     std::cout << "\nOtimização final na melhor solução...\n";
//     run_vnd_ds_operator(bestOverall, irp, aco_params, true);

//     return bestOverall;
// }


// Individual run_genetic_algorithm(const IRP& irp, const GA_Params& ga_params, const ACO_Params& aco_params, bool verbose) {
    
//     const double LARGE_COST = 1e18;
//     vector<Individual> pop;
//     pop.reserve(ga_params.popSize);
    
//     std::cout << "Inicializando população..." << std::endl;
//     for (int i = 0; i < ga_params.popSize; ++i) {
//         Individual ind = make_new_heuristic_individual(irp, aco_params);
//         build_routes_for_individual(ind, irp, aco_params);
//         check_feasibility(ind, irp); 
        
//         if (ind.is_feasible) {
//             calculate_total_cost(ind, irp);
//             // --- MUDANÇA 2: NÃO APLICAR BUSCA LOCAL NA INICIALIZAÇÃO ---
//             // Deixar a população "crua" aumenta a diversidade inicial.
//             // A busca local (Gurobi) é tão forte que faria todos convergirem
//             // para o mesmo ótimo local imediatamente.
//         } else {
//             ind.fitness = LARGE_COST;
//         }
//         pop.push_back(ind);
        
//         if (i % 10 == 0) std::cout << "."; // Progresso
//     }
//     std::cout << "\nPopulação inicial criada." << std::endl;
    
//     std::sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b){
//         return a.fitness < b.fitness;
//     });
//     Individual bestOverall = pop.front();
    
//     std::cout << "\nIniciando loop do GA para " << ga_params.nGen << " gerações..." << std::endl;
//     for (int gen = 0; gen < ga_params.nGen; ++gen) {
//         vector<Individual> newPop;
//         newPop.reserve(ga_params.popSize);

//         // 1. Elitismo (10%)
//         int num_elites = (int)(ga_params.popSize * 0.10);
//         for(int i = 0; i < num_elites && i < pop.size(); ++i) {
//             newPop.push_back(pop[i]);
//         }

//         // 2. Geração de Filhos (Crossover)
//         int children_target_count = num_elites + (int)(ga_params.popSize * 0.70);
//         int crossover_tries = 0; 

//         while (newPop.size() < children_target_count && crossover_tries < ga_params.crossover_max_tries) {
//             Individual parent1 = tournamentSelect(pop, ga_params.tournamentK);
//             Individual parent2 = tournamentSelect(pop, ga_params.tournamentK);
            
//          std::pair<Individual, Individual> children;
//         if (randreal() < ga_params.pCrossover) {
//                         double r = randreal();
//                         if (r < 0.33) {
//                             children = one_point_crossover_customer(p1, p2, irp);
//                         } else if (r < 0.66) {
//                             children = two_point_crossover_customer(p1, p2, irp);
//                         } else {
//                             children = crossover_time_based(p1, p2, irp);
//                         }
//                     } else {
//                         children = {p1, p2};
//                     }
            
//             // Processa Filho 1
//             crossover_tries++;
//             // (Sem mutação, conforme solicitado)
//             build_routes_for_individual(children.first, irp, aco_params);
//             check_feasibility(children.first, irp); 
            
//             if (children.first.is_feasible) {
//                 calculate_total_cost(children.first, irp); 
                
//                 // --- MUDANÇA 3: Probabilidade na Busca Local ---
//                 // Aplica apenas 50% das vezes para economizar tempo e manter diversidade
//                 if (randreal() < 0.50) { 
//                     run_simple_reinsertion_search(children.first, irp, aco_params, false);
//                 }
//                 newPop.push_back(children.first);
//             }
            
//             if (newPop.size() >= children_target_count) break;

//             // Processa Filho 2
//             crossover_tries++;
//             build_routes_for_individual(children.second, irp, aco_params);
//             check_feasibility(children.second, irp);
//             if (children.second.is_feasible) {
//                 calculate_total_cost(children.second, irp); 
//                 if (randreal() < 0.50) {
//                     run_simple_reinsertion_search(children.second, irp, aco_params, false);
//                 }
//                 newPop.push_back(children.second);
//             }
//         }

//         // 3. Imigrantes (Diversidade Pura)
//         while (newPop.size() < ga_params.popSize) {
//             Individual ind = make_new_heuristic_individual(irp, aco_params);
//             build_routes_for_individual(ind, irp, aco_params);
//             check_feasibility(ind, irp);
//             if (ind.is_feasible) {
//                 calculate_total_cost(ind, irp);
//                 // (Sem busca local nos imigrantes para serem "sangue novo")
//             } else {
//                 ind.fitness = LARGE_COST;
//             }
//             newPop.push_back(ind);
//         }
        
//         pop.swap(newPop);

//         // 4. Busca Pesada (Shaking) Periódica
//         if (gen > 0 && gen % 100 == 0) {
//             if (verbose) std::cout << "--- Gen " << gen << ": Executando Busca Pesada (VND-DS) na Elite ---\n";
//             std::sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b){
//                 return a.fitness < b.fitness;
//             });
//             // Aplica apenas nos top 5 para refinar a elite
//             for(int i=0; i<std::min((int)pop.size(), 5); ++i) {
//                 if (pop[i].is_feasible) run_vnd_ds_operator(pop[i], irp, aco_params, false);
//             }
//         }

//         // Atualiza Melhor Global
//         std::sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b){
//             return a.fitness < b.fitness;
//         });
//         if (pop.front().fitness < bestOverall.fitness) {
//             bestOverall = pop.front();
//         }
        
//         // Log
//         if (verbose || (gen % 10 == 0) || (gen == ga_params.nGen - 1)) {
//             double avg_fit = 0.0;
//             int count = 0;
//             for(const auto& ind : pop) if(ind.is_feasible) { avg_fit += ind.fitness; count++; }
//             if(count > 0) avg_fit /= count;

//             std::cout << "Gen " << std::setw(4) << gen + 1 
//                       << " | Melhor Pop: " << std::fixed << std::setprecision(2) << pop.front().fitness
//                       << " | Média: " << avg_fit
//                       << " | Global: " << bestOverall.fitness << "\n";
//         }
//     }
    
//     // Refinamento Final
//     if (bestOverall.is_feasible) {
//         std::cout << "\nOtimização final na melhor solução...\n";
//         run_vnd_ds_operator(bestOverall, irp, aco_params, true); // Verbose na final
//     }

//     return bestOverall;
// }



// ... (Includes e funções anteriores) ...

Individual run_genetic_algorithm(const IRP& irp, const GA_Params& ga_params, const ACO_Params& aco_params, bool verbose) {
    
    const double LARGE_COST = 1e18;
    vector<Individual> pop;
    pop.reserve(ga_params.popSize); 
    
    // --- 1. INICIALIZAÇÃO (COM DIVERSIDADE) ---
    std::cout << "Inicializando população..." << std::endl;
    int initial_pop_size = ga_params.popSize; 
    
    for (int i = 0; i < initial_pop_size; ++i) {
        Individual ind = make_new_heuristic_individual(irp, aco_params);
        
        build_routes_for_individual(ind, irp, aco_params);
        check_feasibility(ind, irp); 
        
        if (ind.is_feasible) {
            calculate_total_cost(ind, irp);
        } else {
            ind.fitness = LARGE_COST;
        }
        
        // --- CORREÇÃO DO LOOP INFINITO ---
        // Adiciona o indivíduo factível (com custo real) ou
        // infactível (com custo LARGE_COST).
        // A lógica 'i--' foi removida.
        pop.push_back(ind);
        // --- FIM DA CORREÇÃO ---

        if (i % 10 == 0 && i > 0) std::cout << "." << std::flush;
    }
    std::cout << "\nPopulação inicial gerada." << std::endl;

    // --- CORREÇÃO NA SELEÇÃO DO MELHOR INICIAL ---
    // Ordena a população para garantir que 'bestOverall' seja
    // o melhor indivíduo factível, e não apenas pop[0].
    std::sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b){
        return a.fitness < b.fitness;
    });
    
    Individual bestOverall = pop.front();
    if (!bestOverall.is_feasible) {
        std::cerr << "AVISO: Nenhum indivíduo factível gerado na população inicial." << std::endl;
        // bestOverall já terá fitness = LARGE_COST, o que é correto.
    }
    // --- FIM DA CORREÇÃO ---

    int num_elites = (int)(ga_params.popSize * 0.10);

    // --- 2. LOOP EVOLUCIONÁRIO (HGS) ---
    std::cout << "\nIniciando loop do GA (HGS) para " << ga_params.nGen << " iterações..." << std::endl;
    
    for (int gen = 0; gen < ga_params.nGen; ++gen) {
        
        // (Log movido para o final do loop para refletir o estado da geração)

        vector<Individual> newPop;
        newPop.reserve(ga_params.popSize);

        // 1. Elitismo
        int num_elites = (int)(ga_params.popSize * 0.10);
        for(int i = 0; i < num_elites && i < pop.size(); ++i) {
            newPop.push_back(pop[i]); // pop já está ordenada
        }

        // 2. Geração de Filhos
        int children_target_count = num_elites + (int)(ga_params.popSize * 0.70);
        int crossover_tries = 0; // Adicionado contador para segurança

        while (newPop.size() < children_target_count && crossover_tries < ga_params.crossover_max_tries) {
            crossover_tries++; // Incrementa tentativas

            Individual parent1 = tournamentSelect(pop, ga_params.tournamentK);
            Individual parent2 = tournamentSelect(pop, ga_params.tournamentK);
            
            std::pair<Individual, Individual> children;
            if (randreal() < ga_params.pCrossover) {
                 double r = randreal();
                 if (r < 0.33) {
                    children = one_point_crossover_customer(parent1, parent2, irp);
                 } else if (r < 0.66) {
                    children = two_point_crossover_customer(parent1, parent2, irp);
                 } else {
                    children = crossover_time_based(parent1, parent2, irp);
                 }
            } else {
                children = {parent1, parent2};
            }
            
            // --- Processa o Filho 1 ---
            build_routes_for_individual(children.first, irp, aco_params);
            
            // Aplica busca local "pesada" (VND-DS)
            if(gen > 0 && gen % 100 == 0) {
                // (Renomeei 'busca_local' para a função que definimos)
                run_vnd_ds_operator(children.first, irp, aco_params, false);
            } else if (randreal() > 0.5f) {
                // Aplica busca local "leve"
                run_simple_reinsertion_search(children.first, irp, aco_params, false);
            }
            
            check_feasibility(children.first, irp); 
            
            if (children.first.is_feasible) {
                calculate_total_cost(children.first, irp); 
            } else {
                children.first.fitness = LARGE_COST;
            }
            newPop.push_back(children.first);
            
            if (newPop.size() >= children_target_count) break;

            // --- Processa o Filho 2 ---
            build_routes_for_individual(children.second, irp, aco_params);
            
            if(gen > 0 && gen % 100 == 0) {
                run_vnd_ds_operator(children.second, irp, aco_params, false);
            } else {
                run_simple_reinsertion_search(children.second, irp, aco_params, false);
            }

            check_feasibility(children.second, irp);
            
            if (children.second.is_feasible) {
                calculate_total_cost(children.second, irp); 
            } else {
                children.second.fitness = LARGE_COST;
            }
            newPop.push_back(children.second);
        } // Fim do while (filhos)

        // 3. Geração de novas soluções (Imigrantes)
        while (newPop.size() < ga_params.popSize) {
            Individual ind = make_new_heuristic_individual(irp, aco_params);
            build_routes_for_individual(ind, irp, aco_params);

            // Aplica busca local "leve" no imigrante
            run_simple_reinsertion_search(ind, irp, aco_params, false);

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
        
        if (verbose || (gen % 10 == 0) || (gen == ga_params.nGen - 1)) {
            std::cout << "Gen " << std::setw(4) << gen + 1 << "/" << ga_params.nGen
                      << " | Fact.: " << std::setw(3) << feasible_count << "/" << (int)pop.size()
                      << " | Melhor: " << std::fixed << std::setprecision(2) << min_fit
                      << " | Média(fact): " << avg_fit
                      << " | Global: " << bestOverall.fitness << std::endl;
        }
    }
    
    // Aplica a busca local "pesada" na melhor solução final
    if(bestOverall.is_feasible) {
        std::cout << "\nOtimização final (Pesada) na melhor solução...\n";
        run_vnd_ds_operator(bestOverall, irp, aco_params, verbose);
    }
    
    return bestOverall;
}