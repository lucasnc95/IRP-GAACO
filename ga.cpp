#include "ga.hpp"
#include "utils.hpp"
#include "evaluation.hpp"
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <utility>

using std::vector;

// Individual make_simple_random_individual(const IRP& irp) {
//     Individual ind(irp.nPeriods, irp.nCustomers);
//     vector<long> current_inv(irp.nCustomers);
//     for(int i = 0; i < irp.nCustomers; ++i) current_inv[i] = irp.customers[i].initialInv;

//     float dynamic_probability = randreal();
//     double depot_inv_cost = irp.depots[0].invCost;

//     for (int t = 0; t < irp.nPeriods; ++t) {
//         // --- NOVA LÓGICA: Representação explícita dos veículos e suas cargas ---
//         vector<long> vehicle_loads(irp.nVehicles, 0);
//         // Mapeia qual cliente está em qual veículo (-1 = não alocado)
//         vector<int> customer_to_vehicle_map(irp.nCustomers, -1);
//         vector<bool> is_mandatory(irp.nCustomers, false);

//         // --- PASSO 1: ATENDER ENTREGAS MÍNIMAS OBRIGATÓRIAS ---
//         for (int c = 0; c < irp.nCustomers; ++c) {
//             long inv_after_demand = current_inv[c] - irp.customers[c].demand[t];
//             if (inv_after_demand < irp.customers[c].minLevelInv) {
//                 is_mandatory[c] = true;
//                 long needed = irp.customers[c].minLevelInv - inv_after_demand;
//                 long delivery_amount = std::min(needed, irp.customers[c].maxLevelInv - current_inv[c]);
//                 if (delivery_amount <= 0) continue;

//                 // Encontra um veículo que possa fazer esta entrega
//                 for (int v = 0; v < irp.nVehicles; ++v) {
//                     if (vehicle_loads[v] + delivery_amount <= irp.Capacity) {
//                         ind.deliveries[t][c] = delivery_amount;
//                         vehicle_loads[v] += delivery_amount;
//                         customer_to_vehicle_map[c] = v; // Associa cliente ao veículo
//                         break; // Veículo encontrado, passa para o próximo cliente
//                     }
//                 }
//             }
//         }

//         // --- PASSO 2: ENTREGAS OPORTUNISTAS E DE "TOP-UP" ---
//         vector<int> customer_order(irp.nCustomers);
//         std::iota(customer_order.begin(), customer_order.end(), 0);
//         std::shuffle(customer_order.begin(), customer_order.end(), rng);

//         for (int c : customer_order) {
//             bool was_mandatory = is_mandatory[c];
//             if (!was_mandatory && randreal() >= dynamic_probability) {
//                 continue;
//             }

//             long inv_before_topup = current_inv[c];
//             long current_delivery = ind.deliveries[t][c];

//             // Tenta encontrar um veículo para este cliente
//             int vehicle_idx = customer_to_vehicle_map[c];
//             if (vehicle_idx == -1) { // Cliente opcional ainda não alocado
//                 for (int v = 0; v < irp.nVehicles; ++v) {
//                     if (vehicle_loads[v] < irp.Capacity) {
//                         vehicle_idx = v;
//                         break;
//                     }
//                 }
//             }
//             if (vehicle_idx == -1) continue; // Nenhum veículo disponível

//             // Calcula o máximo que PODE ser adicionado
//             long max_extra_q = std::min({
//                 (irp.customers[c].maxLevelInv + irp.customers[c].demand[t]) - (inv_before_topup + current_delivery),
//                 (long)irp.Capacity - vehicle_loads[vehicle_idx]
//             });

//             if (max_extra_q > 0) {
//                 int q_extra = 0;
//                 if (irp.customers[c].invCost < depot_inv_cost) {
//                     q_extra = max_extra_q;
//                 } else {
//                     q_extra = was_mandatory ? randint(0, max_extra_q) : randint(1, max_extra_q);
//                 }

//                 if (q_extra > 0) {
//                     ind.deliveries[t][c] += q_extra;
//                     vehicle_loads[vehicle_idx] += q_extra;
//                     customer_to_vehicle_map[c] = vehicle_idx; // Confirma associação
//                 }
//             }
//         }

//         // Atualiza o inventário para o próximo período com as entregas finais
//         for (int c = 0; c < irp.nCustomers; ++c) {
//             current_inv[c] += (long)ind.deliveries[t][c] - irp.customers[c].demand[t];
//         }
//     }
//     return ind;
// }



// Individual make_simple_random_individual_2(const IRP& irp) {
//     Individual ind(irp.nPeriods, irp.nCustomers);
//     vector<long> current_inv(irp.nCustomers);
//     for(int i = 0; i < irp.nCustomers; ++i) current_inv[i] = irp.customers[i].initialInv;
//     float probability = randreal();
//     // Calcula a necessidade líquida total para cada cliente, para a regra de custo
//     vector<long> net_demand_to_satisfy(irp.nCustomers, 0);
//     vector<long> total_delivered_so_far(irp.nCustomers, 0);
//     for (int c = 0; c < irp.nCustomers; ++c) {
//         long total_demand = 0;
//         for (int demand_val : irp.customers[c].demand) total_demand += demand_val;
//         net_demand_to_satisfy[c] = std::max(0L, total_demand - irp.customers[c].initialInv);
//     }

//     long fleet_capacity = (long)irp.nVehicles * irp.Capacity;
//     double depot_inv_cost = irp.depots[0].invCost;

//     for (int t = 0; t < irp.nPeriods; ++t) {
//         long period_load = 0;
//         vector<int> mandatory_custs;
        
//         // --- PASSO 1: IDENTIFICAR E ATENDER ENTREGAS MÍNIMAS OBRIGATÓRIAS ---
//         for (int c = 0; c < irp.nCustomers; ++c) {
//             if (current_inv[c] < irp.customers[c].demand[t]) {
//                 mandatory_custs.push_back(c);
//                 long needed = irp.customers[c].demand[t] - current_inv[c];
//                 ind.deliveries[t][c] = needed;
//                 period_load += needed;
//             }
//         }
        
//         // Se a carga obrigatória já excede a capacidade, não há mais nada a fazer neste período
//         if (period_load > fleet_capacity) {
//              // (Opcional: implementar uma lógica de priorização aqui se isso ocorrer)
//         } else {
//             // --- PASSO 2: ATENDER ENTREGAS ANTECIPADAS/OPORTUNISTAS ---
//             long remaining_fleet_capacity = fleet_capacity - period_load;
            
//             vector<int> customer_order(irp.nCustomers);
//             std::iota(customer_order.begin(), customer_order.end(), 0);
//             std::shuffle(customer_order.begin(), customer_order.end(), rng);

//             for (int c : customer_order) {
//                 if (remaining_fleet_capacity <= 0) break;

//                 // Probabilidade de tentar uma entrega oportunista
//                 if (randreal() < probability) { 
//                     long inv_before_delivery = current_inv[c] + ind.deliveries[t][c];
                    
//                     // Sorteia um horizonte de antecipação N
//                     int N = randint(0, irp.nPeriods - t - 1);
//                     long look_ahead_demand = 0;
//                     for (int p = t; p <= t + N; ++p) {
//                         look_ahead_demand += irp.customers[c].demand[p];
//                     }
                    
//                     long amount_needed = look_ahead_demand - current_inv[c];
//                     if (amount_needed <= 0) continue;

//                     // Calcula a entrega máxima possível
//                     long max_q = std::min({
//                         (long)irp.Capacity - ind.deliveries[t][c], // O que ainda cabe no veículo para este cliente
//                         remaining_fleet_capacity,                  // O que ainda cabe na frota
//                         irp.customers[c].maxLevelInv - inv_before_delivery // Espaço no armazém
//                     });

//                     // Aplica a regra de custo
//                     if (irp.customers[c].invCost >= depot_inv_cost) {
//                         long remaining_net_demand = net_demand_to_satisfy[c] - total_delivered_so_far[c] - ind.deliveries[t][c];
//                         max_q = std::min(max_q, remaining_net_demand);
//                     }

//                     if (max_q > 0) {
//                         int q_extra = randint(1, max_q);
//                         ind.deliveries[t][c] += q_extra;
//                         remaining_fleet_capacity -= q_extra;
//                     }
//                 }
//             }
//         }

//         // Atualiza o inventário e o total entregue para o próximo período
//         for (int c = 0; c < irp.nCustomers; ++c) {
//             current_inv[c] += (long)ind.deliveries[t][c] - irp.customers[c].demand[t];
//             total_delivered_so_far[c] += ind.deliveries[t][c];
//         }
//     }
//     return ind;
// }

Individual make_simple_random_individual(const IRP& irp) {
    Individual ind(irp.nPeriods, irp.nCustomers);
    vector<long> current_inv(irp.nCustomers);
    for(int i = 0; i < irp.nCustomers; ++i) current_inv[i] = irp.customers[i].initialInv;

    double depot_inv_cost = irp.depots[0].invCost;
    float dynamic_probability = randreal();

    for (int t = 0; t < irp.nPeriods; ++t) {
        vector<long> vehicle_loads(irp.nVehicles, 0);
        vector<int> customer_to_vehicle_map(irp.nCustomers, -1);
        vector<bool> is_mandatory(irp.nCustomers, false);

        // --- PASSO 1: ATENDER ENTREGAS MÍNIMAS OBRIGATÓRIAS ---
        for (int c = 0; c < irp.nCustomers; ++c) {
            long inv_after_demand = current_inv[c] - irp.customers[c].demand[t];
            if (inv_after_demand < irp.customers[c].minLevelInv) {
                is_mandatory[c] = true;
                long needed = irp.customers[c].minLevelInv - inv_after_demand;
                
                // Limita a entrega pelo espaço de "pico" e capacidade do veículo
                long space_available = irp.customers[c].maxLevelInv - current_inv[c];
                long delivery_amount = std::min({needed, space_available, (long)irp.Capacity});

                if (delivery_amount <= 0) continue;

                for (int v = 0; v < irp.nVehicles; ++v) {
                    if (vehicle_loads[v] + delivery_amount <= irp.Capacity) {
                        ind.deliveries[t][c] = delivery_amount;
                        vehicle_loads[v] += delivery_amount;
                        customer_to_vehicle_map[c] = v;
                        break;
                    }
                }
            }
        }

        // --- PASSO 2: ENTREGAS OPORTUNISTAS E DE "TOP-UP" ---
        vector<int> customer_order(irp.nCustomers);
        std::iota(customer_order.begin(), customer_order.end(), 0);
        std::shuffle(customer_order.begin(), customer_order.end(), rng);

        for (int c : customer_order) {
            bool was_mandatory = is_mandatory[c];
            if (!was_mandatory && randreal() >= dynamic_probability) {
                continue;
            }

            long inv_before_topup = current_inv[c];
            long current_delivery = ind.deliveries[t][c];

            int vehicle_idx = customer_to_vehicle_map[c];
            if (vehicle_idx == -1) {
                for (int v = 0; v < irp.nVehicles; ++v) {
                    if (vehicle_loads[v] < irp.Capacity) {
                        vehicle_idx = v;
                        break;
                    }
                }
            }
            if (vehicle_idx == -1) continue;

            // Calcula o máximo que PODE ser adicionado, respeitando o "pico"
            long space = irp.customers[c].maxLevelInv - (inv_before_topup + current_delivery);

            long max_extra_q = std::min({
                space,
                (long)irp.Capacity - vehicle_loads[vehicle_idx]
            });

            if (max_extra_q > 0) {
                int q_extra = 0;
                if (irp.customers[c].invCost < depot_inv_cost) {
                    q_extra = max_extra_q;
                } else {
                    q_extra = was_mandatory ? randint(0, max_extra_q) : randint(1, max_extra_q);
                }

                if (q_extra > 0) {
                    ind.deliveries[t][c] += q_extra;
                    vehicle_loads[vehicle_idx] += q_extra;
                    customer_to_vehicle_map[c] = vehicle_idx;
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


std::pair<Individual, Individual> one_point_crossover_customer(const Individual& a, const Individual& b, const IRP& irp) {
    Individual child1(irp.nPeriods, irp.nCustomers);
    Individual child2(irp.nPeriods, irp.nCustomers);
    
    int crosspoint = randint(0, irp.nCustomers);

    for (int c = 0; c < irp.nCustomers; ++c) {
        // Se c está antes do ponto de corte
        if (c < crosspoint) {
            for (int t = 0; t < irp.nPeriods; ++t) {
                child1.deliveries[t][c] = a.deliveries[t][c]; // Filho 1 herda do Pai A
                child2.deliveries[t][c] = b.deliveries[t][c]; // Filho 2 herda do Pai B
            }
        } 
        // Se c está depois do ponto de corte
        else {
            for (int t = 0; t < irp.nPeriods; ++t) {
                child1.deliveries[t][c] = b.deliveries[t][c]; // Filho 1 herda do Pai B
                child2.deliveries[t][c] = a.deliveries[t][c]; // Filho 2 herda do Pai A
            }
        }
    }
    return {child1, child2};
}

// <-- MUDANÇA: Crossover de 2 pontos agora gera dois filhos complementares -->
std::pair<Individual, Individual> two_point_crossover_customer(const Individual& a, const Individual& b, const IRP& irp) {
    Individual child1(irp.nPeriods, irp.nCustomers);
    Individual child2(irp.nPeriods, irp.nCustomers);

    int p1 = randint(0, irp.nCustomers);
    int p2 = randint(0, irp.nCustomers);
    if (p1 > p2) std::swap(p1, p2);

    for (int c = 0; c < irp.nCustomers; ++c) {
        // Se c está dentro do segmento de troca
        if (c >= p1 && c < p2) {
            for (int t = 0; t < irp.nPeriods; ++t) {
                child1.deliveries[t][c] = b.deliveries[t][c]; // Filho 1 herda do Pai B
                child2.deliveries[t][c] = a.deliveries[t][c]; // Filho 2 herda do Pai A
            }
        } 
        // Se c está fora do segmento
        else {
            for (int t = 0; t < irp.nPeriods; ++t) {
                child1.deliveries[t][c] = a.deliveries[t][c]; // Filho 1 herda do Pai A
                child2.deliveries[t][c] = b.deliveries[t][c]; // Filho 2 herda do Pai B
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

// Individual run_genetic_algorithm(const IRP& irp, const GA_Params& ga_params, const ACO_Params& aco_params, bool verbose) {
//     vector<Individual> pop;
//     pop.reserve(ga_params.popSize);
    
//     std::cout << "Inicializando população..." << std::endl;
//     while(pop.size() < ga_params.popSize/2) {
//         Individual ind = make_simple_random_individual(irp);
//         if (check_feasibility(ind, irp)) {
//             pop.push_back(ind);
//         }
//     }
//     while(pop.size() < ga_params.popSize) {
//         Individual ind = make_simple_random_individual_2(irp);
//         if (check_feasibility(ind, irp)) {
//             pop.push_back(ind);
//         }
//     }
//     std::cout << "Avaliando população inicial..." << std::endl;
//     for (auto& ind : pop) {
//         evaluate_and_fill(ind, irp, aco_params);
//     }

//     Individual bestOverall = *std::min_element(pop.begin(), pop.end(), 
//         [](const Individual& a, const Individual& b){ return a.fitness < b.fitness; });
    
//     int stagnation_counter = 0;
//     double last_best_fitness = bestOverall.fitness;
    
//     std::cout << "\nIniciando loop do GA para " << ga_params.nGen << " gerações..." << std::endl;
//     for (int gen = 0; gen < ga_params.nGen; ++gen) {
//         vector<Individual> newPop;
//         newPop.reserve(ga_params.popSize);

//         while (newPop.size() < ga_params.popSize) {
//             Individual parent1 = tournamentSelect(pop, ga_params.tournamentK);
//             Individual parent2 = tournamentSelect(pop, ga_params.tournamentK);
            
//             std::pair<Individual, Individual> children;

//             if (randreal() < ga_params.pCrossover) {
//                 if (randreal() < 0.5) {
//                     children = one_point_crossover_customer(parent1, parent2, irp);
//                 } else {
//                     children = two_point_crossover_customer(parent1, parent2, irp);
//                 }
//             } else {
//                 children = {parent1, parent2};
//             }
            
//             // --- MUDANÇA: A chamada de mutação agora é incondicional ---
//             // A probabilidade é verificada dentro da própria função
//             advance_portion_mutation(children.first, irp, ga_params.pMutation);
//             advance_portion_mutation(children.second, irp, ga_params.pMutation);

//             // Processa o primeiro filho
//             if (check_feasibility(children.first, irp)) {
//                 evaluate_and_fill(children.first, irp, aco_params);
//                 newPop.push_back(children.first);
//             }

//             // Processa o segundo filho, se ainda houver espaço
//             if (newPop.size() < ga_params.popSize) {
//                 if (check_feasibility(children.second, irp)) {
//                     evaluate_and_fill(children.second, irp, aco_params);
//                     newPop.push_back(children.second);
//                 }
//             }
//         }

//         // Elitismo
//         Individual& bestPrev = *std::min_element(pop.begin(), pop.end(), 
//             [](const Individual& a, const Individual& b){ return a.fitness < b.fitness; });
//         auto worst_it = std::max_element(newPop.begin(), newPop.end(),
//             [](const Individual& a, const Individual& b){ return a.fitness < b.fitness; });
//         if(bestPrev.fitness < worst_it->fitness) *worst_it = bestPrev;
//         pop.swap(newPop);

//         Individual& genBest = *std::min_element(pop.begin(), pop.end(), 
//             [](const Individual& a, const Individual& b){ return a.fitness < b.fitness; });
//         if (genBest.fitness < bestOverall.fitness) bestOverall = genBest;
        
//         // Estagnação e Reinicialização
//         if (bestOverall.fitness < last_best_fitness - 1e-9) {
//             last_best_fitness = bestOverall.fitness;
//             stagnation_counter = 0;
//         } else {
//             stagnation_counter++;
//         }
//         if (ga_params.stagnation_threshold > 0 && stagnation_counter >= ga_params.stagnation_threshold) {
//              if (ga_params.stagnation_threshold > 0 && stagnation_counter >= ga_params.stagnation_threshold) {
//             std::cout << "\n*** ESTAGNAÇÃO DETECTADA APÓS " << stagnation_counter << " GERAÇÕES. REINICIALIZANDO POPULAÇÃO... ***\n" << std::endl;
            
//             vector<Individual> restarted_pop;
//             restarted_pop.reserve(ga_params.popSize);
            
//             // 1. Mantém o melhor indivíduo de todos
//             restarted_pop.push_back(bestOverall);
            
//             // 2. Preenche o resto com novos indivíduos aleatórios e factíveis
//             while(restarted_pop.size() < ga_params.popSize/2) {
//                 Individual ind = make_simple_random_individual(irp);
//                 if (check_feasibility(ind, irp)) {
//                     restarted_pop.push_back(ind);
//                 }
//             }
//             while(restarted_pop.size() < ga_params.popSize) {
//                 Individual ind = make_simple_random_individual_2(irp);
//                 if (check_feasibility(ind, irp)) {
//                     restarted_pop.push_back(ind);
//                 }
//             }

//             // 3. Avalia apenas os novos indivíduos (o melhor já está avaliado)
//             for(size_t i = 1; i < restarted_pop.size(); ++i) {
//                  evaluate_and_fill(restarted_pop[i], irp, aco_params);
//             }
            
//             pop = restarted_pop; // Substitui a população estagnada
//             stagnation_counter = 0; // Reseta o contador
//         }
//         }
        
//         // Log da Geração
//         double min_fit = (*std::min_element(pop.begin(), pop.end(), [](const auto& a, const auto& b){ return a.fitness < b.fitness; })).fitness;
//         double max_fit = (*std::max_element(pop.begin(), pop.end(), [](const auto& a, const auto& b){ return a.fitness < b.fitness; })).fitness;
//         double avg_fit = std::accumulate(pop.begin(), pop.end(), 0.0, [](double sum, const auto& ind){ return sum + ind.fitness; }) / pop.size();
        
//         std::cout << "Gen " << std::setw(4) << gen + 1 << "/" << ga_params.nGen
//                   << " | Melhor: " << std::fixed << std::setprecision(2) << min_fit
//                   << " | Média: " << avg_fit
//                   << " | Pior: " << max_fit
//                   << " | Global: " << bestOverall.fitness << "\n";
//     }
//     return bestOverall;
// }







// Individual run_genetic_algorithm(const IRP& irp, const GA_Params& ga_params, const ACO_Params& aco_params, bool verbose) {
//     vector<Individual> pop;
//     pop.reserve(ga_params.popSize);
    
//     std::cout << "Inicializando população..." << std::endl;
//     while(pop.size() < ga_params.popSize/2) {
//         Individual ind = make_simple_random_individual(irp);
//         if (check_feasibility(ind, irp)) {
//             pop.push_back(ind);
//         }
//     }
//     while(pop.size() < ga_params.popSize) {
//         Individual ind = make_simple_random_individual_2(irp);
//         if (check_feasibility(ind, irp)) {
//             pop.push_back(ind);
//         }
//     }
//     std::cout << "Avaliando população inicial..." << std::endl;
//     for (auto& ind : pop) {
//         evaluate_and_fill(ind, irp, aco_params);
//     }

//     Individual bestOverall = *std::min_element(pop.begin(), pop.end(), 
//         [](const Individual& a, const Individual& b){ return a.fitness < b.fitness; });
    
//     int stagnation_counter = 0;
//     double last_best_fitness = bestOverall.fitness;
    
//     std::cout << "\nIniciando loop do GA para " << ga_params.nGen << " gerações..." << std::endl;
//     for (int gen = 0; gen < ga_params.nGen; ++gen) {
//         vector<Individual> newPop;
//         newPop.reserve(ga_params.popSize);

//         while (newPop.size() < ga_params.popSize) {
//             Individual parent1 = tournamentSelect(pop, ga_params.tournamentK);
//             Individual parent2 = tournamentSelect(pop, ga_params.tournamentK);
            
//             std::pair<Individual, Individual> children;

//             if (randreal() < ga_params.pCrossover) {
//                 if (randreal() < 0.5) {
//                     children = one_point_crossover_customer(parent1, parent2, irp);
//                 } else {
//                     children = two_point_crossover_customer(parent1, parent2, irp);
//                 }
//             } else {
//                 children = {parent1, parent2};
//             }
            
//             advance_portion_mutation(children.first, irp, ga_params.pMutation);
//             advance_portion_mutation(children.second, irp, ga_params.pMutation);

//             if (check_feasibility(children.first, irp)) {
//                 evaluate_and_fill(children.first, irp, aco_params);
//                 newPop.push_back(children.first);
//             }

//             if (newPop.size() < ga_params.popSize) {
//                 if (check_feasibility(children.second, irp)) {
//                     evaluate_and_fill(children.second, irp, aco_params);
//                     newPop.push_back(children.second);
//                 }
//             }
//         }

//         // --- MUDANÇA: Estratégia de Seleção (Pais + Filhos) ---
//         // 1. Junta a população de pais (pop) com a de filhos (newPop)
//         vector<Individual> combined_pop = pop;
//         combined_pop.insert(combined_pop.end(), newPop.begin(), newPop.end());

//         // 2. Ordena a população combinada (tamanho 2 * popSize) pelo fitness
//         std::sort(combined_pop.begin(), combined_pop.end(), [](const Individual& a, const Individual& b) {
//             return a.fitness < b.fitness;
//         });
        
//         // 3. Os 'popSize' melhores sobrevivem para a próxima geração
//         pop.assign(combined_pop.begin(), combined_pop.begin() + ga_params.popSize);
//         // --- FIM DA MUDANÇA ---


//         // Atualiza o melhor global (que agora é sempre o primeiro indivíduo)
//         if (pop[0].fitness < bestOverall.fitness) {
//             bestOverall = pop[0];
//         }
        
//         // Estagnação e Reinicialização
//         if (bestOverall.fitness < last_best_fitness - 1e-9) {
//             last_best_fitness = bestOverall.fitness;
//             stagnation_counter = 0;
//         } else {
//             stagnation_counter++;
//         }
//         if (ga_params.stagnation_threshold > 0 && stagnation_counter >= ga_params.stagnation_threshold) {
//             std::cout << "\n*** ESTAGNAÇÃO DETECTADA APÓS " << stagnation_counter << " GERAÇÕES. REINICIALIZANDO POPULAÇÃO... ***\n" << std::endl;
            
//             vector<Individual> restarted_pop;
//             restarted_pop.reserve(ga_params.popSize);
//             restarted_pop.push_back(bestOverall);
            
//             while(restarted_pop.size() < ga_params.popSize/2) {
//                 Individual ind = make_simple_random_individual(irp);
//                 if (check_feasibility(ind, irp)) {
//                     restarted_pop.push_back(ind);
//                 }
//             }
//             while(restarted_pop.size() < ga_params.popSize) {
//                 Individual ind = make_simple_random_individual_2(irp);
//                 if (check_feasibility(ind, irp)) {
//                     restarted_pop.push_back(ind);
//                 }
//             }

//             for(size_t i = 1; i < restarted_pop.size(); ++i) {
//                  evaluate_and_fill(restarted_pop[i], irp, aco_params);
//             }
            
//             pop = restarted_pop;
//             stagnation_counter = 0;
//         }
        
//         // Log da Geração (agora mais eficiente, pois a população está ordenada)
//         double min_fit = pop[0].fitness;
//         double max_fit = pop.back().fitness;
//         double avg_fit = std::accumulate(pop.begin(), pop.end(), 0.0, [](double sum, const auto& ind){ return sum + ind.fitness; }) / pop.size();
        
//         std::cout << "Gen " << std::setw(4) << gen + 1 << "/" << ga_params.nGen
//                   << " | Melhor: " << std::fixed << std::setprecision(2) << min_fit
//                   << " | Média: " << avg_fit
//                   << " | Pior: " << max_fit
//                   << " | Global: " << bestOverall.fitness << "\n";
//     }
//     return bestOverall;
// }


Individual run_genetic_algorithm(const IRP& irp, const GA_Params& ga_params, const ACO_Params& aco_params, bool verbose) {
    vector<Individual> pop;
    pop.reserve(ga_params.popSize);
    
    std::cout << "Inicializando população..." << std::endl;
    while(pop.size() < ga_params.popSize) {
        Individual ind = make_simple_random_individual(irp);
        if (check_feasibility(ind, irp)) {
            pop.push_back(ind);
        }
    }

    std::cout << "Avaliando população inicial..." << std::endl;
    for (auto& ind : pop) {
        evaluate_and_fill(ind, irp, aco_params);
    }

    Individual bestOverall = *std::min_element(pop.begin(), pop.end(), 
        [](const Individual& a, const Individual& b){ return a.fitness < b.fitness; });
    
    std::cout << "\nIniciando loop do GA para " << ga_params.nGen << " gerações..." << std::endl;
    for (int gen = 0; gen < ga_params.nGen; ++gen) {
        vector<Individual> newPop;
        newPop.reserve(ga_params.popSize);

        // --- ESTRATÉGIA DE SELEÇÃO 10% / 70% / 20% ---

        // 1. (10%) Elitismo
        std::sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b){
            return a.fitness < b.fitness;
        });
        
        int num_elites = (int)(ga_params.popSize * 0.10);
        for(int i = 0; i < num_elites && i < pop.size(); ++i) {
            newPop.push_back(pop[i]);
        }

        // 2. (70%) Geração de Filhos (Crossover e Mutação)
        int children_target_count = num_elites + (int)(ga_params.popSize * 0.70);
        while (newPop.size() < children_target_count) {
            Individual parent1 = tournamentSelect(pop, ga_params.tournamentK);
            Individual parent2 = tournamentSelect(pop, ga_params.tournamentK);
            
            std::pair<Individual, Individual> children;

            if (randreal() < ga_params.pCrossover) {
                if (randreal() < 0.5) {
                    children = one_point_crossover_customer(parent1, parent2, irp);
                } else {
                    children = two_point_crossover_customer(parent1, parent2, irp);
                }
            } else {
                children = {parent1, parent2};
            }
            
            advance_portion_mutation(children.first, irp, ga_params.pMutation);
            
            if (check_feasibility(children.first, irp)) {
                evaluate_and_fill(children.first, irp, aco_params);
                newPop.push_back(children.first);
            }

            if (newPop.size() < children_target_count) {
                advance_portion_mutation(children.second, irp, ga_params.pMutation);
                if (check_feasibility(children.second, irp)) {
                    evaluate_and_fill(children.second, irp, aco_params);
                    newPop.push_back(children.second);
                }
            }
        }

        // 3. (20%) Geração de Imigrantes (Novos Aleatórios)
        while (newPop.size() < ga_params.popSize) {
            Individual ind = make_simple_random_individual(irp);
            if (check_feasibility(ind, irp)) {
                evaluate_and_fill(ind, irp, aco_params);
                newPop.push_back(ind);
            }
        }
        
        pop.swap(newPop);

        Individual& genBest = *std::min_element(pop.begin(), pop.end(), 
            [](const Individual& a, const Individual& b){ return a.fitness < b.fitness; });
        if (genBest.fitness < bestOverall.fitness) {
            bestOverall = genBest;
        }
        
        // Log da Geração
        double min_fit = (*std::min_element(pop.begin(), pop.end(), [](const auto& a, const auto& b){ return a.fitness < b.fitness; })).fitness;
        double max_fit = (*std::max_element(pop.begin(), pop.end(), [](const auto& a, const auto& b){ return a.fitness < b.fitness; })).fitness;
        double avg_fit = std::accumulate(pop.begin(), pop.end(), 0.0, [](double sum, const auto& ind){ return sum + ind.fitness; }) / pop.size();
        
        std::cout << "Gen " << std::setw(4) << gen + 1 << "/" << ga_params.nGen
                  << " | Melhor: " << std::fixed << std::setprecision(2) << min_fit
                  << " | Média: " << avg_fit
                  << " | Pior: " << max_fit
                  << " | Global: " << bestOverall.fitness << "\n";
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

// void advance_portion_mutation(Individual& ind, const IRP& irp, double pMutation) {
//     // Percorre todos os clientes
//     for (int c = 0; c < irp.nCustomers; ++c) {
//         // Percorre os períodos a partir do segundo
//         for (int t_source = 1; t_source < irp.nPeriods; ++t_source) {
//             // Verifica a probabilidade de aplicar a mutação para este cliente/período
//             if (randreal() < pMutation) {
//                 int q_available = ind.deliveries[t_source][c];
//                 if (q_available <= 0) continue; // Pula se não houver entrega para mover

//                 // Sorteia um período de destino anterior e uma quantidade
//                 int t_dest = randint(0, t_source - 1);
//                 int q_move = randint(1, q_available);

//                 // Cria um indivíduo temporário para testar a factibilidade
//                 Individual temp_ind = ind;
//                 temp_ind.deliveries[t_source][c] -= q_move;
//                 temp_ind.deliveries[t_dest][c] += q_move;

//                 // Se a solução modificada for factível, atualiza o indivíduo original
//                 // A função check_feasibility já faz a verificação completa (estoque e frota)
//                 if (check_feasibility(temp_ind, irp)) {
//                     ind = temp_ind;
//                 }
//             }
//         }
//     }
// }




void advance_portion_mutation(Individual& ind, const IRP& irp, double pMutation) {
    // Percorre todos os clientes
    for (int c = 0; c < irp.nCustomers; ++c) {
        // Percorre os períodos a partir do segundo
        for (int t_source = 1; t_source < irp.nPeriods; ++t_source) {
            // Verifica a probabilidade de aplicar a mutação para este cliente/período
            if (randreal() < pMutation) {
                int q_available = ind.deliveries[t_source][c];
                if (q_available <= 0) continue; // Pula se não houver entrega para mover

                // Sorteia um período de destino anterior e uma quantidade
                int t_dest = randint(0, t_source - 1);
                int q_move = randint(1, q_available);

                // Cria um indivíduo temporário para testar a factibilidade
                Individual temp_ind = ind;
                temp_ind.deliveries[t_source][c] -= q_move;
                temp_ind.deliveries[t_dest][c] += q_move;

                // Usa a função de checagem robusta que verifica o "pico de estoque"
                if (check_feasibility(temp_ind, irp)) {
                    ind = temp_ind; // Se for factível, aplica a mudança
                }
                // Se não for factível, a mudança é simplesmente descartada
            }
        }
    }
}

void exportAndPlotRoutes(const IRP& irp, const Individual& best, const ACO_Params& acoParams,
                         const std::string& dataFilename, const std::string& pyScript) {
    std::ofstream out(dataFilename);
    if (!out.is_open()) return;

    int T = irp.nPeriods, nCust = irp.nCustomers, nVeh = irp.nVehicles;
    
    // <-- MUDANÇA AQUI: Adicionamos irp.Capacity ao cabeçalho do arquivo de dados -->
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

        ACO_Result res;
        long sumD = std::accumulate(del.begin(), del.end(), 0L);
        if (sumD > 0) res = runACO_for_period(irp, del, acoParams, false);

        out << "ROUTES\n";
        if (!res.bestRoutes.empty()) {
            for (const auto &route : res.bestRoutes) {
                for (size_t i = 0; i < route.size(); ++i) {
                    out << route[i] << (i + 1 < route.size() ? " " : "");
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