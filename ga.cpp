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
#include <set>

using std::vector;


bool check_aco_consistency(const vector<int>& deliveries, const vector<Route>& routes) {
    std::set<int> customers_requiring_delivery;
    for (int c = 0; c < deliveries.size(); ++c) {
        if (deliveries[c] > 0) {
            customers_requiring_delivery.insert(c + 1); // 1-based
        }
    }

    std::set<int> customers_in_routes;
    for (const Route& route : routes) {
        for (int node : route.visits) {
            if (node > 0) { // Ignora o depósito (nó 0)
                customers_in_routes.insert(node);
            }
        }
    }
    return customers_requiring_delivery == customers_in_routes;
}


Individual make_new_heuristic_individual(const IRP& irp, const ACO_Params& aco_params) {
    Individual ind(irp.nPeriods, irp.nCustomers);
    vector<long> current_inv(irp.nCustomers);
    for(int i = 0; i < irp.nCustomers; ++i) current_inv[i] = irp.customers[i].initialInv;

    long fleet_capacity = (long)irp.nVehicles * irp.Capacity;
    const double LARGE_COST = 1e18;
    
    // Simulação dia a dia
    for (int t = 0; t < irp.nPeriods; ++t) {
        vector<int> min_deliveries_today(irp.nCustomers, 0);
        long period_load_min = 0;

        // --- PASSO 1: CALCULAR ENTREGAS MÍNIMAS OBRIGATÓRIAS ---
        for (int c = 0; c < irp.nCustomers; ++c) {
            long inv_after_demand = current_inv[c] - irp.customers[c].demand[t];
            
            if (inv_after_demand < irp.customers[c].minLevelInv) {
                // Cálculo da quantidade mínima exata necessária
                long needed = irp.customers[c].minLevelInv - inv_after_demand;
                
                // Verifica o espaço disponível no cliente (não pode exceder U_i)
                long space_available = irp.customers[c].maxLevelInv - current_inv[c];
                
                // A entrega é o mínimo entre o necessário, o espaço e a capacidade do veículo
                long delivery_amount = std::min({needed, space_available, (long)irp.Capacity});
                
                if (delivery_amount > 0) {
                    min_deliveries_today[c] = delivery_amount;
                    period_load_min += delivery_amount;
                }
            }
        }
        
        // Define o plano de entregas do dia como o plano mínimo
        ind.deliveries[t] = min_deliveries_today;
        long current_period_load = period_load_min;

        // --- PASSO 2: ENTREGAS OPORTUNISTAS (30% de chance) ---
        long remaining_fleet_cap = fleet_capacity - current_period_load;
        
        vector<int> customer_order(irp.nCustomers);
        std::iota(customer_order.begin(), customer_order.end(), 0);
        std::shuffle(customer_order.begin(), customer_order.end(), rng);

        for (int c : customer_order) {
            if (remaining_fleet_cap <= 0) break; // Frota do dia está cheia

            // Tenta adicionar com 30% de chance
            if (randreal() < 0.30) {
                // Calcula a demanda de N períodos à frente (N=1, 2 ou 3)
                int N = randint(1, 3);
                long q_extra = 0;
                for (int t_future = t + 1; t_future <= t + N && t_future < irp.nPeriods; ++t_future) {
                    q_extra += irp.customers[c].demand[t_future];
                }
                
                if (q_extra == 0) continue; // Sem demanda futura para adiantar

                // --- Checagem Tripla de Viabilidade ---
                
                long current_delivery = ind.deliveries[t][c];
                
                // 1. O cliente tem espaço para o extra? (Capacidade Cliente)
                long space = irp.customers[c].maxLevelInv - (current_inv[c] + current_delivery);
                
                // 2. A entrega total (mínima + extra) cabe em UM veículo? (Capacidade Veículo)
                long max_q_vehicle = (long)irp.Capacity - current_delivery;

                // 3. A frota do dia tem capacidade para o extra? (Capacidade Frota)
                // (remaining_fleet_cap já considera a carga mínima)

                // Limita o extra por todas as restrições
                q_extra = std::min({q_extra, space, max_q_vehicle, remaining_fleet_cap});

                // Adiciona a entrega oportunista se for válida
                if (q_extra > 0) {
                    ind.deliveries[t][c] += q_extra;
                    remaining_fleet_cap -= q_extra;
                    current_period_load += q_extra; // Atualiza a carga total do dia
                }
            }
        } // Fim do Passo 2

        // --- PASSO 3: ROTEAMENTO (COM ACO) ---
        ACO_Result ares;
        bool aco_success = false;
        
        // Tenta roteirizar o plano completo (Mínimo + Oportunista)
        if (current_period_load > 0) {
            ares = runACO_for_period(irp, ind.deliveries[t], aco_params, false);
            
            if (ares.bestCost < LARGE_COST / 10.0) {
                // Checa se o ACO roteirizou todos os clientes que precisavam
                aco_success = check_aco_consistency(ind.deliveries[t], ares.bestRoutes);
            }
        } else {
            // Sem entregas hoje, sucesso trivial
            aco_success = true;
            ares.bestRoutes.clear();
        }

        // --- PASSO 4: REPARO (se o ACO falhou) ---
        if (!aco_success && period_load_min > 0) {
            // "Discarte o inventário [oportunista]... e refaça"
            ind.deliveries[t] = min_deliveries_today; // Volta para o plano mínimo
            
            ares = runACO_for_period(irp, ind.deliveries[t], aco_params, false);
            
            if (ares.bestCost < LARGE_COST / 10.0) {
                aco_success = check_aco_consistency(ind.deliveries[t], ares.bestRoutes);
            } else {
                aco_success = false;
            }
        }
        
        // Armazena as rotas se o ACO foi bem-sucedido
        if (aco_success) {
            ind.routes_per_period[t] = ares.bestRoutes;
        } else {
            // A solução será marcada como infactível no final,
            // pois existe uma entrega obrigatória que não pôde ser roteirizada.
            ind.routes_per_period[t].clear(); 
            // Se min_deliveries_today > 0 e aco_success é false, a checagem
            // de factibilidade final irá falhar.
        }

        // --- PASSO 5: ATUALIZA O INVENTÁRIO (para o próximo período t+1) ---
        // Atualiza o inventário com base no plano de entregas FINAL (seja ele o
        // completo ou o reparado/mínimo)
        for (int c = 0; c < irp.nCustomers; ++c) {
            current_inv[c] += (long)ind.deliveries[t][c] - irp.customers[c].demand[t];
        }
    } // Fim do loop 't'
    
    // --- PASSO 6: CÁLCULO DE CUSTOS E FACTIBILIDADE ---
    // A função 'check_feasibility' fará a simulação final e
    // verificará se o 'current_inv' nunca violou os limites.
    if (check_feasibility(ind, irp)) {
        calculate_solution_costs(ind, irp); // Calcula os custos finais
    } else {
        // Se a simulação do check_feasibility falhar (o que não deveria
        // acontecer se a lógica aqui estiver correta), marque como infactível.
        ind.is_feasible = false;
        ind.fitness = LARGE_COST;
    }
    
    return ind;
}

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



Individual run_genetic_algorithm(const IRP& irp, const GA_Params& ga_params, const ACO_Params& aco_params, bool verbose) {
    vector<Individual> pop;
    pop.reserve(ga_params.popSize);
    std::cout << "criando população inicial..." << std::endl;
    
    while(pop.size() < ga_params.popSize) {
        Individual ind = make_new_heuristic_individual(irp, aco_params);
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
            Individual ind = make_new_heuristic_individual(irp, aco_params);
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
    
    // Adiciona Capacidade ao cabeçalho
    out << T << " " << nCust << " " << nVeh << " " << irp.Capacity << "\n";
    
    out << "DEPOT " << irp.depots[0].x << " " << irp.depots[0].y << "\n";

    vector<int> inv(nCust);
    for (int c = 0; c < nCust; ++c) inv[c] = irp.customers[c].initialInv;

    for (int t = 0; t < T; ++t) {
        out << "PERIOD " << t << "\n";
        out << "CUSTOMERS\n";
        
        // <-- MUDANÇA: Pega o plano de entregas do indivíduo 'best' -->
        const vector<int>& del = best.deliveries[t];
        
        for (int c = 0; c < nCust; ++c) {
            long inv_after = (long)inv[c] + (long)del[c] - (long)irp.customers[c].demand[t];
            out << (c+1) << " " << irp.customers[c].x << " " << irp.customers[c].y << " "
                << inv[c] << " " << del[c] << " " << irp.customers[c].demand[t] << " " << inv_after << "\n";
        }

        // As rotas já estão armazenadas no indivíduo 'best'
        out << "ROUTES\n";
        // <-- MUDANÇA: Itera sobre a nova struct de Rota -->
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