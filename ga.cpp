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



static vector<int> repair_period_row(const vector<int>& candidate,
                                     const vector<int>& inv_before,
                                     vector<long>& net_demand_to_satisfy, // Passado por referência
                                     int t,
                                     const IRP& irp) {
    int nCust = irp.nCustomers;
    long totalCap = (long)irp.nVehicles * (long)irp.Capacity;
    vector<int> deliver(nCust, 0);

    long min_needed_total = 0;

    // Passo 1: Calcular e alocar as entregas mínimas obrigatórias
    for (int c = 0; c < nCust; ++c) {
        long inv_after_demand = (long)inv_before[c] - irp.customers[c].demand[t];
        if (inv_after_demand < irp.customers[c].minLevelInv) {
            long needed = irp.customers[c].minLevelInv - inv_after_demand;
            long delivery_amount = std::min(needed, net_demand_to_satisfy[c]);
            deliver[c] = delivery_amount;
            min_needed_total += delivery_amount;
        }
    }
    
    long remainingCap = totalCap - min_needed_total;
    if (remainingCap <= 0) {
        for(int c = 0; c < nCust; ++c) net_demand_to_satisfy[c] -= deliver[c];
        return deliver;
    }
    
    // Passo 2: Distribuir a capacidade restante de forma não determinística
    struct CustInfo { int id; long space; };
    vector<CustInfo> customers_with_space;

    for (int c = 0; c < nCust; ++c) {
        long current_inv_after_min_del = (long)inv_before[c] + deliver[c];
        long space_available = irp.customers[c].maxLevelInv - current_inv_after_min_del;
        if (space_available > 0 && net_demand_to_satisfy[c] > deliver[c]) {
            customers_with_space.push_back({c, space_available});
        }
    }
    
    std::shuffle(customers_with_space.begin(), customers_with_space.end(), rng);

    for(const auto& cust_info : customers_with_space) {
        if (remainingCap <= 0) break;
        int c = cust_info.id;
        long space = cust_info.space;
        
        long max_extra_can_deliver = std::min({
            space, 
            (long)irp.Capacity - deliver[c], 
            remainingCap,
            net_demand_to_satisfy[c] - deliver[c]
        });
        
        if (max_extra_can_deliver > 0) {
            long random_delivery = randint(0, max_extra_can_deliver);
            deliver[c] += random_delivery;
            remainingCap -= random_delivery;
        }
    }
    
    for(int c = 0; c < nCust; ++c) {
        net_demand_to_satisfy[c] -= deliver[c];
    }

    return deliver;
}

void repair_individual(Individual& ind, const IRP& irp) {
    vector<long> net_demand_to_satisfy(irp.nCustomers, 0);
    for (int c = 0; c < irp.nCustomers; ++c) {
        long total_demand = 0;
        for (int t = 0; t < irp.nPeriods; ++t) {
            total_demand += irp.customers[c].demand[t];
        }
        net_demand_to_satisfy[c] = std::max(0L, total_demand - irp.customers[c].initialInv);
    }
    
    vector<int> current_inv(irp.nCustomers);
    for(int i = 0; i < irp.nCustomers; ++i) current_inv[i] = irp.customers[i].initialInv;

    for (int t = 0; t < irp.nPeriods; ++t) {
        ind.deliveries[t] = repair_period_row(ind.deliveries[t], current_inv, net_demand_to_satisfy, t, irp);
        for (int c = 0; c < irp.nCustomers; ++c) {
            current_inv[c] += ind.deliveries[t][c] - irp.customers[c].demand[t];
        }
    }
}

Individual make_simple_random_individual(const IRP& irp) {
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



// NOVA E ÚNICA MUTAÇÃO
void simple_random_mutation(Individual& ind, const IRP& irp) {
    int t = randint(0, irp.nPeriods - 1);
    int c = randint(0, irp.nCustomers - 1);
    
    long current_inv = irp.customers[c].initialInv;
    for(int p = 0; p < t; ++p) {
        current_inv += (long)ind.deliveries[p][c] - irp.customers[c].demand[p];
    }

    long space = irp.customers[c].maxLevelInv - current_inv;
    if (space > 0) {
        ind.deliveries[t][c] = randint(0, std::min({space, (long)irp.Capacity}));
    }
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
    
    int stagnation_counter = 0;
    double last_best_fitness = bestOverall.fitness;
    
    std::cout << "\nIniciando loop do GA para " << ga_params.nGen << " gerações..." << std::endl;
    for (int gen = 0; gen < ga_params.nGen; ++gen) {
        vector<Individual> newPop;
        newPop.reserve(ga_params.popSize);

        while (newPop.size() < ga_params.popSize) {
            Individual parent1 = tournamentSelect(pop, ga_params.tournamentK);
            Individual parent2 = tournamentSelect(pop, ga_params.tournamentK);
            
            std::pair<Individual, Individual> children;
            if (randreal() < 0.5) {
                children = one_point_crossover_customer(parent1, parent2, irp);
            } else {
                children = two_point_crossover_customer(parent1, parent2, irp);
            }

            // Processa o primeiro filho
            if (randreal() < ga_params.pMutation) {
                simple_random_mutation(children.first, irp);
            }
            if (check_feasibility(children.first, irp)) {
                evaluate_and_fill(children.first, irp, aco_params);
                newPop.push_back(children.first);
            }

            // Processa o segundo filho, se ainda houver espaço
            if (newPop.size() < ga_params.popSize) {
                if (randreal() < ga_params.pMutation) {
                    advance_portion_mutation(children.second, irp);
                }
                if (check_feasibility(children.second, irp)) {
                    evaluate_and_fill(children.second, irp, aco_params);
                    newPop.push_back(children.second);
                }
            }
        }

        // Elitismo
        Individual& bestPrev = *std::min_element(pop.begin(), pop.end(), 
            [](const Individual& a, const Individual& b){ return a.fitness < b.fitness; });
        auto worst_it = std::max_element(newPop.begin(), newPop.end(),
            [](const Individual& a, const Individual& b){ return a.fitness < b.fitness; });
        if(bestPrev.fitness < worst_it->fitness) *worst_it = bestPrev;
        pop.swap(newPop);

        Individual& genBest = *std::min_element(pop.begin(), pop.end(), 
            [](const Individual& a, const Individual& b){ return a.fitness < b.fitness; });
        if (genBest.fitness < bestOverall.fitness) bestOverall = genBest;
        
        // Estagnação e Reinicialização
        if (bestOverall.fitness < last_best_fitness - 1e-9) {
            last_best_fitness = bestOverall.fitness;
            stagnation_counter = 0;
        } else {
            stagnation_counter++;
        }
        if (ga_params.stagnation_threshold > 0 && stagnation_counter >= ga_params.stagnation_threshold) {
             if (ga_params.stagnation_threshold > 0 && stagnation_counter >= ga_params.stagnation_threshold) {
            std::cout << "\n*** ESTAGNAÇÃO DETECTADA APÓS " << stagnation_counter << " GERAÇÕES. REINICIALIZANDO POPULAÇÃO... ***\n" << std::endl;
            
            vector<Individual> restarted_pop;
            restarted_pop.reserve(ga_params.popSize);
            
            // 1. Mantém o melhor indivíduo de todos
            restarted_pop.push_back(bestOverall);
            
            // 2. Preenche o resto com novos indivíduos aleatórios e factíveis
            while(restarted_pop.size() < ga_params.popSize) {
                Individual ind = make_simple_random_individual(irp);
                if (check_feasibility(ind, irp)) {
                    restarted_pop.push_back(ind);
                }
            }

            // 3. Avalia apenas os novos indivíduos (o melhor já está avaliado)
            for(size_t i = 1; i < restarted_pop.size(); ++i) {
                 evaluate_and_fill(restarted_pop[i], irp, aco_params);
            }
            
            pop = restarted_pop; // Substitui a população estagnada
            stagnation_counter = 0; // Reseta o contador
        }
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


void advance_portion_mutation(Individual& ind, const IRP& irp) {
    int c = randint(0, irp.nCustomers - 1);

    vector<int> possible_source_periods;
    for (int t = 1; t < irp.nPeriods; ++t) {
        if (ind.deliveries[t][c] > 0) {
            possible_source_periods.push_back(t);
        }
    }
    if (possible_source_periods.empty()) return;

    int t_source = possible_source_periods[randint(0, possible_source_periods.size() - 1)];
    
    if (t_source == 0) return;
    int t_dest = randint(0, t_source - 1);

    int q_available = ind.deliveries[t_source][c];
    if (q_available <= 0) return;
    int q_move = randint(1, q_available);

    long current_dest_load = 0;
    for (int cust = 0; cust < irp.nCustomers; ++cust) {
        current_dest_load += ind.deliveries[t_dest][cust];
    }
    if (current_dest_load + q_move > (long)irp.nVehicles * irp.Capacity) {
        return;
    }

    long inv_at_start_of_dest = irp.customers[c].initialInv;
    for (int p = 0; p < t_dest; ++p) {
        inv_at_start_of_dest += (long)ind.deliveries[p][c] - irp.customers[c].demand[p];
    }
    if (inv_at_start_of_dest + ind.deliveries[t_dest][c] + q_move > irp.customers[c].maxLevelInv) {
        return;
    }

    ind.deliveries[t_source][c] -= q_move;
    ind.deliveries[t_dest][c] += q_move;
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