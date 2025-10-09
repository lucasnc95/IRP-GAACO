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

double calculate_population_variability(const std::vector<Individual>& pop) {
    if (pop.size() < 2) {
        return 0.0;
    }

    double total_distance = 0.0;
    int num_pairs = 0;
    
    for (size_t i = 0; i < pop.size(); ++i) {
        for (size_t j = i + 1; j < pop.size(); ++j) {
            double distance = 0.0;
            // Calcula a soma das diferenças absolutas (Distância de Manhattan) entre duas matrizes
            for (size_t t = 0; t < pop[i].deliveries.size(); ++t) {
                for (size_t c = 0; c < pop[i].deliveries[t].size(); ++c) {
                    distance += std::abs(pop[i].deliveries[t][c] - pop[j].deliveries[t][c]);
                }
            }
            total_distance += distance;
            num_pairs++;
        }
    }
    return total_distance / num_pairs;
}

double evaluateIndividual(
    const Individual& ind, 
    const IRP& irp, 
    const ACO_Params& aco_params, 
    bool verbose
) {
    EvaluationResult result = simulate_and_evaluate(ind, irp, aco_params);
    double final_fitness = result.operational_cost();
    
    if(verbose) {
        std::cout << "  - Fitness final avaliado: " << final_fitness << std::endl;
    }

    return final_fitness;
}


Individual crossover(const Individual& a, const Individual& b, const IRP& irp) {
    Individual child(irp.nPeriods, irp.nCustomers); 
    int crosspoint = randint(0, irp.nCustomers - 1);
    for (int t = 0; t < irp.nPeriods; ++t) {
        for (int c = 0; c < irp.nCustomers; ++c) {
            child.deliveries[t][c] = (c < crosspoint) ? a.deliveries[t][c] : b.deliveries[t][c];
        }
    }
    // As rotas e custos do filho estão vazios/inválidos por padrão
    return child;
}

void raw_mutate(Individual& ind, const IRP& irp, double pMut) {
    for (int t = 0; t < irp.nPeriods; ++t) {
        for (int c = 0; c < irp.nCustomers; ++c) {
            if (randreal() < pMut) {
                int max_q = irp.customers[c].maxLevelInv;
                int current_q = ind.deliveries[t][c];
                int change = randint(-current_q, max_q / 4);
                int new_q = std::max(0, current_q + change);
                ind.deliveries[t][c] = new_q;
            }
        }
    }
}

void advance_portion_mutation(Individual& ind, const IRP& irp) {
    // 1. Escolhe um cliente aleatório
    int c = randint(0, irp.nCustomers - 1);

    // 2. Encontra os períodos onde uma entrega ocorre (e pode ser adiantada)
    vector<int> possible_source_periods;
    for (int t = 1; t < irp.nPeriods; ++t) {
        if (ind.deliveries[t][c] > 0) {
            possible_source_periods.push_back(t);
        }
    }
    if (possible_source_periods.empty()) return; // Nada a fazer

    int t_source = possible_source_periods[randint(0, possible_source_periods.size() - 1)];
    
    // 3. Encontra um período de destino anterior
    if (t_source == 0) return; // Não pode adiantar do período 0
    int t_dest = randint(0, t_source - 1);

    // 4. Sorteia uma quantidade a ser movida
    int q_available = ind.deliveries[t_source][c];
    if (q_available <= 0) return;
    int q_move = randint(1, q_available);

    // 5. Verifica a factibilidade do movimento
    // 5.1 Checagem de capacidade da frota no período de destino
    long current_dest_load = 0;
    for (int cust = 0; cust < irp.nCustomers; ++cust) {
        current_dest_load += ind.deliveries[t_dest][cust];
    }
    if (current_dest_load + q_move > (long)irp.nVehicles * irp.Capacity) {
        return; // Falha na checagem de capacidade da frota
    }

    // 5.2 Checagem de capacidade do armazém do cliente no período de destino
    long inv_at_start_of_dest = irp.customers[c].initialInv;
    for (int p = 0; p < t_dest; ++p) {
        inv_at_start_of_dest += (long)ind.deliveries[p][c] - irp.customers[c].demand[p];
    }
    if (inv_at_start_of_dest + ind.deliveries[t_dest][c] + q_move > irp.customers[c].maxLevelInv) {
        return; // Falha na checagem de capacidade do armazém
    }

    // 6. Se factível, aplica o movimento
    ind.deliveries[t_source][c] -= q_move;
    ind.deliveries[t_dest][c] += q_move;
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

void delivery_delay_mutation(Individual& ind, const IRP& irp) {
    vector<int> customer_indices(irp.nCustomers);
    std::iota(customer_indices.begin(), customer_indices.end(), 0);
    std::shuffle(customer_indices.begin(), customer_indices.end(), rng);
    
    for (int c : customer_indices) {
        for (int t = 0; t < irp.nPeriods - 1; ++t) {
            int q_to_move = ind.deliveries[t][c];
            if (q_to_move == 0) continue;

            int existing_delivery_t_plus_1 = ind.deliveries[t+1][c];

            // Checa se a soma das entregas cabe em um único veículo
            if (q_to_move + existing_delivery_t_plus_1 > irp.Capacity) continue;

            // Simula a evolução do estoque para checar factibilidade
            bool is_feasible = true;
            long current_inv = (t > 0) ? (long)irp.customers[c].initialInv : irp.customers[c].initialInv;
            if (t > 0) {
                for(int p=0; p < t; ++p) current_inv += (long)ind.deliveries[p][c] - irp.customers[c].demand[p];
            }
            
            // Checa o período de origem (t)
            if (current_inv - irp.customers[c].demand[t] < irp.customers[c].minLevelInv) {
                is_feasible = false;
            }
            
            // Checa o período de destino (t+1)
            if(is_feasible) {
                long inv_at_start_t1 = current_inv - irp.customers[c].demand[t];
                if(inv_at_start_t1 + q_to_move + existing_delivery_t_plus_1 > irp.customers[c].maxLevelInv) {
                    is_feasible = false;
                }
            }

            if (is_feasible) {
                ind.deliveries[t][c] = 0;
                ind.deliveries[t+1][c] += q_to_move;
                return; // First improvement
            }
        }
    }
}




Individual make_stock_up_individual(const IRP& irp) {
    Individual ind(irp.nPeriods, irp.nCustomers);
    ind.deliveries.assign(irp.nPeriods, vector<int>(irp.nCustomers, 0));
    repair_individual(ind, irp); // Repara um indivíduo vazio, forçando uma solução factível
    return ind;
}

// Estratégia "Just-in-Time"
Individual make_just_in_time_individual(const IRP& irp) {
 Individual ind(irp.nPeriods, irp.nCustomers); 
    ind.deliveries.assign(irp.nPeriods, vector<int>(irp.nCustomers, 0));
    ind.fitness = 1e18;

    vector<long> current_inv(irp.nCustomers);
    for(int i = 0; i < irp.nCustomers; ++i) {
        current_inv[i] = irp.customers[i].initialInv;
    }

    long fleet_capacity = (long)irp.nVehicles * irp.Capacity;

    for (int t = 0; t < irp.nPeriods; ++t) {
        vector<int> deliveries_t(irp.nCustomers, 0);
        long current_period_load = 0;

        for (int c = 0; c < irp.nCustomers; ++c) {
            long inv_after_consumption = current_inv[c] - irp.customers[c].demand[t];
            
            if (inv_after_consumption < irp.customers[c].minLevelInv) {
                long needed = irp.customers[c].minLevelInv - inv_after_consumption;
                long space_available = irp.customers[c].maxLevelInv - current_inv[c];
                long delivery_amount = std::min(needed, space_available);

                if (delivery_amount > irp.Capacity) delivery_amount = irp.Capacity;

                if (current_period_load + delivery_amount <= fleet_capacity) {
                    deliveries_t[c] = delivery_amount;
                    current_period_load += delivery_amount;
                }
            }
        }
        
        ind.deliveries[t] = deliveries_t;
        for (int c = 0; c < irp.nCustomers; ++c) {
            current_inv[c] += (long)ind.deliveries[t][c] - irp.customers[c].demand[t];
        }
    }
    return ind;
}


Individual run_genetic_algorithm(const IRP& irp, const GA_Params& ga_params, const ACO_Params& aco_params, bool verbose) {
    vector<Individual> pop;
    pop.reserve(ga_params.popSize);
    
    // Geração da População Inicial
    int half_pop = ga_params.popSize / 2;
    std::cout << "Inicializando " << half_pop << " indivíduos com a estratégia 'Stock-Up'..." << std::endl;
    for (int i = 0; i < half_pop; ++i) {
        pop.push_back(make_stock_up_individual(irp));
    }
    std::cout << "Inicializando " << (ga_params.popSize - half_pop) << " indivíduos com a estratégia 'Just-in-Time'..." << std::endl;
    for (int i = 0; i < (ga_params.popSize - half_pop); ++i) {
        pop.push_back(make_just_in_time_individual(irp));
    }

    // Avaliação Inicial
    std::cout << "Avaliando população inicial..." << std::endl;
    for (size_t i = 0; i < pop.size(); ++i) {
        if (verbose) {
            std::cout << "[Gen 0] Avaliando Indivíduo " << i + 1 << "/" << pop.size() << "..." << std::endl;
        }
        EvaluationResult res = simulate_and_evaluate(pop[i], irp, aco_params);
        pop[i].fitness = res.total_fitness();
    }

    Individual bestOverall = *std::min_element(pop.begin(), pop.end(), 
        [](const Individual& a, const Individual& b){ return a.fitness < b.fitness; });
    
    int stagnation_counter = 0;
    double last_best_fitness = bestOverall.fitness;
    
    std::cout << "\nIniciando loop do GA para " << ga_params.nGen << " gerações..." << std::endl;
    for (int gen = 0; gen < ga_params.nGen; ++gen) {
        vector<Individual> newPop;
        newPop.reserve(ga_params.popSize);

        while ((int)newPop.size() < ga_params.popSize) {
            Individual parent1 = tournamentSelect(pop, ga_params.tournamentK);
            Individual parent2 = tournamentSelect(pop, ga_params.tournamentK);
            Individual child = crossover(parent1, parent2, irp);
            
            if (randreal() < ga_params.pMutation) {
                advance_portion_mutation(child, irp);
            }
            
            // A verificação de factibilidade da frota é feita antes da avaliação completa
            bool fleet_capacity_ok = true;
            long total_fleet_capacity = (long)irp.nVehicles * irp.Capacity;
            for (int t = 0; t < irp.nPeriods; ++t) {
                long period_load = std::accumulate(child.deliveries[t].begin(), child.deliveries[t].end(), 0L);
                if (period_load > total_fleet_capacity) {
                    fleet_capacity_ok = false;
                    break;
                }
            }

            if (fleet_capacity_ok) {
                if (verbose) {
                    std::cout << "[Gen " << gen + 1 << "] Avaliando Indivíduo " << newPop.size() + 1 << "/" << ga_params.popSize << "..." << std::endl;
                }
                EvaluationResult res = simulate_and_evaluate(child, irp, aco_params);
                child.fitness = res.total_fitness();
            } else {
                child.fitness = 1e18; // Pena de morte
            }
            
            newPop.push_back(child);
        }

        Individual& bestPrev = *std::min_element(pop.begin(), pop.end(), 
            [](const Individual& a, const Individual& b){ return a.fitness < b.fitness; });
        auto worst_it = std::max_element(newPop.begin(), newPop.end(),
            [](const Individual& a, const Individual& b){ return a.fitness < b.fitness; });
        if(bestPrev.fitness < worst_it->fitness){
            *worst_it = bestPrev;
        }
        pop.swap(newPop);

        Individual& genBest = *std::min_element(pop.begin(), pop.end(), 
            [](const Individual& a, const Individual& b){ return a.fitness < b.fitness; });
        if (genBest.fitness < bestOverall.fitness) {
            bestOverall = genBest;
        }
        
        // Lógica de Detecção de Estagnação e Reinicialização
        if (bestOverall.fitness < last_best_fitness - 1e-9) {
            last_best_fitness = bestOverall.fitness;
            stagnation_counter = 0;
        } else {
            stagnation_counter++;
        }

        if (ga_params.stagnation_threshold > 0 && stagnation_counter >= ga_params.stagnation_threshold) {
            std::cout << "\n*** ESTAGNAÇÃO DETECTADA APÓS " << stagnation_counter << " GERAÇÕES. REINICIALIZANDO POPULAÇÃO... ***\n" << std::endl;
            
            vector<Individual> restarted_pop;
            restarted_pop.reserve(ga_params.popSize);
            restarted_pop.push_back(bestOverall);
            
            int half_remaining = (ga_params.popSize - 1) / 2;
            for(int i = 0; i < half_remaining; ++i) restarted_pop.push_back(make_stock_up_individual(irp));
            for(int i = 0; i < (ga_params.popSize - 1 - half_remaining); ++i) restarted_pop.push_back(make_just_in_time_individual(irp));

            for(size_t i = 1; i < restarted_pop.size(); ++i) {
                 EvaluationResult res = simulate_and_evaluate(restarted_pop[i], irp, aco_params);
                 restarted_pop[i].fitness = res.total_fitness();
            }
            
            pop = restarted_pop;
            stagnation_counter = 0;
        }
        
        // Cálculo de Estatísticas e Log
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