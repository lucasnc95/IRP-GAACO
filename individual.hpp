/*
 * ARQUIVO MODIFICADO: individual.hpp
 *
 * Adicionados campos para suportar o mecanismo de Biased Fitness (HGS).
 */
#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

#include <vector>
#include "route.hpp"

struct Individual {
    std::vector<std::vector<int>> deliveries;
    std::vector<std::vector<Route>> routes_per_period;
    
    double routing_cost;
    double customer_holding_cost;
    double depot_holding_cost;
    double final_inventory_penalty;
    double fitness; // Custo total (penalizado se infactível)
    bool is_feasible;

    // --- NOVOS CAMPOS PARA SELEÇÃO DE SOBREVIVENTES (HGS) ---
    double diversity_contribution; // Distância média para os vizinhos mais próximos
    int rank_cost;                 // Classificação baseada apenas no custo (0 = melhor)
    int rank_diversity;            // Classificação baseada na diversidade (0 = mais diverso)
    double biased_fitness;         // Fitness combinado (Custo + Diversidade)

    // Construtor
    Individual(int nPeriods = 0, int nCustomers = 0) {
        deliveries.assign(nPeriods, std::vector<int>(nCustomers, 0));
        routes_per_period.resize(nPeriods);
        
        routing_cost = 0.0;
        customer_holding_cost = 0.0;
        depot_holding_cost = 0.0;
        final_inventory_penalty = 0.0;
        fitness = 1e18;
        is_feasible = false;

        // Inicializa novos campos
        diversity_contribution = 0.0;
        rank_cost = 0;
        rank_diversity = 0;
        biased_fitness = 0.0;
    }
};

#endif // INDIVIDUAL_HPP