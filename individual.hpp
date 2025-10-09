#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

#include <vector>

struct Individual {
    // A genética: O plano de entregas
    std::vector<std::vector<int>> deliveries;

    // --- A solução decodificada e seus custos ---
    // (Estes campos são preenchidos pela função de avaliação)
    
    // As rotas exatas para cada período
    std::vector<std::vector<std::vector<int>>> routes_per_period;
    
    // Custos detalhados
    double routing_cost = -1.0;
    double customer_holding_cost = -1.0;
    double depot_holding_cost = -1.0;
    double final_inventory_penalty = -1.0;

    // O fitness final, calculado a partir dos custos acima
    double fitness = -1.0;
    bool is_feasible = false;

    // Construtor para inicializar com o tamanho correto
    Individual(int nPeriods, int nCustomers) {
        deliveries.assign(nPeriods, std::vector<int>(nCustomers, 0));
        routes_per_period.resize(nPeriods);
    }
    // Construtor padrão
    Individual() = default;
};

#endif // INDIVIDUAL_HPP