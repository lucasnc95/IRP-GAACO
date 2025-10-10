#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

#include <vector>

struct Individual {
    // A genética: O plano de entregas
    std::vector<std::vector<int>> deliveries;

    // --- A solução decodificada e seus custos (preenchidos pela avaliação) ---
    std::vector<std::vector<std::vector<int>>> routes_per_period;
    double routing_cost;
    double customer_holding_cost;
    double depot_holding_cost;
    double final_inventory_penalty;
    double fitness;
    bool is_feasible;

    // Construtor para inicializar com o tamanho correto e valores padrão
    Individual(int nPeriods = 0, int nCustomers = 0) {
        deliveries.assign(nPeriods, std::vector<int>(nCustomers, 0));
        routes_per_period.resize(nPeriods);
        // Inicializa custos com -1 para indicar que não foi avaliado
        routing_cost = -1.0;
        customer_holding_cost = -1.0;
        depot_holding_cost = -1.0;
        final_inventory_penalty = -1.0;
        fitness = -1.0;
        is_feasible = false;
    }
};

#endif // INDIVIDUAL_HPP