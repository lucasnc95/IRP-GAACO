/*
 * ARQUIVO MODIFICADO: individual.hpp
 * Atualiza a struct Individual para usar a nova struct Route.
 */
#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

#include <vector>
#include <map> // Necessário para std::map
#include "route.hpp" // <-- MUDANÇA: Inclui a nova struct de Rota

struct Individual {
    // A genética: O plano de entregas
    std::vector<std::vector<int>> deliveries;

    // --- A solução decodificada e seus custos (preenchidos pela avaliação) ---
    
    // <-- MUDANÇA: A representação das rotas foi alterada
    std::vector<std::vector<Route>> routes_per_period; 
    
    double routing_cost;
    double customer_holding_cost;
    double depot_holding_cost;
    double final_inventory_penalty;
    double fitness;
    bool is_feasible;

    // --- NOVAS MATRIZES AUXILIARES ---
    // Matriz de entrega mínima para garantir a demanda (preenchida pelo reparo/avaliação)
    std::vector<std::vector<int>> min_delivery_matrix;
    
    // Matriz de inventário final (após entrega e demanda)
    std::vector<std::vector<long>> final_inventory_matrix;


    // Construtor para inicializar com o tamanho correto e valores padrão
    Individual(int nPeriods = 0, int nCustomers = 0) {
        deliveries.assign(nPeriods, std::vector<int>(nCustomers, 0));
        
        // <-- MUDANÇA: Redimensiona o vetor externo de rotas
        routes_per_period.resize(nPeriods); 

        // <-- MUDANÇA: Inicializa as novas matrizes com o tamanho correto
        min_delivery_matrix.assign(nPeriods, std::vector<int>(nCustomers, 0));
        final_inventory_matrix.assign(nPeriods, std::vector<long>(nCustomers, 0));

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