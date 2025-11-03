/*
 * NOVO ARQUIVO: route.hpp
 * Define a estrutura de dados central para uma única rota de veículo.
 */
#ifndef ROUTE_HPP
#define ROUTE_HPP

#include <vector>

/**
 * @struct Route
 * @brief Representa uma única rota de veículo em um dia.
 *
 * Contém a sequência de visitas, a capacidade restante e o custo total da rota.
 */
struct Route {
    // Vetor de inteiros indicando a rota (ex: 0 -> 3 -> 5 -> 0)
    std::vector<int> visits; 
    
    // A capacidade restante do veículo
    int remaining_capacity; 
    
    // O custo total (distância) desta rota
    double cost; 

    // Construtor padrão
    Route() : remaining_capacity(0), cost(0.0) {
        // Inicializa a rota vazia, mas pronta para ser preenchida
        visits.clear();
    }
};

#endif // ROUTE_HPP