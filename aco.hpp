/*
 * ARQUIVO MODIFICADO: aco.hpp
 * Atualiza a struct ACO_Result para usar a nova struct Route.
 */
#ifndef ACO_HPP
#define ACO_HPP

#include "irp.hpp"
#include "parameters.hpp"
#include "route.hpp" // <-- MUDANÇA: Inclui a nova struct de Rota
#include <vector>

/**
 * @struct ACO_Result
 * @brief Armazena o resultado da execução do ACO para um período.
 */
struct ACO_Result {
    double bestCost;
    
    // <-- MUDANÇA: Agora armazena um vetor de structs Route
    std::vector<Route> bestRoutes; 
};

/**
 * @brief Executa o algoritmo ACO para um único período.
 * @param irp O problema IRP (para distâncias e capacidades).
 * @param deliveries_for_period O plano de entregas para este período.
 * @param params Parâmetros do ACO.
 * @param verbose Flag para imprimir logs.
 * @return Um ACO_Result contendo as melhores rotas encontradas e seu custo.
 */
ACO_Result runACO_for_period(const IRP& irp,
                             const std::vector<int>& deliveries_for_period,
                             const ACO_Params& params,
                             bool verbose);

#endif // ACO_HPP