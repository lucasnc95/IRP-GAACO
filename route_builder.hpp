/*
 * NOVO ARQUIVO: route_builder.hpp
 *
 * Declara um construtor de rotas rápido (Algoritmo Sweep)
 * para substituir o ACO, que é muito lento.
 */

#ifndef ROUTE_BUILDER_HPP
#define ROUTE_BUILDER_HPP

#include "irp.hpp"
#include "individual.hpp"
#include "parameters.hpp" // Para ACO_Params (para a flag de busca local)
#include "route.hpp"
#include <vector>

/**
 * @brief Constrói um conjunto de rotas para um dia usando o algoritmo Sweep (Varredura Polar).
 *
 * 1. Calcula o ângulo de cada cliente que precisa de entrega.
 * 2. Ordena os clientes por ângulo.
 * 3. "Varre" os clientes, adicionando-os a uma rota até que a capacidade do
 * veículo seja atingida.
 * 4. Inicia uma nova rota e continua a varredura.
 * 5. Se o número de rotas geradas exceder a frota (irp.nVehicles), a
 * construção falha e retorna um vetor vazio.
 * 6. (Opcional) Aplica a busca local ('improve_routes') às rotas geradas.
 *
 * @param irp O problema IRP.
 * @param deliveries_for_period O plano de entregas (genótipo) para este dia.
 * @param aco_params Usado apenas para obter a flag 'pLocalSearch'.
 * @return Um vetor de 'Route' (fenótipo). Retorna vazio se for infactível.
 */
std::vector<Route> build_routes_with_sweep(
    const IRP& irp, 
    const std::vector<int>& deliveries_for_period,
    const ACO_Params& aco_params // Usado para a flag pLocalSearch
);


#endif // ROUTE_BUILDER_HPP