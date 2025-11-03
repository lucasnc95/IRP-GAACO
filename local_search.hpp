/*
 * NOVO ARQUIVO: local_search.hpp
 * Declara as funções de busca local que operam na nova struct Route.
 */
#ifndef LOCAL_SEARCH_HPP
#define LOCAL_SEARCH_HPP

#include "route.hpp" // Inclui a nova definição de Rota
#include "irp.hpp"
#include <vector>

// Assinaturas das funções atualizadas em local_search.cpp
// Note que elas agora recebem um vetor de `Route`

/**
 * @brief Calcula o custo total (distância) de um conjunto de rotas.
 */
static double calculate_total_cost(const std::vector<Route>& routes);

/**
 * @brief Aplica o 2-Opt (Best Improvement) em uma única rota.
 * Atualiza o custo da rota se uma melhoria for encontrada.
 */
bool apply_2opt_on_route(Route& route, const std::vector<std::vector<double>>& dist);

/**
 * @brief Aplica o Swap (troca) inter-rotas, usando a estratégia Best Improvement.
 * Atualiza custos e capacidades das rotas envolvidas.
 */
bool apply_inter_route_swap(std::vector<Route>& routes, const IRP& irp, 
                            const std::vector<int>& deliveries, const std::vector<std::vector<double>>& dist);

/**
 * @brief Aplica o Relocate (realocação) inter-rotas, usando a estratégia Best Improvement.
 * Atualiza custos e capacidades das rotas envolvidas.
 */
bool apply_inter_route_relocate(std::vector<Route>& routes, const IRP& irp, 
                                const std::vector<int>& deliveries, const std::vector<std::vector<double>>& dist);

/**
 * @brief Orquestrador da busca local, aplica operadores até atingir um ótimo local.
 */
double improve_routes(std::vector<Route>& routes, const IRP& irp, 
                      const std::vector<int>& deliveries_for_period, 
                      const std::vector<std::vector<double>>& dist, bool verbose);

/**
 * @brief Função auxiliar para recalcular o custo e a carga de uma rota do zero.
 */
void recalculate_route_metrics(Route& route, const IRP& irp, 
                               const std::vector<int>& deliveries, 
                               const std::vector<std::vector<double>>& dist);

#endif // LOCAL_SEARCH_HPP