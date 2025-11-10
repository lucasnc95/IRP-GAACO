/*
 * ARQUIVO MODIFICADO: evaluation.hpp
 *
 * A interface foi limpa conforme solicitado.
 * - 'evaluate_and_fill' foi removida.
 * - Adicionadas 'check_feasibility' e 'calculate_total_cost' com responsabilidades únicas.
 */

#ifndef EVALUATION_HPP
#define EVALUATION_HPP

#include "irp.hpp"
#include "individual.hpp"
#include "parameters.hpp" // Para ACO_Params (se check_feasibility precisar dele)
#include "route.hpp"

/**
 * @brief (NOVA FUNÇÃO) Verifica rigorosamente a factibilidade completa de um indivíduo.
 *
 * Esta função assume que 'ind.deliveries' e 'ind.routes_per_period'
 * já estão preenchidos. Ela valida consistência, capacidade de rotas/frota
 * e simula o inventário para checar stock-outs ou excessos.
 *
 * @param ind O indivíduo (com entregas e rotas).
 * @param irp O problema.
 * @return true se o indivíduo é 100% factível, false caso contrário.
 */
void check_feasibility(Individual& ind, const IRP& irp);


/**
 * @brief (NOVA FUNÇÃO) Calcula os custos de um indivíduo factível.
 *
 * Esta função assume que 'check_feasibility' já retornou 'true'.
 * Ela preenche os campos:
 * - ind.routing_cost
 * - ind.customer_holding_cost
 * - ind.depot_holding_cost
 * - ind.fitness
 *
 * @param ind O indivíduo (passado por referência) para preencher os custos.
 * @param irp O problema.
 */
void calculate_total_cost(Individual& ind, const IRP& irp);


#endif // EVALUATION_HPP